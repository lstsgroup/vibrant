!
!   Copyright 2025 Ekin E. Winogradow, Johannes Scheffler, Moritz Leucke, Dorothea Golze
!
!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at
!
!       http://www.apache.org/licenses/LICENSE-2.0
!
!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.
!

!> @brief Module containing all procedures involving the calculation of the final
!> spectrum requested by the user.
MODULE calc_spectra

    USE kinds, ONLY: dp, str_len
    USE output_io, ONLY: write_spectra_data, write_mol_file, append_column, check_file_open
    USE ISO_FORTRAN_ENV, ONLY: output_unit, error_unit
    USE vib_types, ONLY: global_settings, systems, molecular_dynamics, static, dipoles, raman

    USE constants, ONLY: pi, fs2s, debye, speed_light, const_planck, const_boltz, const_permit, cm2m, a3_to_debye_per_e, &
                         hartreebohr2evang, hessian_factor, bohr2ang, reccm2ev, am_u, debye2cm, avo_num, au2vm, ang, at_u, &
                         speed_light_au, debye, reccm2au, sinc_factor
    USE read_traj, ONLY: read_coord_frame!, check_file_open
    USE fin_diff, ONLY: central_diff, forward_diff
    USE vel_cor, ONLY: cvv, cvv_iso, cvv_aniso, cvv_only_x, cvv_resraman
    USE dipole_calc!, ONLY: compute_dipole, check_jumps
    USE pade, ONLY: interpolate

    USE, INTRINSIC                              :: ISO_C_BINDING
    USE OMP_LIB

    IMPLICIT NONE

    INCLUDE 'fftw3.f03'

    PUBLIC :: spec_power, normal_mode_analysis, spec_static_ir, spec_static_raman, spec_ir, spec_raman, spec_abs, spec_static_resraman, spec_abs_md, spec_resraman

CONTAINS

    !> @brief Computes the vibrational power spectrum from an MD trajectory.
    !>
    !> This subroutine evaluates the velocity autocorrelation function (VACF) and
    !> its Fourier transform to obtain the vibrational power spectrum (or power
    !> spectral density). It supports trajectories defined either by positions or
    !> velocities.
    !>
    !> @param[inout] gs  -- Global settings structure (FFT, constants, control parameters)
    !> @param[inout] sys -- System data structure (provides trajectory and atomic data)
    !> @param[inout] md  -- MD structure with coordinates, velocities, and correlation data.
    !>
    SUBROUTINE spec_power(gs, sys, md)

        TYPE(global_settings), INTENT(INOUT)        :: gs
        TYPE(systems), INTENT(INOUT)                :: sys
        TYPE(molecular_dynamics), INTENT(INOUT)     :: md

        INTEGER                                                  :: stat, i, runit
        INTEGER(kind=dp)                                          :: plan
        CHARACTER(LEN=str_len)                                          :: msg, spectra_file_name, c_label
        REAL(kind=dp)                               :: freq_range, freq_res, power_const
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE      :: power_int, freq

        !Allocate
        ALLOCATE (md%zhat(0:md%t_cor*2 - 1), power_int(0:md%t_cor*2 - 1), freq(0:md%t_cor*2 - 1))
        !Initialize
        md%zhat = COMPLEX(0._dp, 0.0_dp)

        CALL read_coord_frame(sys%natom, sys%filename, md%coord_v, sys)

        !!If it is from positions, do finite differences first
        IF (sys%type_traj=='pos') THEN
            CALL central_diff(sys%natom, md%coord_v, md%v, sys, md)

            CALL cvv(sys%natom, md%v, sys, gs, md)

        !!If it is from velocities, compute autocorrelation directly
        ELSEIF (sys%type_traj=='vel') THEN

            CALL cvv(sys%natom, md%coord_v, sys, gs, md)

        END IF

        !Determine the maximum frequency range in cm^{-1} based on `md%dt`
        freq_range = REAL((1.0_dp/(md%dt*fs2s))/speed_light, kind=dp)
        !Determine the frequency resolution
        freq_res = REAL(freq_range/(2.0_dp*md%t_cor), kind=dp)

        !!Call FFT
        CALL dfftw_plan_dft_r2c_1d(plan, 2*md%t_cor, md%z, md%zhat, FFTW_ESTIMATE) !!!FFT
        CALL dfftw_execute_dft_r2c(plan, md%z, md%zhat)
        CALL dfftw_destroy_plan(plan)

        !!Unit conversion to K.cm
        power_const = (md%dt*fs2s*am_u*speed_light/const_boltz)*2.0_dp/(sys%natom*3.0_dp)

        !Generate power spectrum
        DO i = 0, 2*md%t_cor - 1
            freq(i) = i*freq_res
            power_int(i) = REAL(md%zhat(i), KIND=dp)*power_const
            IF (freq(i).GE.5000_dp) CYCLE
        END DO

        !!Write the results to a file
        c_label = "INT"
        spectra_file_name = "power_spec.txt"
        CALL write_spectra_data(spectra_file_name, c_label, freq, power_int, 5000.0_dp)

        DEALLOCATE (power_int, freq)

    END SUBROUTINE spec_power
!***********************************************************************************************!
!***********************************************************************************************!

    !> @brief Computes the infrared (IR) absorption spectrum from dipole trajectories.
    !>
    !> This subroutine calculates the IR spectrum using time-correlation functions of
    !> dipole moment derivatives (dipole velocities) obtained from molecular dynamics
    !> simulations. It supports both total-system and fragment-based spectra, and
    !> two dipole types: **Wannier** (from localized orbitals) and **Berry-phase** (from
    !> polarization).
    !>
    !> @param[inout] gs   -- Global settings (constants, FFT parameters, etc.)
    !> @param[inout] sys  -- System data structure (provides atomic, cell, and fragment information)
    !> @param[inout] md   -- MD structure (coordinates, velocities, correlation arrays)
    !> @param[inout] dips -- Dipole trajectory data structure (dipole type, dipole arrays, filenames)
    !>
    SUBROUTINE spec_ir(gs, sys, md, dips)
        TYPE(global_settings), INTENT(INOUT) :: gs
        TYPE(systems), INTENT(INOUT) :: sys
        TYPE(molecular_dynamics), INTENT(INOUT) :: md
        TYPE(dipoles), INTENT(INOUT) :: dips

        CHARACTER(LEN=str_len) :: msg
        CHARACTER(LEN=40) :: output_fname, f_name, idx, spectra_file_name, c_label
        INTEGER :: stat, i, i_group, runit
        INTEGER(KIND=dp) :: plan
        REAL(KIND=dp) :: freq_range, freq_res, sinc_const, ir_const
        REAL(KIND=dp), ALLOCATABLE :: ir_int(:), freq(:), z_frag(:, :)
        REAL(KIND=dp), ALLOCATABLE :: v_frag(:, :, :, :), dipole_frag(:, :, :, :)
        COMPLEX(KIND=dp), ALLOCATABLE :: zhat_frag(:, :)

        !!Allocate and initialize
        ALLOCATE (md%zhat(0:2*md%t_cor - 1), ir_int(0:2*md%t_cor - 1), freq(0:2*md%t_cor - 1))
        md%zhat = CMPLX(0._dp, 0.0_dp)
        !!If there are no fragments
        IF (.NOT. sys%fragments%frag) THEN
            sys%fragments%ngroup = 1
        END IF

        !Determine the maximum frequency range in cm^{-1} based on `md%dt`
        freq_range = REAL((1._dp/(md%dt*fs2s))/speed_light, dp)
        !Determine the frequency resolution
        freq_res = REAL(freq_range/(2._dp*md%t_cor), dp)
        !Apply sinc function
        sinc_const = freq_res*md%dt*fs2s*2._dp*pi*speed_light
        !!final IR units in K*cm*km*mol^-1
        ir_const = avo_num*md%dt*fs2s*2._dp*10._dp/(12._dp*const_permit*speed_light*const_boltz)

        !!Start with looping over fragments if there are any
        DO i_group = 1, sys%fragments%ngroup
            ! Dipole handling
            IF (dips%type_dipole=='wannier') THEN
                CALL read_coord_frame(sys%natom, dips%dip_file, md%coord_v, sys)
                IF (.NOT. sys%fragments%frag) THEN
                    CALL compute_dipole(dips%dipole, sys, md)
                    ALLOCATE (md%v(sys%framecount, sys%natom, 3))
                    CALL central_diff(sys%mol_num, dips%dipole, md%v, sys, md)
                !!If there are fragments, assign their wannier centers first
                ELSE
                    CALL assign_wannier(sys, md)
                    CALL compute_dipole_frag(dipole_frag, sys, md)
                    sys%mol_num = sys%fragments%type_frag(i_group)%nfrag
                    CALL central_diff(sys%mol_num, dipole_frag(i_group, :, 1:sys%mol_num, :), md%v, sys, md)
                END IF
            ELSEIF (dips%type_dipole=='berry') THEN
                CALL read_coord_frame(sys%mol_num, dips%dip_file, md%coord_v, sys)
                !!Correct discontinuities due to polarization quantum, if there are any
                CALL check_jumps(md%coord_v, sys%mol_num, sys)
                ALLOCATE (md%v(sys%framecount, sys%natom, 3))
                CALL central_diff(sys%mol_num, md%coord_v, md%v, sys, md)
            END IF
            !!Call dipole autocorrelation function
            CALL cvv(sys%mol_num, md%v, sys, gs, md)

            !Call FFT
            CALL dfftw_plan_dft_r2c_1d(plan, 2*md%t_cor, md%z, md%zhat, FFTW_ESTIMATE)
            CALL dfftw_execute_dft_r2c(plan, md%z, md%zhat)
            CALL dfftw_destroy_plan(plan)

            !! Convert from debye to cm^{-1}
            md%zhat = REAL(md%zhat*debye2cm*debye2cm, dp)

            !! Generate the IR spectrum
            DO i = 0, 2*md%t_cor - 1
                freq(i) = i*freq_res
                ir_int(i) = md%zhat(i)*ir_const*(sinc_const*i/SIN(sinc_const*i))**2._dp
                IF (freq(i)>=5000._dp) CYCLE
                ir_int(0) = 0._dp
                ir_int(i) = -1.0_dp*REAL(ir_int(i))
            END DO

            !!Write the results to a file
            c_label = "INT"
            !!If there are no fragments
            IF (.NOT. sys%fragments%frag) THEN
                spectra_file_name = "IR_spectrum.txt"
            !!If there are fragments, generate separate output files for each of them
            ELSE
                WRITE (spectra_file_name, '(A,I0,A)') "IR_spectrum_fragment_", i_group, ".txt"
            END IF
            CALL write_spectra_data(spectra_file_name, c_label, freq, ir_int, 5000.0_dp)
            DEALLOCATE (md%z)
        END DO
        DEALLOCATE (ir_int, freq)
    END SUBROUTINE spec_ir

!****************************************************************************************!
!****************************************************************************************!
    !> @brief Computes Raman spectra from time-dependent polarizability derivatives.
    !>
    !> This subroutine calculates  Raman spectra from MD trajectories, using the autocorrelation
    !> functions of polarizability tensor derivatives obtained via finite electric field perturbations.
    !> It supports both **Wannier** and **Berry-phase** dipole formulations, direct DFPT
    !> polarizabilities, and fragment-based analysis. The resulting Raman intensities are
    !> computed for multiple incident laser frequencies and written to output files.
    !>
    !> @param[inout] gs   --  Global settings (temperature, constants, and FFT control)
    !> @param[inout] sys  --  System information (atomic, cell, and fragment data)
    !> @param[inout] md   --  MD trajectory data
    !> @param[inout] dips --  Dipole data (dipole type, dipole files, etc.)
    !> @param[inout] rams --  Raman data structure (provides polarizability and field information)
    !>
    SUBROUTINE spec_raman(gs, sys, md, dips, rams)
        TYPE(global_settings), INTENT(INOUT)        :: gs
        TYPE(systems), INTENT(INOUT)        :: sys
        TYPE(molecular_dynamics), INTENT(INOUT)        :: md
        TYPE(dipoles), INTENT(INOUT)        :: dips
        TYPE(raman), INTENT(INOUT)        :: rams

        CHARACTER(LEN=str_len)                               :: output_fname, f_name, idx, outfile, c_label
        CHARACTER(LEN=str_len)                               :: msg, fname
        INTEGER                                                  :: stat, i, j, xyz, runit, i_laser, i_group
        INTEGER(kind=dp)                                          :: plan
        INTEGER, DIMENSION(:), ALLOCATABLE                         :: natom_frag_free
        INTEGER, DIMENSION(:, :, :), ALLOCATABLE                     :: fragment_free
        REAL(kind=dp)                                             :: f, freq_res, freq_range
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE                    :: zhat_unpol_x, zhat_depol_x, zhat_para_all, zhat_depol
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE                    :: raman_ortho, raman_para, raman_depol, raman_unpol
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE                    :: raman_const, sinc_func, freq
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE                :: dip_free
        REAL(kind=dp), DIMENSION(:, :, :, :), ALLOCATABLE                :: v_frag, dipole_frag, dip_frag_free
        COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE                 :: zhat_iso, zhat_aniso

        !Allocate
        ALLOCATE (rams%e_field(1)%alpha_xyz(sys%framecount, sys%mol_num, 3))
        ALLOCATE (rams%e_field(2)%alpha_xyz(sys%framecount, sys%mol_num, 3))
        ALLOCATE (rams%e_field(3)%alpha_xyz(sys%framecount, sys%mol_num, 3))
        ALLOCATE (freq(0:md%t_cor*2), raman_const(0:md%t_cor*2))

!    IF (rams%averaging=='1') THEN

        !!If there are no fragments
        IF (.NOT. sys%fragments%frag) THEN
            sys%fragments%ngroup = 1
        END IF

        !Determine the maximum frequency range in cm^{-1} based on `md%dt`
        freq_range = REAL((1._dp/(md%dt*fs2s))/speed_light, dp)
        !Determine the frequency resolution
        freq_res = REAL(freq_range/(2.0_dp*md%t_cor), kind=dp)
        !Apply sinc function
        f = freq_res*md%dt*sinc_factor

        !!Start with looping over fragments if there are any
        DO i_group = 1, sys%fragments%ngroup
        !!First field-free (unperturbed dipole moments)
            IF (dips%type_dipole=='wannier') THEN
                CALL read_coord_frame(sys%natom, dips%dip_file, md%coord_v, sys)
                IF (.NOT. sys%fragments%frag) THEN
                    CALL compute_dipole(dip_free, sys, md)
                !!If there are fragments, assign their wannier centers first
                ELSE
                    CALL assign_wannier(sys, md)
                    CALL compute_dipole_frag(dip_frag_free, sys, md)
                END IF
            ELSEIF (dips%type_dipole=='berry') THEN
                CALL read_coord_frame(sys%natom, dips%dip_file, dip_free, sys)
                !!Correct discontinuities due to polarization quantum, if there are any
                CALL check_jumps(dip_free, sys%mol_num, sys)
            END IF

        !!Then electric field (perturbed dipole moments)
            DO xyz = 1, 3 !X-FIELD Y-FIELD Z-FIELD
                IF (xyz==1) THEN
                    f_name = dips%dip_x_file
                ELSE IF (xyz==2) THEN
                    f_name = dips%dip_y_file
                ELSE IF (xyz==3) THEN
                    f_name = dips%dip_z_file
                END IF
                CALL read_coord_frame(sys%natom, f_name, md%coord_v, sys)
                IF (dips%type_dipole=='wannier') THEN
                    IF (.NOT. sys%fragments%frag) THEN
                        CALL compute_dipole(dips%dipole, sys, md)
                        !!Calculate polarizabilities
                        CALL forward_diff(sys%mol_num, rams%e_field(xyz)%alpha_xyz, dip_free, dips%dipole, gs, sys, dips)
                        DEALLOCATE (dips%dipole)
                    !!If there are fragments, assign their wannier centers first
                    ELSE
                        CALL assign_wannier(sys, md)
                        CALL compute_dipole_frag(dipole_frag, sys, md)
                        sys%mol_num = sys%fragments%type_frag(i_group)%nfrag
                        !!Calculate polarizabilities
                        CALL forward_diff(sys%mol_num, rams%e_field(xyz)%alpha_xyz, &
                                          dip_frag_free(i_group, :, 1:sys%mol_num, :), dipole_frag(i_group, :, 1:sys%mol_num, :), gs, sys, dips)
                    END IF
                ELSEIF (dips%type_dipole=='berry') THEN
                    CALL check_jumps(md%coord_v, sys%mol_num, sys)
                    !!Calculate polarizabilities
                    CALL forward_diff(sys%mol_num, rams%e_field(xyz)%alpha_xyz, dip_free, md%coord_v, gs, sys, dips)
                ELSEIF (dips%type_dipole=='dfpt') THEN
                    rams%e_field(xyz)%alpha_xyz = REAL(md%coord_v*a3_to_debye_per_e, kind=dp)
                END IF
                !!Time derivatives of the polarizabilities
                CALL central_diff(sys%mol_num, rams%e_field(xyz)%alpha_xyz, rams%e_field(xyz)%alpha_diff_xyz, sys, md)
            END DO

            !!Allocate
            ALLOCATE (zhat_iso(0:md%t_cor*2), zhat_aniso(0:md%t_cor*2))
            ALLOCATE (raman_ortho(0:md%t_cor*2), raman_para(0:md%t_cor*2), raman_depol(0:md%t_cor*2), raman_unpol(0:md%t_cor*2))

            !!Initialize
            zhat_iso = COMPLEX(0._dp, 0.0_dp)
            zhat_aniso = COMPLEX(0._dp, 0.0_dp)
            raman_ortho = 0.0d0
            raman_para = 0.0d0
            raman_unpol = 0.0d0
            raman_depol = 0.0d0
            !!Isotropic polarizability autocorrelation function
            CALL cvv_iso(sys%mol_num, rams%z_iso, rams%e_field(1)%alpha_diff_xyz, rams%e_field(2)%alpha_diff_xyz, rams%e_field(3)%alpha_diff_xyz, sys, md)
            !!Call FFT
            CALL dfftw_plan_dft_r2c_1d(plan, 2*md%t_cor, rams%z_iso, zhat_iso, FFTW_ESTIMATE)
            CALL dfftw_execute_dft_r2c(plan, rams%z_iso, zhat_iso)
            CALL dfftw_destroy_plan(plan)

            !!Anisotropic polarizability autocorrelation function
            CALL cvv_aniso(sys%mol_num, rams%z_aniso, rams%e_field(1)%alpha_diff_xyz, rams%e_field(2)%alpha_diff_xyz, rams%e_field(3)%alpha_diff_xyz, sys, md)
            !!Call FFT
            CALL dfftw_plan_dft_r2c_1d(plan, 2*md%t_cor, rams%z_aniso, zhat_aniso, FFTW_ESTIMATE)
            CALL dfftw_execute_dft_r2c(plan, rams%z_aniso, zhat_aniso)
            CALL dfftw_destroy_plan(plan)

            !!Unit conversion of Debye^2/(E^2*s^2) into C^4*s^2/kg^2
            zhat_iso(:) = REAL(zhat_iso(:), KIND=dp)*debye2cm*debye2cm/(au2vm*au2vm)
            zhat_aniso(:) = REAL(zhat_aniso(:), KIND=dp)*debye2cm*debye2cm/(au2vm*au2vm)

            !!Calculate the Raman constants for a range of laser wavelengths requested by user
            DO i_laser = 1, SIZE(rams%laser_in)
                DO i = 0, 2*md%t_cor - 2
                    freq(i) = i*freq_res
            !!conversion of the Raman intensities into m^2*K*cm*10^-30!!
                    raman_const(i) = const_planck/(8.0_dp*const_boltz*const_permit*const_permit) &
                                     *1.e+30*md%dt*fs2s*((((rams%laser_in(i_laser)/reccm2ev - freq(i))/cm2m)**4)/freq(i))* &
                                     (1.0_dp/(1.0_dp - EXP(-1._dp*const_planck*speed_light*cm2m*freq(i)/ &
                                                           (const_boltz*gs%temp))))*2.0_dp
                END DO

                !!!Orthogonal Raman intensities, written out only if the user requested!!!

                DO i = 0, 2*md%t_cor - 2
                    IF (i_laser==1) THEN
                        zhat_aniso(i + 1) = REAL(zhat_aniso(i + 1), kind=dp)*(f*(i + 1)/SIN(f*(i + 1)))**2._dp
                    END IF
                    raman_ortho(i) = ((REAL(zhat_aniso(i), kind=dp))/15.0_dp)*raman_const(i)
                    raman_ortho(0) = 0.0_dp
                    IF (freq(i).GE.5000.0_dp) CYCLE
                END DO
                IF (gs%spectra_verbosity=='high') THEN
                    IF (.NOT. sys%fragments%frag) THEN
                        outfile = "raman_orthogonal.txt"
                    ELSE
                        WRITE (outfile, '(A,I0,A)') "raman_orthogonal_fragment_", i_group, ".txt"
                    END IF

                    WRITE (c_label, '("INT ",F10.6, " eV")') rams%laser_in(i_laser)
                    IF (i_laser==1) THEN
                        CALL write_spectra_data(outfile, c_label, freq, raman_ortho(:), 5000.0_dp)
                    ELSE
                        CALL append_column(outfile, c_label, raman_ortho(:), freq, 5000.0_dp)
                    END IF
                END IF

                !!!Parallel Raman intensities, written out only if the user requested!!!

                DO i = 0, 2*md%t_cor - 2
                    IF (i_laser==1) THEN
                        zhat_iso(i + 1) = REAL(zhat_iso(i + 1), kind=dp)*(f*(i + 1)/SIN(f*(i + 1)))**2._dp
                    END IF
                    raman_para(i) = (zhat_iso(i) + (zhat_aniso(i)*4.0_dp/45.0_dp))*raman_const(i)
                    raman_para(0) = 0.0_dp
                    IF (freq(i).GE.5000.0_dp) CYCLE
                    !WRITE (runit, *) freq(i), raman_para(i)
                END DO
                IF (gs%spectra_verbosity=='high') THEN
                    IF (.NOT. sys%fragments%frag) THEN
                        outfile = "raman_parallel.txt"
                    ELSE
                        WRITE (outfile, '(A,I0,A)') "raman_parallel_fragment_", i_group, ".txt"
                    END IF
                    WRITE (c_label, '("INT ",F10.6, " eV")') rams%laser_in(i_laser)
                    IF (i_laser==1) THEN
                        CALL write_spectra_data(outfile, c_label, freq, raman_para(:), 5000.0_dp)
                    ELSE
                        CALL append_column(outfile, c_label, raman_para(:), freq, 5000.0_dp)
                    END IF
                END IF

                !!Unpolarized Raman intensities, printed out by default
                DO i = 0, 2*md%t_cor - 2
                    raman_unpol(i) = raman_ortho(i) + raman_para(i)
                    raman_unpol(0) = 0.00_dp
                    IF (freq(i).GE.5000.0_dp) CYCLE
                    !WRITE (runit, *) freq(i), raman_unpol(i)
                END DO
                IF (gs%spectra_verbosity=='high' .OR. gs%spectra_verbosity=='normal') THEN
                    IF (.NOT. sys%fragments%frag) THEN
                        outfile = "raman_unpolarized.txt"
                    ELSE
                        WRITE (outfile, '(A,I0,A)') "raman_unpolarized_fragment_", i_group, ".txt"
                    END IF
                    WRITE (c_label, '("INT ",F10.6, " eV")') rams%laser_in(i_laser)
                    IF (i_laser==1) THEN
                        CALL write_spectra_data(outfile, c_label, freq, raman_unpol(:), 5000.0_dp)
                    ELSE
                        CALL append_column(outfile, c_label, raman_unpol(:), freq, 5000.0_dp)
                    END IF
                END IF

                !!!Depolarization ratios, written out only if the user requested!!!
                DO i = 0, 2*md%t_cor - 2
                    raman_depol(i) = REAL(raman_ortho(i), kind=dp)/REAL(raman_para(i), kind=dp)
                    raman_depol(0) = 0.00_dp
                    IF (freq(i).GE.5000.0_dp) CYCLE
                    !WRITE (runit, *) freq(i), raman_depol(i)
                END DO
                IF (gs%spectra_verbosity=='high') THEN
                    IF (.NOT. sys%fragments%frag) THEN
                        outfile = "raman_depolarization_ratio.txt"
                    ELSE
                        WRITE (outfile, '(A,I0,A)') "raman_depolarization_ratio_fragment_", i_group, ".txt"
                    END IF
                    WRITE (c_label, '("INT ",F10.6, " eV")') rams%laser_in(i_laser)
                    IF (i_laser==1) THEN
                        CALL write_spectra_data(outfile, c_label, freq, raman_depol(:), 5000.0_dp)
                    ELSE
                        CALL append_column(outfile, c_label, raman_depol(:), freq, 5000.0_dp)
                    END IF
                END IF

            END DO
            DEALLOCATE (rams%z_iso, rams%z_aniso)
            DEALLOCATE (raman_depol, raman_para, raman_unpol, raman_ortho, zhat_aniso, zhat_iso)
        END DO
        DEALLOCATE (freq, raman_const)
!    ELSEIF (rams%averaging=='2') THEN

!        IF (dips%type_dipole=='wannier') THEN
!            CALL read_coord_frame(sys%natom, rams%wannier_free, md%coord_v, sys)
        !sys%filename = rams%wannier_free
        !CALL read_coord_frame(sys, md)
!            CALL wannier(sys%filename, dip_free, sys, md)
!        ELSEIF (dips%type_dipole=='berry') THEN

        !           CALL read_coord_frame(sys%natom, rams%wannier_free, dip_free, sys)
        !       END IF

        !!!!ELECTRIC FIELD!!!
!        DO xyz = 1, 3 !X-FIELD Y-FIELD Z-FIELD
!            IF (rams%direction==CHAR(48 + xyz)) THEN !<-- maybe change direction to INTEGER input
!                CALL read_coord_frame(sys%natom, rams%e_field(xyz)%wannier_xyz, md%coord_v, sys)
!                IF (dips%type_dipole=='wannier') THEN
!                    CALL wannier(rams%e_field(xyz)%wannier_xyz,rams%e_field(xyz)%dip_xyz, sys, md)
!                    CALL forward_diff(sys%mol_num, rams%e_field(xyz)%alpha_xyz, dip_free,rams%e_field(xyz)%dip_xyz, gs, sys)
!                ELSEIF (dips%type_dipole=='berry') THEN
!                    CALL forward_diff(sys%mol_num, rams%e_field(xyz)%alpha_xyz, dip_free, md%coord_v, gs, sys)
!                END IF
!                CALL central_diff(sys%natom, rams%e_field(xyz)%alpha_xyz,rams%e_field(xyz)%alpha_diff_xyz, sys, md)
!            END IF
!        END DO

        !ALLOCATE (zhat_para(0:md%t_cor*2), zhat_unpol_x(0:md%t_cor*2), zhat_ortho(0:md%t_cor*2), zhat_depol_x(0:md%t_cor*2))

        ! zhat_para = COMPLEX(0._dp, 0.0_dp)
        ! zhat_ortho = COMPLEX(0._dp, 0.0_dp)
        ! zhat_unpol_x = COMPLEX(0._dp, 0.0_dp)
        ! zhat_depol_x = COMPLEX(0._dp, 0.0_dp)

!!IF ONLY ISOTROPIC AVERAGING IS CONSIDERED!!
!        CALL cvv_only_x(sys%mol_num, sys%framecount, rams%z_para, rams%z_ortho,rams%e_field(1)%alpha_diff_xyz, &
!                       rams%e_field(2)%alpha_diff_xyz,rams%e_field(3)%alpha_diff_xyz, rams%direction)

        !CALL dfftw_plan_dft_r2c_1d(plan, 2*md%t_cor, rams%z_para, zhat_para, FFTW_ESTIMATE)
        !CALL dfftw_execute_dft_r2c(plan, rams%z_para, zhat_para)

        ! CALL dfftw_plan_dft_r2c_1d(plan, 2*md%t_cor, rams%z_ortho, zhat_ortho, FFTW_ESTIMATE)
        ! CALL dfftw_execute_dft_r2c(plan, rams%z_ortho, zhat_ortho)

!!ORTHOGONAL!!
        ! OPEN (UNIT=68, FILE='result_fft_water_lib_ortho_iso.txt', STATUS='unknown', IOSTAT=stat)
        ! zhat_ortho = REAL(zhat_ortho, kind=dp)
        ! freq_res = REAL(md%freq_range/md%t_cor, kind=dp)
        ! f = freq_res*md%dt*1.883652d-4

        !DO i = 0, 2*md%t_cor - 2
        !  zhat_ortho(i + 1) = REAL(zhat_ortho(i + 1), kind=dp)*(f*(i + 1)/SIN(f*(i + 1)))**2._dp

        !   zhat_ortho(i) = REAL(zhat_ortho(i), kind=dp)*((const_planck)/(8.0_dp*const_boltz*const_permit*const_permit) &
        !                                                 *1d-29*0.421_dp*md%dt &
        !                                                 *(((rams%laser_in(1) - ((i)*freq_res))**4)/((i)*freq_res)) &
        !                                                 *(1.0_dp/(1.0_dp - EXP((-1.438777_dp*((i)*freq_res)) &
        !                                                                         /gs%temp))))*2.0_dp*pi*2.0_dp

        !     zhat_ortho(0) = 0.0_dp
        !     IF ((i*freq_res).GE.5000_dp) CYCLE
        !     WRITE (68, *) i*freq_res, (REAL(zhat_ortho(i), kind=dp))
        ! END DO
        ! CLOSE (68)

!!PARALLEL!!
        !OPEN (UNIT=67, FILE='result_fft_water_lib_para_iso.txt', STATUS='unknown', IOSTAT=stat)
        !zhat_para = REAL(zhat_para, kind=dp)
        !freq_res = REAL(md%freq_range/md%t_cor, kind=dp)
        !f = freq_res*md%dt*1.883652d-4

        !DO i = 0, 2*md%t_cor - 2
        !zhat_para(i + 1) = REAL(zhat_para(i + 1), kind=dp)*(f*(i + 1)/SIN(f*(i + 1)))**2._dp

        ! zhat_para(i) = REAL(zhat_para(i), kind=dp)*((const_planck)/(8.0_dp*const_boltz*const_permit*const_permit) &
        !                                            *1d-29*0.421_dp*md%dt &
        !                                            *(((rams%laser_in(1) - ((i)*freq_res))**4)/((i)*freq_res)) &
        !                                            *(1.0_dp/(1.0_dp - EXP((-1.438777_dp*((i)*freq_res)) &
        !                                                                    /gs%temp))))*2.0_dp*pi*2.0_dp
        !  zhat_para(0) = 0.0_dp
        !   IF ((i*freq_res).GE.5000_dp) CYCLE
        !    WRITE (67, *) (i)*freq_res, (REAL(zhat_para(i), kind=dp))!,REAL(integral(i),kind=dp)
        ! END DO
        ! CLOSE (67)

!!UNPOL!!
        ! OPEN (UNIT=69, FILE='result_fft_water_lib_unpol_iso.txt', STATUS='unknown', IOSTAT=stat)
        ! freq_res = REAL(md%freq_range/md%t_cor, kind=dp)

        !  DO i = 0, 2*md%t_cor - 2
        !    zhat_unpol_x(i) = zhat_para(i) + zhat_ortho(i)
        !     IF ((i*freq_res).GE.5000_dp) CYCLE
        !      WRITE (69, *) i*freq_res, REAL(zhat_unpol_x(i), kind=dp)
        !   END DO
        !   CLOSE (69)

!!DEPOL RATIO!!
        ! OPEN (UNIT=70, FILE='result_fft_water_lib_depol_iso.txt', STATUS='unknown', IOSTAT=stat)

        ! DO i = 0, 2*md%t_cor - 2
        !    zhat_depol_x(i) = REAL(zhat_ortho(i), kind=dp)/REAL(zhat_para(i), kind=dp)
        !    IF ((i*freq_res).GE.5000_dp) CYCLE
        !     WRITE (70, *) i*freq_res, REAL(zhat_depol_x(i), kind=dp)!REAL(zhat_ortho(i),kind=dp)/REAL(zhat_para(i),kind=dp)
        !  END DO
        !   CLOSE (70)

        !    DEALLOCATE (rams%z_para, rams%z_ortho, zhat_unpol_x, zhat_depol_x, zhat_para, zhat_ortho)
        ! END IF

        ! CALL dfftw_destroy_plan(plan)

        !DEALLOCATE(dip_free,dip_xyz(1),dip_xyz(2),dip_xyz(3))
!    IF (sys%system=='1') THEN
!        DEALLOCATE ( fragment_free)
!        !  DEALLOCATE(natom_frag_x,natom_frag_y,natom_frag_z,natom_frag_free)
!    END IF

        !DEALLOCATE (alpha_xyz)
        !DEALLOCATE (alpha_diff_xyz)

    END SUBROUTINE spec_raman

!***********************************************************************************************!
!***********************************************************************************************!
    !> @brief Performs normal mode analysis via finite-difference evaluation of the Hessian.
    !>
    !> This subroutine constructs the mass-weighted Hessian matrix using finite differences
    !> of atomic forces, diagonalizes it via LAPACK’s `DSYEV` routine, and extracts
    !> the vibrational frequencies and normal mode displacement vectors.
    !>
    !> @param[inout] sys   -- System structure (provides atomic and mass information).
    !> @param[inout] stats -- Static calculation structure (provides force data, step size,
    !>                         and arrays for frequencies and displacements).
    !>
    SUBROUTINE normal_mode_analysis(sys, stats, gs)
        TYPE(global_settings), INTENT(INOUT)        :: gs
        TYPE(systems), INTENT(INOUT)        :: sys
        TYPE(static), INTENT(INOUT)        :: stats

        CHARACTER(len=str_len)                                         :: msg, fname
        LOGICAL, DIMENSION(9)                                        :: mk = .TRUE.
        INTEGER                                                     :: stat, i, j, m, n, p, k, info, lwork, lwmax, lda, runit
        REAL(kind=dp)                                               :: fin_diff_factor
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE                       :: work
        REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE                     :: hessian, atomic_displacements

        ! Conservative guess; LAPACK recommends this for DSYEV
        lwmax = MAX(1, 3*sys%natom*64)
        lda = sys%natom*3
        stats%nmodes = lda

        !Allocate
        ALLOCATE (work(lwmax), hessian(0:sys%natom*3 - 1, 0:sys%natom*3 - 1))
        ALLOCATE (stats%freq(stats%nmodes), atomic_displacements(stats%nmodes, sys%natom*3), stats%disp(stats%nmodes, sys%natom, 3))

        !!Finite difference factor
        fin_diff_factor = 1.0_dp/(2.0_dp*stats%dx)

        !!Construct hessian
        p = 0
        DO i = 0, sys%natom - 1
            DO m = 0, 2
                k = 0
                DO j = 0, sys%natom - 1
                    DO n = 0, 2
                        hessian(i + m + p, j + n + k) = sys%mass_mat(i + 1, j + 1)*hartreebohr2evang*fin_diff_factor*hessian_factor* &
                                                        (stats%force(j + 1, n + 1)%atom(i + 1)%displacement(2)%XYZ(m + 1)%frame(1) - stats%force(j + 1, n + 1)%atom(i + 1)%displacement(1)%XYZ(m + 1)%frame(1))
                    END DO
                    k = k + 2
                END DO
            END DO
            p = p + 2
        END DO
        !!Symmetrize
        hessian(:, :) = REAL((hessian(:, :) + TRANSPOSE(hessian(:, :)))/2.0_dp, kind=dp)
        n = SIZE(hessian, 1)

        ! work size query
        lwork = -1
        CALL DSYEV('V', 'U', n, hessian, lda, stats%freq, work, lwork, info)
        lwork = MIN(lwmax, INT(work(1)))

        ! get eigenvalues and eigenvectors
        CALL dsyev('V', 'U', n, hessian, lda, stats%freq, work, lwork, info)

        hessian = TRANSPOSE(hessian)

        stats%freq = REAL(stats%freq*SQRT(ABS(stats%freq))/ABS(stats%freq), kind=dp)
        !!Unit conversion for the frequencies
        stats%freq = REAL(stats%freq/(2.0_dp*pi*speed_light), kind=dp)


        atomic_displacements(1:stats%nmodes, 1:sys%natom*3) = hessian(0:3*sys%natom - 1, :)

        m = 0
        DO j = 0, sys%natom - 1
            DO k = 0, 2 !dims
                stats%disp(1:stats%nmodes, j + 1, k + 1) = atomic_displacements(1:stats%nmodes, j + k + 1 + m)
            END DO
            m = m + 2
        END DO

        !!Write the normal mode frequencies to a file
        OPEN (FILE='normal_mode_freq.txt', STATUS='unknown', ACTION='write', IOSTAT=stat, IOMSG=msg, NEWUNIT=runit)
        !Check if file exists
        CALL check_file_open(stat, msg, 'normal_mode_freq.txt')
        DO i = 1, stats%nmodes
            WRITE (runit, *) stats%freq(i)
        END DO
        CLOSE (runit)

        !!Write the normal mode displacements to a file
        OPEN (FILE='normal_mode_displ.txt', STATUS='unknown', ACTION='write', IOSTAT=stat, IOMSG=msg, NEWUNIT=runit)
        !Check if file exists
        CALL check_file_open(stat, msg, 'normal_mode_displ.txt')
        DO i = 1, stats%nmodes
            DO j = 1, sys%natom 
                WRITE (runit, *) stats%disp(i, j, 1:3)
            END DO
        END DO
        !!Write .mol output
        IF (gs%spectral_type%read_function=='NMA') THEN
            IF (stats%write_mol) THEN
                fname = "vibrations.mol"
                CALL write_mol_file(sys, stats, gs, fname)
            ENDIF
        ENDIF
        CLOSE (runit)
        
        DEALLOCATE(work, hessian, atomic_displacements)
    END SUBROUTINE normal_mode_analysis

!***********************************************************************************************!
!***********************************************************************************************!

    !> @brief Computes the static infrared (IR) spectrum.
    !>
    !> This subroutine calculates the IR spectrum using precomputed normal mode
    !> frequencies and dipole derivatives (∂μ/∂Q) obtained from finite-difference
    !> static calculations. It converts dipole derivatives to intensities, applies
    !> Gaussian broadening, and writes the resulting continuous IR spectrum to file.
    !>
    !> @param[inout] gs    -- Global settings (contains FWHM and constants)
    !> @param[inout] sys   -- System data (atomic and molecular information)
    !> @param[inout] stats -- Static calculation results (frequencies, displacements)
    !> @param[inout] dips  -- Dipole derivatives along normal modes (∂μ/∂Q)
    !>
    SUBROUTINE spec_static_ir(gs, sys, stats, dips)

        TYPE(global_settings), INTENT(INOUT)        :: gs
        TYPE(systems), INTENT(INOUT)        :: sys
        TYPE(static), INTENT(INOUT)        :: stats
        TYPE(dipoles), INTENT(INOUT)        :: dips

        CHARACTER(len=str_len)                                         :: msg, fname, c_label, spectra_file_name
        INTEGER                                                  :: stat, i, j, k, x, freq_res, runit
        INTEGER                                                  :: start_freq, end_freq
        REAL(kind=dp)                                             :: broad, ir_factor
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE        :: ir_int
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE                    :: gamma_sq, data2, freq!,broad

        !!Allocate
        ALLOCATE (gamma_sq(stats%nmodes), ir_int(stats%nmodes))
        !!Estimate frequency range
        start_freq = 1
        end_freq = INT(MAXVAL(stats%freq) + 1000.0_dp)
        freq_res = INT(end_freq - start_freq)

        !!Allcate based on the freq range
        ALLOCATE (data2(freq_res + 1))
        ALLOCATE (freq(freq_res + 1))

        data2 = 0.0_dp

        !!Dipole derivatives
        DO k = 1, stats%nmodes
            gamma_sq(k) = SQRT(DOT_PRODUCT(dips%dip_dq(k, :), dips%dip_dq(k, :)))
        END DO

        ! first convert debye²angstrom⁻²amu⁻¹ to C^2/kg then to km/mol
        ir_factor = (debye2cm/ang)**2.0_dp*avo_num*1.0e-3_dp/(12.0_dp*const_permit*(speed_light*cm2m)**2.0_dp)/am_u
        ir_int(:) = (gamma_sq(:)**2.0_dp)*ir_factor

        !Broadening the spectrum
        freq(:) = 0.0_dp
        DO i = start_freq, end_freq
            broad = 0.0_dp
            DO x = 1, stats%nmodes
                broad = broad + (ir_int(x)*(1.0_dp/(gs%fwhm*SQRT(2.0_dp*pi)))*EXP(-0.50_dp*((i - stats%freq(x))/gs%fwhm)**2.0_dp))
            END DO
            data2(i) = data2(i) + broad
            freq(:) = i
        END DO

        !!Write the results to a file
        c_label = "INT"
        spectra_file_name = "result_static_ir.txt"
        CALL write_spectra_data(spectra_file_name, c_label, freq, data2(:))

        !!Write .mol output
        IF (stats%write_mol) THEN
            fname = "IR.mol"
            CALL write_mol_file(sys, stats, gs, fname, ir_int)
        END IF
        DEALLOCATE (gamma_sq, data2, ir_int, dips%dip_dq, stats%freq, stats%disp)

    END SUBROUTINE spec_static_ir

!***********************************************************************************************!
!***********************************************************************************************!

    !> @brief Computes the static Raman spectrum from polarizability derivatives.
    !>
    !> This subroutine calculates the Raman intensities for each normal mode using the
    !> isotropic and anisotropic components of the polarizability derivative tensor
    !> (∂α/∂Q), applies Gaussian broadening, and writes the resulting Raman spectra
    !> to file.
    !>
    !> @param[inout] gs    -- Global settings (constants, temperature, FWHM)
    !> @param[inout] sys   -- System structure (atomic data, coordinates)
    !> @param[inout] stats -- Static calculation results (normal mode frequencies and displacements)
    !> @param[inout] dips  -- Dipole structure (included for interface consistency)
    !> @param[inout] rams  -- Raman structure (polarizability derivatives and laser frequencies)
    !>
    SUBROUTINE spec_static_raman(gs, sys, stats, dips, rams)

        TYPE(global_settings), INTENT(INOUT)        :: gs
        TYPE(systems), INTENT(INOUT)        :: sys
        TYPE(static), INTENT(INOUT)        :: stats
        TYPE(dipoles), INTENT(INOUT)        :: dips
        TYPE(raman), INTENT(INOUT)        :: rams

        CHARACTER(len=str_len)                                      :: msg, fname, outfile, c_label
        LOGICAL                                                   :: first_column
        INTEGER                                                  :: stat, i, j, x, freq_res, runit, i_laser
        INTEGER                                                  :: start_freq, end_freq!, recl
        REAL(kind=dp)                                             :: broad
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE                    :: iso_sq, aniso_sq, ram_const, data2, freq!,broad

        !!Allocate
        ALLOCATE (iso_sq(stats%nmodes), aniso_sq(stats%nmodes))
        ALLOCATE (rams%raman_int(stats%nmodes), ram_const(stats%nmodes))

        !!Estimate frequency range
        start_freq = 1
        end_freq = INT(MAXVAL(stats%freq) + 1000.0_dp)
        freq_res = INT(end_freq - start_freq)
        !!Allcate based on the freq range
        ALLOCATE (data2(freq_res + 1))
        ALLOCATE (freq(freq_res + 1))
        data2 = 0.0_dp
        ! recl = 4096

        outfile = 'result_static_raman.txt'

        !!!Isotropic and anisotropic contributions!!
        iso_sq(:) = REAL((rams%pol_dq(:, 1, 1) + rams%pol_dq(:, 2, 2) + rams%pol_dq(:, 3, 3))/3.0_dp, kind=dp)**2.0_dp

        aniso_sq(:) = (0.50_dp*(((rams%pol_dq(:, 1, 1) - rams%pol_dq(:, 2, 2))**2.0_dp) + &
                                ((rams%pol_dq(:, 2, 2) - rams%pol_dq(:, 3, 3))**2.0_dp) &
                                + ((rams%pol_dq(:, 3, 3) - rams%pol_dq(:, 1, 1))**2.0_dp))) &
                      + (3.0_dp*((rams%pol_dq(:, 1, 2)**2.0_dp) + (rams%pol_dq(:, 2, 3)**2.0_dp) &
                                 + (rams%pol_dq(:, 3, 1)**2.0_dp)))

        !!!Conversion from angstrom^4 amu⁻¹ to m^4 kg -1
        iso_sq = iso_sq*(ang**4._dp)/am_u
        aniso_sq = aniso_sq*(ang**4._dp)/am_u
        !!Different laser wavelengths
        DO i_laser = 1, SIZE(rams%laser_in)

            !!! Conversion of static Raman units into 10^{-30}*cm^2/sr
            ram_const(:) = (const_planck/(8.0_dp*speed_light*cm2m*const_permit*const_permit)*1.e+30* &
                            REAL(((rams%laser_in(i_laser)/reccm2ev - stats%freq(:))**4.0_dp)/(stats%freq(:)*cm2m**3.0_dp), kind=dp)* &
                            (1.0_dp/(1.0_dp - EXP(-1._dp*const_planck*speed_light*cm2m*stats%freq(:)/ &
                                                  (const_boltz*gs%temp)))))/(cm2m**2._dp)

            !!! Unpolarized Raman intensities
            rams%raman_int(:) = REAL(((7.0_dp*aniso_sq(:)) + (45.0_dp*iso_sq(:)))/45.0_dp, kind=dp)*ram_const(:)

            data2(:) = 0.0_dp
            freq(:) = 0.0_dp
            !!! Apply Gaussian broadening
            DO i = start_freq, end_freq
                broad = 0.0_dp
                DO x = 1, stats%nmodes
                    broad = broad + (rams%raman_int(x)*(1.0_dp/(gs%fwhm*SQRT(2.0_dp*pi)))* &
                                     EXP(-0.50_dp*((i - stats%freq(x))/gs%fwhm)**2.0_dp))
                END DO
                freq(i) = i
                data2(i) = broad
            END DO

            !!Write the results to a file
            WRITE (c_label, '("INT ",F10.6, " eV")') rams%laser_in(i_laser)
            IF (i_laser==1) THEN
                !CALL write_spectra_data(outfile, c_label, start_freq, end_freq, data2(:))
                CALL write_spectra_data(outfile, c_label, freq, data2(:))
            ELSE
                !WRITE(c_label,'("INT ",F10.6, " eV")') rams%laser_in(i_laser)
                CALL append_column(outfile, c_label, data2(:), freq)
            END IF

        !!Write .mol output
        IF (stats%write_mol) THEN
            IF (SIZE(rams%laser_in)==1) THEN
                fname = "Raman.mol"
            ELSE
                WRITE(fname, '(A,G0.3,A)') 'Raman_', rams%laser_in(i_laser), '_eV.mol'
            ENDIF
            CALL write_mol_file(sys, stats, gs, fname, rams%raman_int)
        END IF

        END DO
        DEALLOCATE (iso_sq, aniso_sq, data2, ram_const)

    END SUBROUTINE spec_static_raman

!***********************************************************************************************!
!***********************************************************************************************!

    !> @brief Computes the absorption spectrum from RT-TDDFT polarizabilities.
    !>
    !> This subroutine processes time-dependent polarizability data obtained from
    !> real-time TDDFT (RT-TDDFT) simulations to produce absorption spectra.
    !> It performs FFTs on the polarizability tensor components, and optionally
    !> applies Padé interpolation to enhance the spectral resolution.
    !>
    !> @param[inout] gs    -- Global settings structure (defines simulation mode and constants)
    !> @param[inout] sys   -- System data structure (contains atomic information)
    !> @param[inout] dips  -- Dipole data (provides applied electric field)
    !> @param[inout] rams  -- Raman/RT-TDDFT data (contains polarizabilities and FFT settings)
    !>
    SUBROUTINE spec_abs(gs, sys, dips, rams)

        TYPE(global_settings), INTENT(INOUT)        :: gs
        TYPE(systems), INTENT(INOUT)        :: sys
        TYPE(dipoles), INTENT(INOUT)        :: dips
        TYPE(raman), INTENT(INOUT)        :: rams

        CHARACTER(LEN=256)                                     :: filename, msg, c_label, spectra_file_name
        INTEGER                                                       :: stat, i, j, k, m, x, o, dims, dir, runit
        INTEGER(KIND=dp)                                               :: plan
        REAL(KIND=dp)                                                  :: rtp_freq_res, freq_au
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE                            :: freq, abs_int
        REAL(KIND=dp), DIMENSION(:, :, :, :), ALLOCATABLE                   :: trace, abs_intens
        COMPLEX(KIND=dp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE            :: y_out

        !! Assign dimensions and displacement directions, if absorption spectrum is requested
        !! do not perform Pade or FFT for the all shifted structures.
        IF (gs%spectral_type%read_function=='RR') THEN
            dims = 3
            dir = 2
        ELSEIF (gs%spectral_type%read_function=='ABS') THEN
            dims = 1
            dir = 1
            sys%natom = 1
        END IF

        !Allocate
        ALLOCATE (rams%RR%zhat_pol_rtp(sys%natom, dims, dir, 3, 3, rams%RR%framecount_rtp))
        !Initialize
        rams%RR%zhat_pol_rtp = COMPLEX(0._dp, 0.0_dp)

        !!!FFT of the RTP polarizabilities
        DO j = 1, sys%natom
            DO i = 1, dims
                DO k = 1, dir
                    DO m = 1, 3
                        DO o = 1, 3
                            CALL dfftw_plan_dft_r2c_1d(plan, rams%RR%framecount_rtp, &
                                                       rams%RR%pol_rtp(m, o)%atom(j)%displacement(k)%XYZ(i)%frame(1:rams%RR%framecount_rtp), &
                                                       rams%RR%zhat_pol_rtp(j, i, k, m, o, 1:rams%RR%framecount_rtp), FFTW_ESTIMATE)
                            CALL dfftw_execute_dft_r2c(plan, &
                                                       rams%RR%pol_rtp(m, o)%atom(j)%displacement(k)%XYZ(i)%frame(1:rams%RR%framecount_rtp), &
                                                       rams%RR%zhat_pol_rtp(j, i, k, m, o, 1:rams%RR%framecount_rtp))
                            CALL dfftw_destroy_plan(plan)
                        END DO
                    END DO
                END DO
            END DO
        END DO

        !!If Pade interpolation is requested
        IF (rams%RR%check_pade=='y') THEN

            ALLOCATE (y_out(sys%natom, dims, dir, 3, 3, rams%RR%framecount_rtp_pade))
        !!Call Pade
!$OMP PARALLEL DO COLLAPSE(5)
            DO j = 1, sys%natom
                DO i = 1, dims
                    DO k = 1, dir
                        DO m = 1, 3
                            DO o = 1, 3
                                CALL interpolate(rams%RR%framecount_rtp, rams%RR%zhat_pol_rtp(j, i, k, m, o, 1:rams%RR%framecount_rtp), &
                                                 rams%RR%framecount_rtp_pade, y_out(j, i, k, m, o, :))
                            END DO
                        END DO
                    END DO
                END DO
            END DO
!$OMP END PARALLEL DO
            !!Reassign the polarizability arrays
            rams%RR%framecount_rtp = rams%RR%framecount_rtp_pade
            rams%RR%zhat_pol_rtp = y_out
            DEALLOCATE (y_out)
        END IF

!!!Dividing by electric field and multiplying by rams%RR%dt_rtp which is coming from FFT
        rams%RR%zhat_pol_rtp = rams%RR%zhat_pol_rtp*(rams%RR%dt_rtp*fs2s)/dips%e_field

        !!Find the maximum frequency range in cm^{-1} based on rams%RR%dt_rtp
        rams%RR%freq_range_rtp = REAL((1.0_dp/(rams%RR%dt_rtp*fs2s))/speed_light, kind=dp)
!!!Finding frequency range
        rtp_freq_res = REAL(rams%RR%freq_range_rtp/rams%RR%framecount_rtp, kind=dp)

!!!Calculate absorption spectra

        ALLOCATE (trace(sys%natom, dims, dir, rams%RR%framecount_rtp))
        ALLOCATE (abs_intens(sys%natom, dims, dir, rams%RR%framecount_rtp))
        trace = 0.0_dp
        trace(:, :, :, :) = DIMAG(rams%RR%zhat_pol_rtp(:, :, :, 1, 1, :)) + DIMAG(rams%RR%zhat_pol_rtp(:, :, :, 2, 2, :)) &
                            + DIMAG(rams%RR%zhat_pol_rtp(:, :, :, 3, 3, :))

      !!Conversion of absorption spectrum units into a.u.
        abs_intens(:, :, :, :) = (4.0_dp*pi*debye*trace(:, :, :, :))/(3.0_dp*speed_light_au*at_u)

      !! Conversion from cm-1 to a.u.
        freq_au = rtp_freq_res*(-1.0_dp)*reccm2au

        ALLOCATE (freq(rams%RR%framecount_rtp))
        ALLOCATE (abs_int(rams%RR%framecount_rtp))

        freq = 0.0_dp; abs_int = 0.0_dp
        !!Generate the absorption spectrum
        DO o = 1, rams%RR%framecount_rtp
            freq(o) = o*rtp_freq_res*reccm2ev
            abs_int(o) = abs_intens(1, 1, 1, o)*o*freq_au
        END DO
        !!Write the results to a file
        c_label = "INT"
        spectra_file_name = "absorption_spectrum.txt"
        CALL write_spectra_data(spectra_file_name, c_label, freq, abs_int, MAXVAL(freq) + 1)

        !CLOSE (runit)
        DEALLOCATE (trace, abs_intens)

    END SUBROUTINE spec_abs

!***********************************************************************************************!
!***********************************************************************************************!

    !> @brief Computes the static resonance Raman spectrum from frequency-dependent polarizabilities.
    !>
    !> This subroutine calculates the resonance Raman (RR) spectrum by evaluating
    !> frequency-dependent polarizability derivatives with respect to mass-weighted
    !> normal coordinates. The derivatives are obtained via finite differences of
    !> RT-TDDFT polarizabilities, projected along each vibrational normal mode, and
    !> used to compute isotropic and anisotropic contributions.
    !>
    !> @param[inout] gs    -- Global settings (temperature, broadening width, constants)
    !> @param[inout] sys   -- System structure (atoms, masses)
    !> @param[inout] stats -- Static results (normal mode frequencies, displacements)
    !> @param[inout] rams  -- Raman/RT-TDDFT data (polarizabilities, laser frequencies, time grids)
    !>
    SUBROUTINE spec_static_resraman(gs, sys, stats, rams)

        TYPE(global_settings), INTENT(INOUT)        :: gs
        TYPE(systems), INTENT(INOUT)        :: sys
        TYPE(static), INTENT(INOUT)        :: stats
        TYPE(raman), INTENT(INOUT)        :: rams

        INTEGER                                                       :: stat, i, j, k, m, x, freq_res, l, o, n, r, runit, i_laser
        INTEGER                                                       :: start_freq, end_freq, rtp_point
        INTEGER(kind=dp)                                               :: plan
        CHARACTER(len=str_len)                                         :: msg, fname, outfile, c_label
        REAL(kind=dp)                                                  :: rtp_freq_res, pade_freq_res, fin_diff_factor, broad
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE                         :: data2, ram_const, freq
        REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE                       :: iso_sq, aniso_sq, raman_int
        COMPLEX(kind=dp), DIMENSION(:, :, :, :), ALLOCATABLE                   :: zhat_pol_dq_rtp
        COMPLEX(kind=dp), DIMENSION(:, :, :, :, :), ALLOCATABLE                 :: zhat_pol_dxyz_rtp

        !!Estimate frequency range
        start_freq = 1.0_dp
        end_freq = INT(MAXVAL(stats%freq) + 1000.0_dp)
        freq_res = INT(end_freq - start_freq)

        !!Allocate
        ALLOCATE (data2(freq_res*stats%nmodes))
        ALLOCATE (freq(freq_res + 1))
        ALLOCATE (zhat_pol_dxyz_rtp(sys%natom, 3, 3, 3, rams%RR%framecount_rtp))
        ALLOCATE (zhat_pol_dq_rtp(stats%nmodes, 3, 3, rams%RR%framecount_rtp))
        ALLOCATE (iso_sq(stats%nmodes, rams%RR%framecount_rtp), aniso_sq(stats%nmodes, rams%RR%framecount_rtp))
        ALLOCATE (raman_int(stats%nmodes, rams%RR%framecount_rtp), ram_const(stats%nmodes))

        outfile = 'result_static_resraman.txt'

        !!Finite difference factor
        fin_diff_factor = 1.0_dp/(2.0_dp*stats%dx)
        !Initialize
        zhat_pol_dq_rtp = (0.0_dp, 0.0_dp)
        data2 = 0.0_dp; freq = 0.0_dp; broad = 0.0_dp

        !!Find the maximum frequency range in cm^{-1} based on rams%RR%dt_rtp
        rams%RR%freq_range_rtp = REAL((1.0_dp/(rams%RR%dt_rtp*fs2s))/speed_light, kind=dp)
!!!Finding laser frequency
        rtp_freq_res = REAL(rams%RR%freq_range_rtp/rams%RR%framecount_rtp, kind=dp)

!!!Finite differences
        zhat_pol_dxyz_rtp(:, :, :, :, :) = (rams%RR%zhat_pol_rtp(:, :, 2, :, :, :) &
                                            - rams%RR%zhat_pol_rtp(:, :, 1, :, :, :))*fin_diff_factor

!!!Derivatives w.r.t. mass weighted normal coordinates
        DO i = 1, stats%nmodes
            DO o = 1, rams%RR%framecount_rtp
                DO j = 1, sys%natom
                    zhat_pol_dq_rtp(i, :, :, o) = zhat_pol_dq_rtp(i, :, :, o) &
                                                  + (zhat_pol_dxyz_rtp(j, 1, :, :, o)*stats%disp(i, j, 1)*sys%atom_mass_inv_sqrt(j)) &
                                                  + (zhat_pol_dxyz_rtp(j, 2, :, :, o)*stats%disp(i, j, 2)*sys%atom_mass_inv_sqrt(j)) &
                                                  + (zhat_pol_dxyz_rtp(j, 3, :, :, o)*stats%disp(i, j, 3)*sys%atom_mass_inv_sqrt(j))
                END DO
            END DO
        END DO

!!!Isotropic and anisotropic contributions!!
        iso_sq(:, :) = ABS((zhat_pol_dq_rtp(:, 1, 1, :) + zhat_pol_dq_rtp(:, 2, 2, :) + zhat_pol_dq_rtp(:, 3, 3, :))/3.0_dp)**2

        aniso_sq(:, :) = 0.5_dp*(ABS(zhat_pol_dq_rtp(:, 1, 1, :) - zhat_pol_dq_rtp(:, 2, 2, :))**2 + &
                                 ABS(zhat_pol_dq_rtp(:, 2, 2, :) - zhat_pol_dq_rtp(:, 3, 3, :))**2 + &
                                 ABS(zhat_pol_dq_rtp(:, 3, 3, :) - zhat_pol_dq_rtp(:, 1, 1, :))**2) &
                         + 3.0_dp*(ABS(zhat_pol_dq_rtp(:, 1, 2, :))**2 + ABS(zhat_pol_dq_rtp(:, 2, 3, :))**2 + ABS(zhat_pol_dq_rtp(:, 3, 1, :))**2)

!!!Conversion from (debye/E)^2 angstrom^-2 amu⁻¹ to angstrom^6 angstrom^-2 amu⁻¹

        iso_sq = iso_sq/(a3_to_debye_per_e*a3_to_debye_per_e)
        aniso_sq = aniso_sq/(a3_to_debye_per_e*a3_to_debye_per_e)

        !!!Conversion from angstrom^4 amu⁻¹ to m^4 kg^-1
        iso_sq = iso_sq*(ang**4._dp)/am_u
        aniso_sq = aniso_sq*(ang**4._dp)/am_u

        !!!Finding laser frequency
        rtp_freq_res = REAL(rams%RR%freq_range_rtp/rams%RR%framecount_rtp, kind=dp)

        DO i_laser = 1, SIZE(rams%laser_in)

            rtp_point = ANINT(rams%laser_in(i_laser)/(rtp_freq_res*reccm2ev), kind=dp)

            WRITE (*, '(4X,"rams%laser_in", T60, G0)') rams%laser_in(i_laser)
            !!! Conversion of static resonance Raman units into 10^{-30}*cm^2/sr
            ram_const(:) = (const_planck/(8.0_dp*speed_light*cm2m*const_permit*const_permit)*1.e+30* &
                            REAL(((rams%laser_in(i_laser)/reccm2ev - stats%freq(:))**4.0_dp)/(stats%freq(:)*cm2m**3.0_dp), kind=dp)* &
                            (1.0_dp/(1.0_dp - EXP(-1._dp*const_planck*speed_light*cm2m*stats%freq(:)/ &
                                                  (const_boltz*gs%temp)))))/(cm2m**2._dp)

            !!!Calculation of the unpolarized resonance Raman intensities!!
            raman_int(:, rtp_point) = REAL(((7.0_dp*aniso_sq(:, rtp_point)) + (45.0_dp*iso_sq(:, rtp_point)))/45.0_dp, kind=dp)* &
                                      ram_const(:)

            data2(:) = 0.0_dp
            freq(:) = 0.0_dp
            !!!Broadening the spectrum!!
            DO x = start_freq, end_freq
                broad = 0.0_dp
                DO i = 1, stats%nmodes
                    broad = broad + (raman_int(i, rtp_point)*(1.0_dp/(gs%fwhm*SQRT(2.0_dp*pi))) &
                                     *EXP(-0.50_dp*((x - stats%freq(i))/gs%fwhm)**2.0_dp))
                END DO
                data2(x) = broad
                freq(x) = x
            END DO

            !!Writing out the results
            WRITE (c_label, '("INT ",F10.6, " eV")') rams%laser_in(i_laser)
            IF (i_laser==1) THEN
                !CALL write_spectra_data(outfile, c_label, start_freq, end_freq, data2(:))
                CALL write_spectra_data(outfile, c_label, freq, data2(:))
            ELSE
                !WRITE(c_label,'("INT ",F10.6, " eV")') rams%laser_in(i_laser)
                CALL append_column(outfile, c_label, data2(:), freq)
            END IF
        
            !!Write .mol output
            IF (stats%write_mol) THEN
                IF (SIZE(rams%laser_in)==1) THEN
                    fname = "Resonance_Raman.mol"
                ELSE
                    WRITE(fname, '(A,G0.3,A)') 'Resonance_Raman_', rams%laser_in(i_laser), '_eV.mol'
                ENDIF
                CALL write_mol_file(sys, stats, gs, fname, raman_int(:, rtp_point))
            ENDIF
        END DO
        DEALLOCATE (iso_sq, aniso_sq, data2, ram_const, raman_int, stats%disp, stats%freq)
        DEALLOCATE (rams%RR%zhat_pol_rtp, zhat_pol_dxyz_rtp, zhat_pol_dq_rtp)

    END SUBROUTINE spec_static_resraman

!***********************************************************************************************!
!***********************************************************************************************!

    !> @brief Computes the MD-averaged absorption spectrum from time-dependent dipole data.
    !>
    !> This subroutine evaluates the absorption spectrum by performing Fourier transforms
    !> of time-domain dipole responses along MD trajectories.
    !> The routine computes polarizabilities under X, Y, and Z electric fields, applies
    !> FFT (and optional Padé interpolation), averages over MD snapshots, and converts
    !> the resulting spectra into absorption cross sections.
    !>
    !> @param[inout] gs   --  Global settings (defines FFT and Padé parameters)
    !> @param[inout] sys  --  System structure (MD snapshot data, atom count)
    !> @param[inout] md   --  Molecular dynamics control structure (trajectory information)
    !> @param[inout] rams --  Raman/RT-TDDFT data (time-dependent dipoles and polarizabilities)
    !> @param[inout] dips --  Dipole data (electric field amplitudes and orientation)
    !>
    SUBROUTINE spec_abs_md(gs, sys, md, rams, dips)

        TYPE(global_settings) :: gs
        TYPE(systems)         :: sys
        TYPE(molecular_dynamics)   :: md
        TYPE(dipoles)   :: dips
        TYPE(raman)   :: rams

        CHARACTER(LEN=str_len)                                        :: chara, msg, c_label, spectra_file_name
        INTEGER                                                  :: stat, i, j, k, m, t0, t1, runit
        INTEGER                                                  :: Nw, Nmd
        INTEGER(kind=dp)                                          :: plan, rtp_point
        REAL(kind=dp)                                             :: f, freq_res, rtp_freq_res, pade_freq_res, laser_in, freq_au
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE                     :: raman_const, sinc_func, freq
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE                    :: abs_intens, trace_pade, abs_intens_pade
        REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE                  :: trace
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE                :: alpha_x, alpha_y, alpha_z
        COMPLEX(kind=dp), DIMENSION(:, :, :), ALLOCATABLE             :: zhat_resraman_x, zhat_resraman_y, zhat_resraman_z
        COMPLEX(kind=dp), DIMENSION(:, :, :), ALLOCATABLE             :: y_out

        !!Remove the unperturbed dipole moment from the array
        Nw = rams%RR%framecount_rtp - 1
        Nmd = sys%framecount

        !!Allocate
        ALLOCATE (zhat_resraman_x(sys%framecount, Nw, 3), zhat_resraman_y(sys%framecount, Nw, 3), zhat_resraman_z(sys%framecount, Nw, 3))
        ALLOCATE (rams%RR%alpha_resraman_x_re(sys%framecount, rams%RR%framecount_rtp, 3), rams%RR%alpha_resraman_y_re(sys%framecount, rams%RR%framecount_rtp, 3), rams%RR%alpha_resraman_z_re(sys%framecount, rams%RR%framecount_rtp, 3))
        ALLOCATE (rams%RR%alpha_resraman_x_im(sys%framecount, rams%RR%framecount_rtp, 3), rams%RR%alpha_resraman_y_im(sys%framecount, rams%RR%framecount_rtp, 3), rams%RR%alpha_resraman_z_im(sys%framecount, rams%RR%framecount_rtp, 3))

        ALLOCATE (alpha_x(Nmd, Nw, 3), alpha_y(Nmd, Nw, 3), alpha_z(Nmd, Nw, 3))

        !!Initialize
        zhat_resraman_x = COMPLEX(0._dp, 0.0_dp)
        zhat_resraman_y = COMPLEX(0._dp, 0.0_dp)
        zhat_resraman_z = COMPLEX(0._dp, 0.0_dp)

!!!X-Field!!
        !!Calculate polarizabilities
        CALL forward_diff(rams%RR%framecount_rtp, alpha_x, rams%RR%dip_x_rtp, rams%RR%dip_x_rtp, gs, sys, dips, rams)
        !!Call FFT
        DO i = 1, sys%framecount
            DO j = 1, 3
                CALL dfftw_plan_dft_r2c_1d(plan, Nw, alpha_x(i, 1:Nw, j), zhat_resraman_x(i, 1:Nw, j), FFTW_ESTIMATE)
                CALL dfftw_execute_dft_r2c(plan, alpha_x(i, 1:Nw, j), zhat_resraman_x(i, 1:Nw, j)) !!!important to specify arrays!!
                CALL dfftw_destroy_plan(plan)
            END DO
        END DO
        !!If Pade interpolation is requested
        IF (rams%RR%check_pade=='y') THEN
            ALLOCATE (y_out(Nmd, rams%RR%framecount_rtp_pade, 3))
        !!Call Pade
!$OMP PARALLEL DO COLLAPSE(2)
            DO i = 1, Nmd
                DO j = 1, 3
                    CALL interpolate(Nw, zhat_resraman_x(i, 1:Nw, j), &
                                     rams%RR%framecount_rtp_pade, y_out(i, :, j))
                END DO
            END DO
!$OMP END PARALLEL DO
            zhat_resraman_x = y_out
            DEALLOCATE (y_out)
        END IF

        !!Multiply by dt_rtp coming from the FFT
        zhat_resraman_x = zhat_resraman_x*rams%RR%dt_rtp*fs2s
        !!Define real and imaginary components of the freq-dependent polarizabilities
        rams%RR%alpha_resraman_x_re = REAL(zhat_resraman_x, kind=dp)
        rams%RR%alpha_resraman_x_im = DIMAG(zhat_resraman_x)

!!!Y-Field!!
        CALL forward_diff(rams%RR%framecount_rtp, alpha_y, rams%RR%dip_y_rtp, rams%RR%dip_y_rtp, gs, sys, dips, rams)

        DO i = 1, sys%framecount
            DO j = 1, 3
                CALL dfftw_plan_dft_r2c_1d(plan, Nw, alpha_y(i, 1:Nw, j), zhat_resraman_y(i, 1:Nw, j), FFTW_ESTIMATE)
                CALL dfftw_execute_dft_r2c(plan, alpha_y(i, 1:Nw, j), zhat_resraman_y(i, 1:Nw, j)) !!!important to specify arrays!!
                CALL dfftw_destroy_plan(plan)
            END DO
        END DO

        IF (rams%RR%check_pade=='y') THEN
            ALLOCATE (y_out(Nmd, rams%RR%framecount_rtp_pade, 3))
        !!Call Pade
!$OMP PARALLEL DO COLLAPSE(2)
            DO i = 1, Nmd
                DO j = 1, 3
                    CALL interpolate(Nw, zhat_resraman_y(i, 1:Nw, j), &
                                     rams%RR%framecount_rtp_pade, y_out(i, :, j))
                END DO
            END DO
!$OMP END PARALLEL DO
            zhat_resraman_y = y_out
            DEALLOCATE (y_out)
        END IF

        !!Multiply by dt_rtp coming from the FFT
        zhat_resraman_y = zhat_resraman_y*rams%RR%dt_rtp*fs2s

        rams%RR%alpha_resraman_y_re = REAL(zhat_resraman_y, kind=dp)
        rams%RR%alpha_resraman_y_im = DIMAG(zhat_resraman_y)
!
!!!Z-Field!!
        CALL forward_diff(rams%RR%framecount_rtp, alpha_z, rams%RR%dip_z_rtp, rams%RR%dip_z_rtp, gs, sys, dips, rams)

        DO i = 1, sys%framecount
            DO j = 1, 3
                CALL dfftw_plan_dft_r2c_1d(plan, Nw, alpha_z(i, 1:Nw, j), zhat_resraman_z(i, 1:Nw, j), FFTW_ESTIMATE)
                CALL dfftw_execute_dft_r2c(plan, alpha_z(i, 1:Nw, j), zhat_resraman_z(i, 1:Nw, j)) !!!important to specify arrays!!
                CALL dfftw_destroy_plan(plan)
            END DO
        END DO

        IF (rams%RR%check_pade=='y') THEN
            ALLOCATE (y_out(Nmd, rams%RR%framecount_rtp_pade, 3))
        !!Call Pade
!$OMP PARALLEL DO COLLAPSE(2)
            DO i = 1, Nmd
                DO j = 1, 3
                    CALL interpolate(Nw, zhat_resraman_z(i, 1:Nw, j), &
                                     rams%RR%framecount_rtp_pade, y_out(i, :, j))
                END DO
            END DO
!$OMP END PARALLEL DO
            Nw = rams%RR%framecount_rtp_pade
            zhat_resraman_z = y_out
            DEALLOCATE (y_out)
        END IF

        !!Multiply by dt_rtp coming from the FFT
        zhat_resraman_z = zhat_resraman_z*rams%RR%dt_rtp*fs2s

        rams%RR%alpha_resraman_z_re = REAL(zhat_resraman_z, kind=dp)
        rams%RR%alpha_resraman_z_im = DIMAG(zhat_resraman_z)

!!!!Calculate absorption spectra
        ALLOCATE (trace(Nmd, Nw), abs_intens(Nw))

        trace = 0.0_dp
        trace(:, :) = rams%RR%alpha_resraman_x_im(:, :, 1) + rams%RR%alpha_resraman_y_im(:, :, 2) + rams%RR%alpha_resraman_z_im(:, :, 3)

        abs_intens(:) = SUM(trace(1:Nmd, :), DIM=1)/REAL(Nmd, dp) !!average over MD snapshots

        !!Conversion of absorption spectrum units into a.u.
        abs_intens(:) = abs_intens(:)*(4.0_dp*pi*debye)/(3.0_dp*speed_light_au*at_u)

        !!Find the maximum frequency range in cm^{-1} based on rams%RR%dt_rtp
        rams%RR%freq_range_rtp = REAL((1.0_dp/(rams%RR%dt_rtp*fs2s))/speed_light, kind=dp)
        !! Find the frequency resolution
        rtp_freq_res = REAL(rams%RR%freq_range_rtp/Nw, kind=dp)
        !! Conversion from cm-1 to a.u.
        freq_au = rtp_freq_res*reccm2au

        ALLOCATE (freq(Nw))
        DO i = 1, Nw
            abs_intens(i) = abs_intens(i)*i*freq_au
            freq(i) = i*rtp_freq_res*reccm2ev
        END DO
        c_label = "INT"
        spectra_file_name = "absorption_spectrum_md.txt"
        CALL write_spectra_data(spectra_file_name, c_label, freq, abs_intens, MAXVAL(freq) + 1)

        DEALLOCATE (zhat_resraman_x, zhat_resraman_y, zhat_resraman_z)
        DEALLOCATE (alpha_x, alpha_y, alpha_z)
        DEALLOCATE (trace, abs_intens)

    END SUBROUTINE spec_abs_md

!***********************************************************************************************!
!***********************************************************************************************!

    !> @brief Computes the resonance Raman spectrum from time-dependent polarizabilities.
    !>
    !> This subroutine calculates resonance Raman (RR) spectra using time-dependent
    !> polarizability tensors obtained from real-time TDDFT (RT-TDDFT) simulations,
    !> averaged over molecular dynamics (MD) trajectories. It performs central
    !> differentiation of polarizabilities, constructs isotropic and anisotropic
    !> autocorrelation functions, and evaluates the resulting Raman intensities for
    !> multiple laser excitation energies.
    !>
    !> @param[inout] gs    -- Global settings (temperature, FFT/sinc parameters, constants)
    !> @param[inout] sys   -- System structure (atoms, MD frame count)
    !> @param[inout] md    -- Molecular dynamics data (time step, correlation length)
    !> @param[inout] rams  -- Raman/RT-TDDFT data (polarizabilities, correlation arrays, laser energies)
    !> @param[inout] dips  -- Dipole/electric field information
    !>
    SUBROUTINE spec_resraman(gs, sys, md, rams, dips)

        TYPE(global_settings) :: gs
        TYPE(systems)         :: sys
        TYPE(molecular_dynamics)   :: md
        TYPE(dipoles)   :: dips
        TYPE(raman)   :: rams

        CHARACTER(LEN=str_len)                                      :: chara, msg, c_label, outfile
        INTEGER                                                  :: stat, i, j, k, m, runit, i_laser
        INTEGER                                                  :: Nw, Nmd
        INTEGER(kind=dp)                                          :: plan, rtp_point
        REAL(kind=dp)                                             :: f, freq_res, freq_range, rtp_freq_res, pade_freq_res, laser_in
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE                    :: trace, abs_intens, trace_pade, abs_intens_pade
        REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE                  :: zhat_unpol_resraman
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE                     :: raman_const, sinc_func, freq
        COMPLEX(kind=dp), DIMENSION(:, :), ALLOCATABLE               :: zhat_iso_resraman, zhat_aniso_resraman

        IF (rams%RR%check_pade=='n') THEN
            Nw = rams%RR%framecount_rtp - 1
        ELSEIF (rams%RR%check_pade=='y') THEN
            Nw = rams%RR%framecount_rtp_pade
        END IF

        Nmd = sys%framecount - 2

        !!Allocate
        ALLOCATE (rams%RR%alpha_resraman_x_diff_re(Nmd, Nw, 3), rams%RR%alpha_resraman_y_diff_re(Nmd, Nw, 3), rams%RR%alpha_resraman_z_diff_re(Nmd, Nw, 3))
        ALLOCATE (rams%RR%alpha_resraman_x_diff_im(Nmd, Nw, 3), rams%RR%alpha_resraman_y_diff_im(Nmd, Nw, 3), rams%RR%alpha_resraman_z_diff_im(Nmd, Nw, 3))
        ALLOCATE (freq(0:2*md%t_cor), raman_const(0:2*md%t_cor))

!!Time derivatives of polarizabilities under electric field in the X-direction
        CALL central_diff(Nw, rams%RR%alpha_resraman_x_re, rams%RR%alpha_resraman_x_diff_re, sys, md)
        CALL central_diff(Nw, rams%RR%alpha_resraman_x_im, rams%RR%alpha_resraman_x_diff_im, sys, md)

!!Time derivatives of polarizabilities under electric field in the Y-direction
        CALL central_diff(Nw, rams%RR%alpha_resraman_y_re, rams%RR%alpha_resraman_y_diff_re, sys, md)
        CALL central_diff(Nw, rams%RR%alpha_resraman_y_im, rams%RR%alpha_resraman_y_diff_im, sys, md)

!!Time derivatives of polarizabilities under electric field in the Z-direction
        CALL central_diff(Nw, rams%RR%alpha_resraman_z_re, rams%RR%alpha_resraman_z_diff_re, sys, md)
        CALL central_diff(Nw, rams%RR%alpha_resraman_z_im, rams%RR%alpha_resraman_z_diff_im, sys, md)

        !!Complex auto-correlation functions
        CALL cvv_resraman(sys, md, rams)

        !! Allocate
        ALLOCATE (zhat_iso_resraman(0:md%t_cor*2, Nw), zhat_aniso_resraman(0:md%t_cor*2, Nw))
        ALLOCATE (zhat_unpol_resraman(0:md%t_cor*2, Nw))
        !!Initialize
        zhat_iso_resraman = COMPLEX(0._dp, 0.0_dp)
        zhat_aniso_resraman = COMPLEX(0._dp, 0.0_dp)

        !!Find the maximum frequency range in cm^{-1} based on rams%RR%dt_rtp
        rams%RR%freq_range_rtp = REAL((1.0_dp/(rams%RR%dt_rtp*fs2s))/speed_light, kind=dp)
        !! Find the frequency resolution
        rtp_freq_res = REAL(rams%RR%freq_range_rtp/Nw, kind=dp)

        !!Loop over laser frequencies
        DO i_laser = 1, SIZE(rams%laser_in)
            rtp_point = ANINT(rams%laser_in(i_laser)/(rtp_freq_res*reccm2ev), kind=dp)
            CALL dfftw_plan_dft_1d(plan, 2*md%t_cor, rams%RR%z_iso_resraman(0:2*md%t_cor, rtp_point), zhat_iso_resraman(0:2*md%t_cor, rtp_point), &
                                   FFTW_FORWARD, FFTW_ESTIMATE)
            CALL dfftw_execute_dft(plan, rams%RR%z_iso_resraman(0:2*md%t_cor, rtp_point), zhat_iso_resraman(0:2*md%t_cor, rtp_point)) !!!important to specify arrays!!
            CALL dfftw_destroy_plan(plan)

            CALL dfftw_plan_dft_1d(plan, 2*md%t_cor, rams%RR%z_aniso_resraman(0:2*md%t_cor, rtp_point), zhat_aniso_resraman(0:2*md%t_cor, rtp_point), &
                                   FFTW_FORWARD, FFTW_ESTIMATE)
            CALL dfftw_execute_dft(plan, rams%RR%z_aniso_resraman(0:2*md%t_cor, rtp_point), zhat_aniso_resraman(0:2*md%t_cor, rtp_point)) !!!important to specify arrays!!
            CALL dfftw_destroy_plan(plan)

        END DO

        !Determine the maximum frequency range in cm^{-1} based on `md%dt`
        freq_range = REAL((1.0_dp/(md%dt*fs2s))/speed_light, kind=dp)
        !Determine the frequency resolution
        freq_res = REAL(freq_range/(2.0_dp*md%t_cor), kind=dp)
        !!Sinc factor
        f = freq_res*md%dt*sinc_factor

!!!!!Unpolarized Raman intensities
!Unit conversion of Debye^2/(E^2*s^2) into C^4*s^2/kg^2
        zhat_iso_resraman(:, :) = zhat_iso_resraman(:, :)*debye2cm*debye2cm/(au2vm*au2vm)
        zhat_aniso_resraman(:, :) = zhat_aniso_resraman(:, :)*debye2cm*debye2cm/(au2vm*au2vm)
        !!Loop over different laser frequencies
        DO i_laser = 1, SIZE(rams%laser_in)
            !!calculate RTP point for each laser wavelength
            rtp_point = ANINT(rams%laser_in(i_laser)/(rtp_freq_res*reccm2ev), kind=dp)
            DO i = 0, 2*md%t_cor - 1
                freq(i) = i*freq_res
                IF (i==0) THEN
                    raman_const(i) = 0.0_dp
                ELSE
         !!conversion of the Raman intensities into m^2*K*cm*10^-30!!
                    raman_const(i) = const_planck/(8.0_dp*const_boltz*const_permit*const_permit) &
                                     *1.e+30*md%dt*fs2s*((((rams%laser_in(i_laser)/reccm2ev - freq(i))/cm2m)**4)/freq(i))* &
                                     (1.0_dp/(1.0_dp - EXP(-1._dp*const_planck*speed_light*cm2m*freq(i)/ &
                                                           (const_boltz*gs%temp))))*2.0_dp
                END IF
                !!Apply sinc functions
                zhat_iso_resraman(i + 1, rtp_point) = (zhat_iso_resraman(i + 1, rtp_point))*(f*(i + 1)/SIN(f*(i + 1)))**2._dp
                zhat_aniso_resraman(i + 1, rtp_point) = (zhat_aniso_resraman(i + 1, rtp_point))*(f*(i + 1)/SIN(f*(i + 1)))**2._dp
                !!Generate the spectrum
                zhat_unpol_resraman(i, rtp_point) = REAL(zhat_iso_resraman(i, rtp_point) + (zhat_aniso_resraman(i, rtp_point)*7.0_dp/45.0_dp), KIND=dp)*raman_const(i)

                zhat_unpol_resraman(0, rtp_point) = 0.0_dp

            END DO

            !!Writing out the files
            outfile = "resonance_raman_unpol.txt"
            WRITE (c_label, '("INT ",F10.6, " eV")') rams%laser_in(i_laser)
            IF (i_laser==1) THEN
                CALL write_spectra_data(outfile, c_label, freq, zhat_unpol_resraman(:, rtp_point), 5000.0_dp)
            ELSE
                CALL append_column(outfile, c_label, zhat_unpol_resraman(:, rtp_point), freq, 5000.0_dp)
            END IF
        END DO

        DEALLOCATE (zhat_iso_resraman, zhat_aniso_resraman, zhat_unpol_resraman, freq, raman_const)

    END SUBROUTINE spec_resraman

END MODULE calc_spectra
