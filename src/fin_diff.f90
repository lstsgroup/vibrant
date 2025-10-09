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

!> @brief Module containing all procedures involving finite differences.
MODULE fin_diff

    USE kinds, ONLY: dp
    USE constants, ONLY: pi, debye, bohr2ang, speed_light, fs2s, joule_unit, ev_unit, action_unit
    USE vib_types, ONLY: global_settings, systems, molecular_dynamics, static, dipoles, raman
    USE cell_types, ONLY: build_hmat, pbc, invert3x3, determinant3x3
    IMPLICIT NONE
    PUBLIC :: central_diff, forward_diff, finite_diff_static, finite_diff_static_resraman

CONTAINS

!****************************************************************************!
!****************************************************************************!

    !> @brief Calculates time derivatives of a given parameter using a
    !>        central-difference scheme
    !>
    !> @param[in,out] sys    -- System information (provides framecount).
    !> @param[in,out] md     -- Molecular dynamics information (provides time step `dt`).
    !> @param[in,out] natom  -- Number of atoms for cartesian coordinates, or just `1`
    !>                          if the quantity (e.g. dipole or polarizability) refers
    !>                          to the entire system.
    !> @param[in,out] shifts -- 3D array of atomic displacements (or dipole moments or
    !>                          polarizability tensors), dimensions (nframes, natom, 3).
    !> @param[out]    diff   -- 3D array of computed central differences of the given
    !>                          parameter `shifts`, dimensions (nframes, natom, 3).
    !>                          Allocated if not present.
    !>
    SUBROUTINE central_diff(natom, shifts, diff, sys, md)

        TYPE(systems), INTENT(INOUT)        :: sys
        TYPE(molecular_dynamics), INTENT(INOUT)        :: md
        INTEGER, INTENT(INOUT)                                    :: natom
        REAL(kind=dp), DIMENSION(:, :, :), INTENT(INOUT)  ::  shifts
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE, INTENT(OUT)    ::  diff

        INTEGER                                                  :: stat, i, j, k, m
       
        !! Allocate 
        IF (.NOT. ALLOCATED(diff)) THEN
            ALLOCATE (diff(sys%framecount, natom, 3))
        END IF
        diff = 0.0_dp

        !! Central finite differences over MD frames
        DO j = 1, sys%framecount - 2
            DO i = 1, natom
                DO k = 1, 3
                    diff(j, i, k) = (shifts(j + 2, i, k) - shifts(j, i, k))/(2.0_dp*md%dt)
                END DO
            END DO
        END DO

    END SUBROUTINE central_diff
!**************************************************************************************************************!
!**************************************************************************************************************!

    !> @brief Calculates derivatives between two parameters using a forward finite-difference scheme.
    !>
    !> This routine computes finite differences between two sets of quantities to obtain their
    !> derivatives, typically used to calculate polarizabilities.
    !>
    !> Depending on `gs%spectral_type%read_function`, the routine applies either:
    !> - a simple forward-difference formula, or
    !> - a damped and Hann-windowed version (for time-dependent dipole moments)
    !>
    !> @param[in,out] sys      -- System information (provides `framecount`).
    !> @param[in,out] md       -- Molecular dynamics information (provides time step `dt`).
    !> @param[in,out] gs       -- Global setting (provides `spectral_type%read_function`).
    !> @param[in,out] dips     -- Dipole information (contains applied field strength `e_field`).
    !> @param[in,out] mol_num  -- Number of molecules, for cartesian coordinates it is equal to natom,
    !>                            if the dipole moment of the whole system is considered, it is `1`.
    !> @param[in]     dip_free -- 3D array of unperturbed dipole moments, dimensions (nframes, mol_num, 3).
    !> @param[in]     dip_x    -- 3D array of perturbed dipole moments, dimensions (nframes, mol_num, 3).
    !> @param[out]    alpha    -- 3D array of computed finite-difference derivatives (e.g., polarizabilities),
    !>                            dimensions (nframes, mol_num, 3). Allocated within the routine.
    !> @param[in,out] rams     -- Raman-related parameters (optional; required for spectral damping).
    !>
    SUBROUTINE forward_diff(mol_num, alpha, dip_free, dip_x, gs, sys, dips, rams)
        TYPE(global_settings), INTENT(INOUT) :: gs
        TYPE(systems), INTENT(INOUT)        :: sys
        TYPE(dipoles), INTENT(INOUT)        :: dips
        TYPE(raman), OPTIONAL         :: rams
        INTEGER, INTENT(INOUT)                                    :: mol_num
        REAL(kind=dp), DIMENSION(:, :, :), INTENT(IN)  :: dip_free, dip_x
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE, INTENT(OUT)    :: alpha

        INTEGER                                                  :: stat, i, j, k, m
        REAL(kind=dp)    :: conv_unit, damping_factor

        ALLOCATE (alpha(sys%framecount, mol_num, 3))

        !! IF the read_function is MD-RR or MD-ABS, do the unit conversion for the damping factor
        IF (gs%spectral_type%read_function.EQ.'MD-RR' .OR. gs%spectral_type%read_function.EQ.'MD-ABS') THEN
            conv_unit = rams%RR%damping_constant*joule_unit/ev_unit !! J
            damping_factor = conv_unit/action_unit*rams%RR%dt_rtp*fs2s !! s-1
        END IF

        DO j = 1, sys%framecount
            DO k = 1, 3
                !! IF the read_function is not MD-RR or MD-ABS, compute polarizabilities
                IF (gs%spectral_type%read_function.NE.'MD-RR' .AND. gs%spectral_type%read_function.NE.'MD-ABS') THEN
                    DO i = 1, mol_num  !!! change to mol_num later
                        alpha(j, i, k) = REAL((dip_x(j, i, k) - dip_free(j, i, k))/dips%e_field, kind=dp)
                    END DO
                !! IF the read_function is MD-RR or MD-ABS, compute time-dependent polarizabilities and apply damping
                ELSEIF (gs%spectral_type%read_function=='MD-RR' .OR. gs%spectral_type%read_function=='MD-ABS') THEN
                    DO i = 2, mol_num
                        alpha(j, i - 1, k) = REAL((dip_x(j, i, k) - dip_x(j, 1, k))*EXP(-1.0_dp*damping_factor*(i - 1)), KIND=dp)
                        alpha(j, i - 1, k) = alpha(j, i - 1, k)*((COS((i - 1)/((mol_num - 1) - 1.0_dp)/2.0_dp*pi))**2) !!Hann window function
                    END DO
                END IF
            END DO
        END DO

    END SUBROUTINE forward_diff

!**************************************************************************************************************!
!**************************************************************************************************************!
    !> @brief Calculates static finite-difference derivatives of dipole moments or polarizabilities.
    !>
    !> This routine computes the derivatives of dipole moments (for IR) or polarizabilities (for Raman) using
    !> a central finite-difference scheme along the normal-mode displacements. The resulting derivatives
    !> are used for obtaining Raman and IR intensities within the harmonic approximation.
    !>
    !> @param[in,out] gs     -- Global setting (provides `spectral_type%read_function`).
    !> @param[in,out] sys    -- System information (provides `natom`, inverse mass factors, etc.).
    !> @param[in,out] stats  -- Static spectral data (provides `nmodes`, normal mode displacements`disp`, step size `dx`).
    !> @param[in,out] dips   -- Dipole information (provides displaced dipoles and stores IR derivatives)
    !> @param[in,out] rams   -- Raman information (provides displaced polarizabilities and stores Raman derivatives)
    !>
    SUBROUTINE finite_diff_static(gs, sys, stats, dips, rams)

        TYPE(global_settings), INTENT(INOUT) :: gs
        TYPE(systems), INTENT(INOUT)   :: sys
        TYPE(static), INTENT(INOUT):: stats
        TYPE(dipoles), INTENT(INOUT)    ::  dips
        TYPE(raman), INTENT(INOUT)   :: rams

        INTEGER :: i_pol, j_pol, k, l, j, xyz
        REAL(kind=dp) :: delta, coeff, fin_diff_factor
        
        fin_diff_factor = 1.0_dp/(2.0_dp*stats%dx)
        !! Calculation of static Raman derivatives based on derivatives of polarizabilites wrt normal mode coordinates

        IF (gs%spectral_type%read_function=='R') THEN
            ALLOCATE (rams%pol_dq(stats%nmodes, 3, 3))

            rams%pol_dq = 0.0_dp

            DO i_pol = 1, 3
                DO j_pol = 1, 3
                    DO k = 1, sys%natom
                        DO l = 1, 3
                            coeff = fin_diff_factor*sys%atom_mass_inv_sqrt(k)
                            delta = rams%pol(i_pol, j_pol)%atom(k)%displacement(1)%XYZ(l)%frame(1) &
                                    - rams%pol(i_pol, j_pol)%atom(k)%displacement(2)%XYZ(l)%frame(1)

                            DO j = 1, stats%nmodes
                                rams%pol_dq(j, i_pol, j_pol) = rams%pol_dq(j, i_pol, j_pol) &
                                                               + coeff*stats%disp(j, k, l)*delta
                            END DO
                        END DO
                    END DO
                END DO
            END DO

        !! Calculation of static IR derivatives based on derivatives of dipole moments wrt normal mode coordinates
        ELSEIF (gs%spectral_type%read_function=='IR') THEN

            ALLOCATE (dips%dip_dq(stats%nmodes, 3))

            dips%dip_dq = 0.0_dp

            DO xyz = 1, 3
                DO k = 1, sys%natom
                    DO l = 1, 3
                        coeff = fin_diff_factor*sys%atom_mass_inv_sqrt(k)
                        delta = dips%static_dip(xyz)%atom(k)%displacement(2)%XYZ(l)%frame(1) &
                                - dips%static_dip(xyz)%atom(k)%displacement(1)%XYZ(l)%frame(1)

                        DO j = 1, stats%nmodes
                            dips%dip_dq(j, xyz) = dips%dip_dq(j, xyz) + coeff*stats%disp(j, k, l)*delta

                        END DO
                    END DO
                END DO
            END DO
        END IF

    END SUBROUTINE finite_diff_static

!**************************************************************************************************************!
!**************************************************************************************************************!
    !> @brief Constructs time-dependent polarizability responses for RR spectra using 
    !>        finite differences of static dipole moments.
    !>
    !> This routine computes the time-dependent polarizability tensor `pol_rtp` by taking
    !> the difference between the time-evolved and initial static dipole moments for
    !> each polarization direction. The resulting response is exponentially
    !> damped. The final polarizabilities are then used to compute RR intensities.
    !>
    !> @param[in,out] sys  -- System information (provides `natom`, inverse mass factors, etc.).
    !> @param[in,out] rams -- Raman data structure (provides static time-dependent dipole moments,
    !>                        resulting polarizabilities, RTP time step `rams%RR%dt_rtp`, and 
    !>                          RTP framecount `rams%RR%framecount_rtp`)
    !>
    SUBROUTINE finite_diff_static_resraman(sys, rams)

        TYPE(raman), INTENT(INOUT)        :: rams
        TYPE(systems), INTENT(INOUT)        :: sys

        INTEGER                              ::  x, y, i, j, k, m, l, i_pol, j_pol, stat
        INTEGER                              ::  dir = 2
        REAL(kind=dp)                        :: damping_factor, conv_unit

        !! Allocate
        CALL rams%RR%init_rr_pol(sys%natom, rams%RR%framecount_rtp)

        !! Calculation of the damping constant
        conv_unit = rams%RR%damping_constant*joule_unit/ev_unit !! J
        damping_factor = conv_unit/action_unit*rams%RR%dt_rtp*fs2s !! s-1

        !! Calculate the polarizabilities from time-dependent dipole moments
        DO k = 1, dir
            DO i = 1, sys%natom
                DO j = 1, 3
                    DO j_pol = 1, 3
                        DO l = 2, rams%RR%framecount_rtp + 1
                            rams%RR%pol_rtp(1, j_pol)%atom(i)%displacement(k)%XYZ(j)%frame(l - 1) = &
                            (rams%RR%static_dip_x_rtp(j_pol)%atom(i)%displacement(k)%XYZ(j)%frame(l) - &
                            rams%RR%static_dip_x_rtp(j_pol)%atom(i)%displacement(k)%XYZ(j)%frame(1))
                            rams%RR%pol_rtp(2, j_pol)%atom(i)%displacement(k)%XYZ(j)%frame(l - 1) = &
                            (rams%RR%static_dip_y_rtp(j_pol)%atom(i)%displacement(k)%XYZ(j)%frame(l) - &
                            rams%RR%static_dip_y_rtp(j_pol)%atom(i)%displacement(k)%XYZ(j)%frame(1))
                            rams%RR%pol_rtp(3, j_pol)%atom(i)%displacement(k)%XYZ(j)%frame(l - 1) = &
                            (rams%RR%static_dip_z_rtp(j_pol)%atom(i)%displacement(k)%XYZ(j)%frame(l) - &
                            rams%RR%static_dip_z_rtp(j_pol)%atom(i)%displacement(k)%XYZ(j)%frame(1))
                        END DO
                    END DO
                END DO
            END DO
        END DO

        !! Apply damping to polarizabilities
        DO k = 1, dir
            DO i = 1, sys%natom
                DO j = 1, 3
                    DO j_pol = 1, 3
                        DO l = 1, rams%RR%framecount_rtp
                            rams%RR%pol_rtp(1, j_pol)%atom(i)%displacement(k)%XYZ(j)%frame(l) = &
                            rams%RR%pol_rtp(1, j_pol)%atom(i)%displacement(k)%XYZ(j)%frame(l)*EXP(-1.0_dp*damping_factor*l)
                            rams%RR%pol_rtp(2, j_pol)%atom(i)%displacement(k)%XYZ(j)%frame(l) = &
                            rams%RR%pol_rtp(2, j_pol)%atom(i)%displacement(k)%XYZ(j)%frame(l)*EXP(-1.0_dp*damping_factor*l)
                            rams%RR%pol_rtp(3, j_pol)%atom(i)%displacement(k)%XYZ(j)%frame(l) = &
                            rams%RR%pol_rtp(3, j_pol)%atom(i)%displacement(k)%XYZ(j)%frame(l)*EXP(-1.0_dp*damping_factor*l)
                        END DO
                    END DO
                END DO
            END DO
        END DO

        !! Hann window function, commented-out for now
        ! DO i = 1, framecount_rtp
        !    pol_rtp(:,:,:,:,:,i) = pol_rtp(:,:,:,:,:,i)*((COS(i/(framecount_rtp - 1.0_dp)/2.0_dp*3.14_dp))**2) !!Hann Window function
        !END DO

    END SUBROUTINE finite_diff_static_resraman
END MODULE fin_diff
