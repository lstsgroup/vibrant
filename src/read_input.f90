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

!> @brief Module containing all procedures involving reading of the input file
MODULE read_input

    USE kinds, ONLY: dp, str_len
    USE vib_types!, ONLY: global_settings, systems, static, dipoles, raman, molecular_dynamics
    USE output_io, ONLY: check_file_open
    USE ISO_FORTRAN_ENV, ONLY: output_unit, error_unit
    USE dipole_calc, ONLY: check_lattice_parameters
    IMPLICIT NONE

    PRIVATE

    PUBLIC :: parse_command_line, parse_input, check_input

CONTAINS
!************************************************************************************!
!************************************************************************************!

    !> @brief Parses the command-line arguments to obtain the input file name.
    !>
    !> @param[out] input_file_name -- Name of the input file specified on the command line
    !>
    SUBROUTINE parse_command_line(input_file_name)

        CHARACTER(LEN=str_len), INTENT(OUT) :: input_file_name

        CHARACTER(LEN=str_len)                :: arg
        INTEGER(KIND=4)                                     :: narg
        INTEGER                                             :: i, stat

        narg = iargc()

        IF (narg/=1) THEN
            WRITE (*, '(T3, A37)') "Usage: vibrant_input.x your_input.inp"
            STOP
        END IF

        CALL getarg(1, input_file_name)

    END SUBROUTINE parse_command_line

!************************************************************************************!
!************************************************************************************!

    !> @brief Converts a given string to lowercase.
    !>
    !> @param[in] -- str         Input string to be converted.
    !> @return    -- lower_str   Lowercase version of the input string.
    !>
    FUNCTION to_lower(str) RESULT(lower_str)
        CHARACTER(len=*), INTENT(IN) :: str
        CHARACTER(len=LEN(str))      :: lower_str
        INTEGER :: i

        DO i = 1, LEN(str)
            SELECT CASE (str(i:i))
            CASE ('A':'Z')
                lower_str(i:i) = CHAR(IACHAR(str(i:i)) + 32)
            CASE DEFAULT
                lower_str(i:i) = str(i:i)
            END SELECT
        END DO
    END FUNCTION to_lower

    !> @brief Parses the Vibrant input file and initializes all simulation modules.
    !>
    !> This subroutine reads a user-specified input file and extracts settings for
    !> all simulation components. Furthermore, it identifies section headers in the input file such as
    !> `&global`, `&system`, `&cell`, `&static`, `&raman`, and so on, and
    !> processes all key–value pairs within each block until an `&end <section>`
    !> line is encountered.
    !>
    !> @param[in,out] gs    -- Global settings (provides temperature, FWHM, verbosity, spectra type, etc.).
    !> @param[in,out] sys   -- System information (provides coordinates, cell parameters, fragment groups, etc.).
    !> @param[in,out] md    --  Molecular dynamics information (provides time step `dt`, correlation depth `t_cor`).
    !> @param[in,out] stats -- Static spectral data (provides Hessian info, displacement step size `dx`, etc.).
    !> @param[in,out] dips  -- Dipole-related data (dipole file names, field strength, etc.).
    !> @param[in,out] rams  -- Raman-related parameters (laser frequencies, damping, Resonance Raman setup).
    !> @param[in]           -- input_file_name  Path to the input file provided by the user.
    !>
    SUBROUTINE parse_input(gs, sys, md, stats, dips, rams, input_file_name)

        TYPE(global_settings), INTENT(INOUT) :: gs
        TYPE(systems), INTENT(INOUT)        :: sys
        TYPE(molecular_dynamics), INTENT(INOUT)     :: md
        TYPE(static), INTENT(INOUT)           :: stats
        TYPE(dipoles), INTENT(INOUT)            :: dips
        TYPE(raman), INTENT(INOUT)           :: rams
        CHARACTER(LEN=str_len), INTENT(IN) :: input_file_name
        INTEGER :: ios, i, n, m
        REAL(dp) :: buf(10), sentinel

        CHARACTER(len=:), ALLOCATABLE :: rest
        CHARACTER(len=:), ALLOCATABLE :: s
        CHARACTER(len=str_len) :: tmp(10)
        CHARACTER(len=256) :: iomsg
        CHARACTER(len=str_len) :: dummy, line, msg
        INTEGER :: runit, stat
        INTEGER :: arr(200), nfound, group_id, frag_id
        TYPE(fragment_group), ALLOCATABLE :: tmp_group(:)
        TYPE(fragment_type), ALLOCATABLE  :: tmp_frag(:)

        !! Initialize
        LOGICAL :: in_global = .FALSE.
        LOGICAL :: in_system = .FALSE.
        LOGICAL :: in_cell = .FALSE.
        LOGICAL :: in_fragment_group = .FALSE.
        LOGICAL :: in_coordinates = .FALSE.
        LOGICAL :: in_fragments = .FALSE.
        LOGICAL :: in_md = .FALSE.
        LOGICAL :: in_static = .FALSE.
        LOGICAL :: in_hessian = .FALSE.
        LOGICAL :: in_dipoles = .FALSE.
        LOGICAL :: in_polarizabilities = .FALSE.
        LOGICAL :: in_raman = .FALSE.
        LOGICAL :: in_rtp = .FALSE.
        LOGICAL :: angles_set = .FALSE.
        sys%fragments%frag = .FALSE.
        stats%write_mol = .FALSE.
        md%write_acf = .FALSE.

        sys%fragments%ngroup = 0

        WRITE (*, '(2X, A)') "Input Data:"
        OPEN (FILE=TRIM(input_file_name), STATUS='old', ACTION='read', IOSTAT=stat, IOMSG=msg, NEWUNIT=runit)
        !Check if file exists
        CALL check_file_open(stat, msg, TRIM(input_file_name))
        !! Determine sections based on the keywords
        DO
            READ (runit, '(A)', END=100) line
            line = ADJUSTL(line)

            ! Identify section starts
            IF (INDEX(line, '&global')>0) THEN
                in_global = .TRUE.
                CYCLE
            END IF

            IF (INDEX(line, '&end global')>0) THEN
                in_global = .FALSE.
                CYCLE
            END IF

            IF (INDEX(line, '&system')>0) THEN
                in_system = .TRUE.
                CYCLE
            END IF

            IF (INDEX(line, '&end system')>0) THEN
                in_system = .FALSE.
                CYCLE
            END IF

            IF (INDEX(line, '&cell')>0) THEN
                in_cell = .TRUE.
                angles_set = .FALSE.
                CYCLE
            END IF

            IF (INDEX(line, '&end cell')>0) THEN
                in_cell = .FALSE.
                CYCLE
            END IF

            IF (INDEX(line, '&fragment')>0) THEN
                in_fragment_group = .TRUE.
                sys%fragments%frag = .TRUE.
                !! Count how many times `fragment` is triggered
                sys%fragments%ngroup = sys%fragments%ngroup + 1
                group_id = sys%fragments%ngroup

                ! allocate/expand type_frag
                IF (.NOT. ALLOCATED(sys%fragments%type_frag)) THEN
                    ALLOCATE (sys%fragments%type_frag(1))
                ELSE
                    CALL MOVE_ALLOC(sys%fragments%type_frag, tmp_group)
                    ALLOCATE (sys%fragments%type_frag(SIZE(tmp_group) + 1))
                    sys%fragments%type_frag(1:SIZE(tmp_group)) = tmp_group
                    DEALLOCATE (tmp_group)
                END IF

                sys%fragments%type_frag(group_id)%nfrag = 0
                CYCLE
            END IF

            IF (INDEX(line, '&end fragment')>0) THEN
                in_fragment_group = .FALSE.
                CYCLE
            END IF

            IF (INDEX(line, '&coordinates')>0) THEN
                in_coordinates = .TRUE.
                CYCLE
            END IF

            IF (INDEX(line, '&end coordinates')>0) THEN
                in_coordinates = .FALSE.
                CYCLE
            END IF

            IF (INDEX(line, '&md')>0) THEN
                in_md = .TRUE.
                CYCLE
            END IF

            IF (INDEX(line, '&end md')>0) THEN
                in_md = .FALSE.
                CYCLE
            END IF

            IF (INDEX(line, '&static')>0) THEN
                in_static = .TRUE.
                CYCLE
            END IF

            IF (INDEX(line, '&end static')>0) THEN
                in_static = .FALSE.
                CYCLE
            END IF

            IF (INDEX(line, '&hessian')>0) THEN
                in_hessian = .TRUE.
                CYCLE
            END IF

            IF (INDEX(line, '&end hessian')>0) THEN
                in_hessian = .FALSE.
                CYCLE
            END IF

            IF (INDEX(line, '&dipoles')>0) THEN
                in_dipoles = .TRUE.
                CYCLE
            END IF

            IF (INDEX(line, '&end dipoles')>0) THEN
                in_dipoles = .FALSE.
                CYCLE
            END IF

            IF (INDEX(line, '&polarizabilities')>0) THEN
                in_polarizabilities = .TRUE.
                CYCLE
            END IF

            IF (INDEX(line, '&end polarizabilities')>0) THEN
                in_polarizabilities = .FALSE.
                CYCLE
            END IF

            IF (INDEX(line, '&raman')>0) THEN
                in_raman = .TRUE.
                CYCLE
            END IF

            IF (INDEX(line, '&end raman')>0) THEN
                in_raman = .FALSE.
                CYCLE
            END IF

            IF (INDEX(line, '&rtp')>0) THEN
                in_rtp = .TRUE.
                CYCLE
            END IF

            IF (INDEX(line, '&end rtp')>0) THEN
                in_rtp = .FALSE.
                CYCLE
            END IF

            ! Parse within active sections
            IF (in_global) THEN
                IF (INDEX(to_lower(line), 'temperature')>0) THEN
                    READ (line, *) dummy, gs%temp
                    WRITE (*, '(4X,A, T60, F0.2)') 'Temperature (K):', gs%temp
                    !input%system%cell%present = .TRUE. ! add is present later
                ELSEIF (INDEX(to_lower(line), 'fwhm')>0) THEN
                    READ (line, *) dummy, gs%fwhm
                    WRITE (*, '(4X,A, T60, F0.4)') 'FWHM for Gaussian broadening (cm^-1):', gs%fwhm
                ELSEIF (INDEX(to_lower(line), 'spectra_verbosity')>0) THEN
                    READ (line, *) dummy, gs%spectra_verbosity
                    IF (TRIM(gs%spectra_verbosity)/='normal' .AND. TRIM(gs%spectra_verbosity)/='high') THEN
                        WRITE (error_unit, '(4X,"[WARN] ",A)') 'spectra verbosity incorrectly defined in the input.'
                        WRITE (error_unit, '(4X,"[WARN] ",A)') 'Please use keyword "normal" or "high".'
                        gs%spectra_verbosity = 'normal'
                        WRITE (error_unit, '(4X,"[WARN] ",A)') 'Spectra verbosity is now set to normal'
                    END IF
                    WRITE (*, '(4X,A, T60, A)') 'Spectra output verbosity:', gs%spectra_verbosity
                ELSEIF (INDEX(to_lower(line), 'spectra')>0) THEN
                    READ (line, *) dummy, gs%spectral_type%read_function
                    WRITE (*, '(4X,A, T60, A)') 'spectra:', TRIM(gs%spectral_type%read_function)
                END IF
            END IF

            IF (in_system) THEN
                IF (INDEX(to_lower(line), 'filename')>0) THEN
                    READ (line, *) dummy, sys%filename
                    WRITE (*, '(4X,A, T60, A)') 'Coordinates will be read from:', TRIM(sys%filename)
                ELSEIF (INDEX(to_lower(line), 'type_traj')>0) THEN
                    READ (line, *) dummy, sys%type_traj !'Enter the type of the trajectory (type pos for positions, vel for velocities)'
                    WRITE (*, '(4X,A, T60, A)') 'Type of the trajectory for the power spectrum:', TRIM(sys%type_traj)
                ELSEIF (INDEX(to_lower(line), 'mass_weighting')>0) THEN !'Do you want to apply mass weighting (y/n)?
                    READ (line, *) dummy, sys%input_mass
                    WRITE (*, '(4X,A, T60, A)') 'Mass weighting:', TRIM(sys%input_mass)
                ELSEIF (in_cell) THEN
                    IF (INDEX(to_lower(line), 'cell_type')>0) THEN !! orthorombic, hexagonal or triclinic
                        READ (line, *) dummy, sys%cell%cell_type
                        WRITE (*, '(4X,A, T60, A)') "cell type: ", sys%cell%cell_type
                    END IF
                    IF (INDEX(to_lower(line), 'box')>0) THEN
                        IF (INDEX(to_lower(line), 'box_xyz')>0) THEN
                            READ (line, *) dummy, sys%cell%box_x, sys%cell%box_y, sys%cell%box_z
                            WRITE (*, '(4X,A, T60, F0.6, F0.6, F0.6)') "cell vector: ", sys%cell%box_x, sys%cell%box_y, sys%cell%box_z
                        ELSEIF (INDEX(to_lower(line), 'box_x')>0) THEN
                            READ (line, *) dummy, sys%cell%box_x
                            WRITE (*, '(4X,A, T60, F0.6)') "cell vector x: ", sys%cell%box_x
                        ELSEIF (INDEX(to_lower(line), 'box_y')>0) THEN
                            READ (line, *) dummy, sys%cell%box_y
                            WRITE (*, '(4X,A, T60, F0.6)') "cell vector y: ", sys%cell%box_y
                        ELSEIF (INDEX(to_lower(line), 'box_z')>0) THEN
                            READ (line, *) dummy, sys%cell%box_z
                            WRITE (*, '(4X,A, T60, F0.6)') "cell vector z: ", sys%cell%box_z
                        END IF
                    END IF
                    IF (INDEX(to_lower(line), 'angle')>0) THEN
                        IF (INDEX(to_lower(line), 'angle_alpha_beta_gamma')>0) THEN
                            READ (line, *) dummy, sys%cell%angle_alpha, sys%cell%angle_beta, sys%cell%angle_gamma
                            WRITE (*, '(4X,A, T60, F0.6,1X,F0.6,1X,F0.6)') "Angle alpha: ", sys%cell%angle_alpha, sys%cell%angle_beta, sys%cell%angle_gamma
                        ELSEIF (INDEX(to_lower(line), 'angle_alpha')>0) THEN
                            READ (line, *) dummy, sys%cell%angle_alpha
                            WRITE (*, '(4X,A, T60, F0.6)') "Angle alpha: ", sys%cell%angle_alpha
                        ELSEIF (INDEX(to_lower(line), 'angle_beta')>0) THEN
                            READ (line, *) dummy, sys%cell%angle_beta
                            WRITE (*, '(4X,A, T60, F0.6)') "Angle beta: ", sys%cell%angle_beta
                        ELSEIF (INDEX(to_lower(line), 'angle_gamma')>0) THEN
                            READ (line, *) dummy, sys%cell%angle_gamma
                            WRITE (*, '(4X,A, T60, F0.6)') "Angle gamma: ", sys%cell%angle_gamma
                        END IF
                    END IF
                    IF (INDEX(to_lower(line), 'lattice_x')>0) THEN
                        ALLOCATE (sys%cell%lattice_x(3))
                        READ (line, *) dummy, sys%cell%lattice_x(1), sys%cell%lattice_x(2), sys%cell%lattice_x(3)
                        WRITE (*, '(4X,A, T60, F0.6,1X,F0.6,1X,F0.6)') "lattice vector x: ", sys%cell%lattice_x(1), sys%cell%lattice_x(2), sys%cell%lattice_x(3)
                    END IF
                    IF (INDEX(to_lower(line), 'lattice_y')>0) THEN
                        ALLOCATE (sys%cell%lattice_y(3))
                        READ (line, *) dummy, sys%cell%lattice_y(1), sys%cell%lattice_y(2), sys%cell%lattice_y(3)
                        WRITE (*, '(4X,A, T60, F0.6,1X,F0.6,1X,F0.6)') "lattice vector y: ", sys%cell%lattice_y(1), sys%cell%lattice_y(2), sys%cell%lattice_y(3)
                    END IF
                    IF (INDEX(to_lower(line), 'lattice_z')>0) THEN
                        ALLOCATE (sys%cell%lattice_z(3))
                        READ (line, *) dummy, sys%cell%lattice_z(1), sys%cell%lattice_z(2), sys%cell%lattice_z(3)
                        WRITE (*, '(4X,A, T60, F0.6,1X,F0.6,1X,F0.6)') "lattice vector z: ", sys%cell%lattice_z(1), sys%cell%lattice_z(2), sys%cell%lattice_z(3)
                    END IF
                END IF
                IF (in_fragment_group) THEN
                    IF (INDEX(to_lower(line), 'atom_list')>0) THEN
                        arr = 0
                        READ (line, *, IOSTAT=ios) dummy, (arr(i), i=1, 200)
                        nfound = COUNT(arr/=0)

                        sys%fragments%type_frag(group_id)%nfrag = &
                            sys%fragments%type_frag(group_id)%nfrag + 1
                        frag_id = sys%fragments%type_frag(group_id)%nfrag

                        ! extend fragment array inside this group
                        IF (.NOT. ALLOCATED(sys%fragments%type_frag(group_id)%fragment)) THEN
                            ALLOCATE (sys%fragments%type_frag(group_id)%fragment(1))
                        ELSE
                            CALL MOVE_ALLOC(sys%fragments%type_frag(group_id)%fragment, tmp_frag)
                            ALLOCATE (sys%fragments%type_frag(group_id)%fragment(SIZE(tmp_frag) + 1))
                            sys%fragments%type_frag(group_id)%fragment(1:SIZE(tmp_frag)) = tmp_frag
                            DEALLOCATE (tmp_frag)
                        END IF

                        ! assign atoms
                        ALLOCATE (sys%fragments%type_frag(group_id)%fragment(frag_id)%frag_atoms(nfound))
                        sys%fragments%type_frag(group_id)%fragment(frag_id)%frag_atoms = arr(1:nfound)
                    END IF
                END IF
            END IF

            IF (in_static) THEN
                IF (in_hessian) THEN
                    IF (INDEX(to_lower(line), 'force_file')>0) THEN
                        READ (line, *) dummy, stats%force_file
                        WRITE (*, '(4X,A, T60, A)') 'Forces will be read from:', TRIM(stats%force_file)
                    END IF
                END IF
                IF (INDEX(to_lower(ADJUSTL(line)), 'displacement ')==1) THEN  ! only match if first token
                    READ (line, *) dummy, stats%dx
                    WRITE (*, '(4X,A, T60, F0.6)') 'Finite difference displacement (Angstrom):', stats%dx
                END IF

                IF (INDEX(to_lower(line), 'diag_hessian')>0) THEN !Diagonalize hessian or read the normal mode freqs/disps from A file
                    READ (line, *) dummy, stats%diag_hessian
                    WRITE (*, '(4X,A, T60, A)') 'Hessian diagonalization:', TRIM(stats%diag_hessian)
                END IF
                IF (INDEX(to_lower(line), 'normal_freq_file')>0) THEN !Read normal mode frequencies
                    READ (line, *) dummy, stats%normal_freq_file
                    WRITE (*, '(4X,A, T60, A)') 'Normal mode frequencies will be read from:', TRIM(stats%normal_freq_file)
                END IF
                IF (INDEX(to_lower(line), 'normal_displ_file')>0) THEN !Read normal mode displacements
                    READ (line, *) dummy, stats%normal_displ_file
                    WRITE (*, '(4X,A, T60, A)') 'Normal mode displacements will be read from:', TRIM(stats%normal_displ_file)
                END IF
                IF (INDEX(line, 'write_mol_file')>0) THEN
                    WRITE (*, '(4X,A, T60, A)') 'A ".mol" file will be generated.'
                    stats%write_mol = .TRUE.
                END IF
            END IF

            IF (in_dipoles) THEN
                IF (INDEX(to_lower(line), 'type_dipole')>0) THEN !Type of the dipole moment
                    READ (line, *) dummy, dips%type_dipole
                    WRITE (*, '(4X,A, T60, A)') 'Type of the dipole moments:', TRIM(dips%type_dipole)
                END IF
                IF (INDEX(to_lower(line), 'dip_file')>0) THEN !Type of the dipole moment
                    READ (line, *) dummy, dips%dip_file
                    WRITE (*, '(4X,A, T60, A)') 'Dipole file:', TRIM(dips%dip_file)
                END IF
                IF (INDEX(to_lower(line), 'dip_x_file')>0) THEN !Type of the dipole moment
                    READ (line, *) dummy, dips%dip_x_file
                    WRITE (*, '(4X,A, T60, A)') 'Dipole file under x-field:', TRIM(dips%dip_x_file)
                END IF
                IF (INDEX(to_lower(line), 'dip_y_file')>0) THEN !Type of the dipole moment
                    READ (line, *) dummy, dips%dip_y_file
                    WRITE (*, '(4X,A, T60, A)') 'Dipole file under y-field:', TRIM(dips%dip_y_file)
                END IF
                IF (INDEX(to_lower(line), 'dip_z_file')>0) THEN !Type of the dipole moment
                    READ (line, *) dummy, dips%dip_z_file
                    WRITE (*, '(4X,A, T60, A)') 'Dipole file under z-field:', TRIM(dips%dip_z_file)
                END IF
            END IF

            IF (in_polarizabilities) THEN
                IF (INDEX(to_lower(line), 'type_pol')>0) THEN !Type of the polarizabilities
                    READ (line, *) dummy, rams%type_pol
                    WRITE (*, '(4X,A, T60, A)') 'Type of the polarizabilities:', TRIM(rams%type_pol)
                END IF
                IF (INDEX(to_lower(line), 'pol_x_file')>0) THEN
                    READ (line, *) dummy, rams%pol_x_file
                    WRITE (*, '(4X,A, T60, A)') 'x-components of the polarizabilities:', TRIM(rams%pol_x_file)
                END IF
                IF (INDEX(to_lower(line), 'pol_y_file')>0) THEN
                    READ (line, *) dummy, rams%pol_y_file
                    WRITE (*, '(4X,A, T60, A)') 'y-components of the polarizabilities:', TRIM(rams%pol_y_file)
                END IF
                IF (INDEX(to_lower(line), 'pol_z_file')>0) THEN
                    READ (line, *) dummy, rams%pol_z_file
                    WRITE (*, '(4X,A, T60, A)') 'z-components of the polarizabilities:', TRIM(rams%pol_z_file)
                END IF
                IF (INDEX(to_lower(line), 'static_pol_file')>0) THEN !polarizability file
                    READ (line, *) dummy, rams%static_pol_file
                    WRITE (*, '(4X,A, T60, A)') 'Polarizability file:', TRIM(rams%static_pol_file)
                END IF
                IF (INDEX(to_lower(line), 'field_strength')>0) THEN !Field strength
                    READ (line, *) dummy, dips%e_field
                    WRITE (*, '(4X,A, T60, F0.6)') 'Electric field strength (A.u.):', dips%e_field
                END IF
            END IF

            IF (in_raman) THEN
                IF (INDEX(to_lower(line), 'laser_in')>0) THEN !Type of the dipole moment
                    sentinel = HUGE(1.0_dp)
                    buf = sentinel

                    READ (line, *, iostat=ios) dummy, (buf(i), i=1, 4)
                    m = COUNT(buf/=sentinel)

                    IF (m==0) THEN
                        WRITE (error_unit, '(4X,"[Error]  ",A)') 'No laser frequencies provided after laser_in'
                        STOP
                    ELSE
                        ALLOCATE (rams%laser_in(m))
                        rams%laser_in = buf(:m)
                    END IF

                    IF (m==1) THEN
                        WRITE (*, '(4X,A, T60, F0.6)') 'Incident laser frequency (eV):', rams%laser_in
                    ELSE
                        WRITE (*, '(4X,A)') 'Found incident laser frequencies:'
                        DO i = 1, m
                            WRITE (*, '(6X,I0, A, T60, F0.6)') i, " Incident laser frequency (eV)", rams%laser_in(i)
                        END DO
                    END IF
                END IF
            END IF

            IF (in_rtp) THEN
                IF (INDEX(to_lower(line), 'rtp_time_step')>0) THEN !RTP time step
                    READ (line, *) dummy, rams%RR%dt_rtp
                    WRITE (*, '(4X,A, T60, F0.6)') 'RTP time step (fs):', rams%RR%dt_rtp
                END IF
                IF (INDEX(to_lower(line), 'rtp_framecount')>0) THEN !RTP time step
                    READ (line, *) dummy, rams%RR%framecount_rtp
                    WRITE (*, '(4X,A, T60,I0)') "RTP framecount: ", rams%RR%framecount_rtp
                END IF
                IF (INDEX(to_lower(line), 'check_pade')>0) THEN !RTP time step
                    READ (line, *) dummy, rams%RR%check_pade
                    WRITE (*, '(4X,A, T60, A)') 'Apply Pade:', TRIM(rams%RR%check_pade)
                END IF
                IF (INDEX(to_lower(line), 'pade_framecount')>0) THEN !RTP time step
                    READ (line, *) dummy, rams%RR%framecount_rtp_pade
                    WRITE (*, '(4X,A, T60,I0)') "Requested framecount after Pade: ", rams%RR%framecount_rtp_pade
                END IF
                IF (INDEX(to_lower(line), 'damping_constant')>0) THEN !RTP time step
                    READ (line, *) dummy, rams%RR%damping_constant
                    WRITE (*, '(4X,A, T60, F0.6)') 'Damping constant (eV):', rams%RR%damping_constant
                END IF
            END IF

            IF (in_md) THEN
                IF (INDEX(line, 'time_step')>0) THEN
                    READ (line, *) dummy, md%dt
                    WRITE (*, '(4X,A, T60, F0.6)') 'MD time step (fs):', md%dt
                END IF
                IF (INDEX(line, 'correlation_depth')>0) THEN
                    READ (line, *) dummy, md%t_cor
                    WRITE (*, '(4X,A, T60, I0)') 'Correlation depth:', md%t_cor
                END IF
                IF (INDEX(line, 'write_acf_file')>0) THEN
                    WRITE (*, '(4X,A, T60, A)') 'Output file that contains the autocorrelation data will be printed out.'
                    md%write_acf = .TRUE.
                END IF
            END IF
        END DO
100     CONTINUE
        CLOSE (runit)

    END SUBROUTINE parse_input

!*********************************************************************************************************!
!*********************************************************************************************************!
    !> @brief Verifies that all required input parameters are defined and consistent.
    !>
    !> This routine performs error and sanity checks on all input data after parsing.
    !> It ensures that mandatory parameters (filenames, time steps, displacements,
    !> field strengths, etc.) are provided for each spectral type and if possible, assigns
    !> default values to unidentified parameters.
    !>
    !> @param[in,out] gs     Global settings (provides spectral type, temperature, broadening, etc.).
    !> @param[in,out] sys    System information (provides trajectory, coordinates, mass weighting).
    !> @param[in,out] md     Molecular dynamics information (provides time step, correlation depth).
    !> @param[in,out] stats  Static spectral information (Hessian, displacement, etc.).
    !> @param[in,out] dips   Dipole-related information (files, dipole type, field strength).
    !> @param[in,out] rams   Raman and resonance Raman settings (laser frequencies, damping, etc.).
    !>
    SUBROUTINE check_input(gs, sys, md, stats, dips, rams)

        TYPE(global_settings)       :: gs
        TYPE(systems)               :: sys
        TYPE(molecular_dynamics)    :: md
        TYPE(static)                :: stats
        TYPE(dipoles)               :: dips
        TYPE(raman)                 :: rams

        IF (TRIM(gs%spectral_type%read_function)=='') THEN
            WRITE (error_unit, '(4X,"[ERROR] ",A)') 'Spectra not defined in the input'
            STOP
            !****************************************************************************************************!
            !****************************************************************************************************!
        !!Check for the power spectrum
        ELSEIF (gs%spectral_type%read_function=='P') THEN
            !check for input_type
            IF (TRIM(sys%type_traj)=='') THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'type_traj not defined in the input - provide "type_traj pos" for positions, &
        &                  "type_traj vel" for velocities'
                STOP
            END IF
            !check for filename
            IF (TRIM(sys%filename)=='') THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'Filename not defined in the input'
                STOP
            END IF
            !check for mass_weighting
            IF (TRIM(sys%input_mass)=='') THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'mass_weighting not defined in the input'
                STOP
            END IF
            !check time step
            IF (md%dt<0) THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'time_step not defined in the input'
                STOP
            END IF
            !check t_cor
            IF (md%t_cor<0) THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'correlation depth not defined in the input, we will continue with an estimate' !can be worded differently
                STOP
            END IF
            !check for the temperature
            IF (gs%temp<0) THEN
                WRITE (error_unit, '(4X,"[WARN]  ",A)') 'Temperature is not defined, setting it to 300 K'
                gs%temp = 300
            END IF
            !****************************************************************************************************!
            !****************************************************************************************************!
            !Check for normal mode analysis
        ELSEIF (gs%spectral_type%read_function=='NMA') THEN
            !check for input_dipole not needed for P but set to A default value
      !      IF (TRIM(dips%type_dipole)=='') THEN
      !          WRITE (error_unit, '(4X,"[WARN]  ",A)') 'type_dipole not defined in the input setting it to 1'
      !          dips%type_dipole = '1'
      !      END IF
            !check for filename
            IF (TRIM(sys%filename)=='') THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'Filename not defined in the input'
                STOP
            END IF
            !check for force_file
            IF (TRIM(stats%force_file)=='') THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'Force filename not defined in the input'
                STOP
            END IF
            !check for displacement
            IF (stats%dx<0) THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'Displacement not defined in the input'
                STOP
            END IF
            !****************************************************************************************************!
            !****************************************************************************************************!
            !Check for static IR
        ELSEIF (gs%spectral_type%read_function=='IR') THEN
            !check for filename
            IF (TRIM(sys%filename)=='') THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'Filename not defined in the input'
                STOP
            END IF
            !check for force_file
            IF (stats%diag_hessian=='y') THEN
                IF (TRIM(stats%force_file)=='') THEN
                    WRITE (error_unit, '(4X,"[ERROR] ",A)') 'File name of the forces not defined in the input'
                    STOP
                END IF
            ELSEIF (stats%diag_hessian=='n') THEN
                !check for the external file containing the normal mode frequencies
                IF (TRIM(stats%normal_freq_file)=='') THEN
                    WRITE (error_unit, '(4X,"[ERROR] ",A)') 'File name of the normal mode frequencies not defined in the input'
                    STOP
                END IF
                !check for the external file containing the normal mode displacements
                IF (TRIM(stats%normal_displ_file)=='') THEN
                    WRITE (error_unit, '(4X,"[ERROR] ",A)') 'File name of the normal mode displacements not defined in the input'
                    STOP
                END IF
            END IF
            !check for normal mode displacement
            IF (stats%dx<0) THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'Displacement not defined in the input'
                STOP
            END IF
            !check for dipole file
            IF (TRIM(dips%dip_file)=='') THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'Dipole filename not defined in the input'
                STOP
            END IF
            !check for type_dipole
            IF (TRIM(dips%type_dipole)=='') THEN
                WRITE (error_unit, '(4X,"[WARN]  ",A)') 'type_dipole not defined in the input, setting it to berry'
                dips%type_dipole = 'berry'
            END IF
            !check for the full width at half-maximum
            IF (gs%fwhm<0) THEN
                WRITE (error_unit, '(4X,"[WARN]  ",A)') 'Full width at half-maximum is not defined, setting it to 10 cm^{-1}'
                gs%fwhm = 10
            END IF

            !****************************************************************************************************!
            !****************************************************************************************************!
            !Check for static Raman
        ELSEIF (gs%spectral_type%read_function=='R') THEN
            !check for filename
            IF (TRIM(sys%filename)=='') THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'Filename not defined in the input'
                STOP
            END IF
            !check for force_file
            IF (stats%diag_hessian=='y') THEN
                IF (TRIM(stats%force_file)=='') THEN
                    WRITE (error_unit, '(4X,"[ERROR] ",A)') 'File name of the forces not defined in the input'
                    STOP
                END IF
            ELSEIF (stats%diag_hessian=='n') THEN
                !check for the external file containing the normal mode frequencies
                IF (TRIM(stats%normal_freq_file)=='') THEN
                    WRITE (error_unit, '(4X,"[ERROR] ",A)') 'File name of the normal mode frequencies not defined in the input'
                    STOP
                END IF
                !check for the external file containing the normal mode displacements
                IF (TRIM(stats%normal_displ_file)=='') THEN
                    WRITE (error_unit, '(4X,"[ERROR] ",A)') 'File name of the normal mode displacements not defined in the input'
                    STOP
                END IF
            END IF
            !check for normal mode displacement
            IF (stats%dx<0) THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'Displacement not defined in the input'
                STOP
            END IF
            !check for polarizability file
            IF (TRIM(rams%static_pol_file)=='') THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'Polarizability filename not defined in the input'
                STOP
            END IF
            !check for type_pol
            IF (TRIM(rams%type_pol)=='') THEN
                WRITE (error_unit, '(4X,"[WARN]  ",A)') 'type_pol not defined in the input, setting it to "analytical"'
                rams%type_pol = 'analytical'
            END IF
            !check for the electric field strength
            IF (rams%type_pol.NE.'analytical' .AND. dips%e_field<0) THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'Electric field strength not defined!'
                STOP
            END IF
            !check for the incident laser frequency
            IF (.NOT. ALLOCATED(rams%laser_in)) THEN
                WRITE (error_unit, '(4X,"[WARN]  ",A)') 'Incident laser frequency not defined, setting it to 0.5 eV'
                ALLOCATE (rams%laser_in(1))
                rams%laser_in(1) = 0.5
            END IF
            !check for the temperature
            IF (gs%temp<0) THEN
                WRITE (error_unit, '(4X,"[WARN]  ",A)') 'Temperature is not defined, setting it to 300 K'
                gs%temp = 300
            END IF
            !check for the full width at half-maximum
            IF (gs%fwhm<0) THEN
                WRITE (error_unit, '(4X,"[WARN]  ",A)') 'Full width at half-maximum is not defined, setting it to 10 cm^{-1}'
                gs%fwhm = 10
            END IF

            !****************************************************************************************************!
            !****************************************************************************************************!
            !Check for absorption spectrum
        ELSEIF (gs%spectral_type%read_function=='ABS') THEN
            !check for filename
            IF (TRIM(sys%filename)=='') THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'Filename not defined in the input'
                STOP
            END IF
            !check for the file containing perturbed dipole moments under x_field
            IF (TRIM(dips%dip_x_file)=='') THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'X-field dipole file name not defined in the input'
                STOP
            END IF
            !check for the file containing perturbed dipole moments under y_field
            IF (TRIM(dips%dip_y_file)=='') THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'Y-field dipole file name not defined in the input'
                STOP
            END IF
            !check for the file containing perturbed dipole moments under z_field
            IF (TRIM(dips%dip_z_file)=='') THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'Z-field dipole file name not defined in the input'
                STOP
            END IF
            !check for type_pol
            IF (TRIM(rams%type_pol)=='') THEN
                WRITE (error_unit, '(4X,"[WARN]  ",A)') 'type_pol not defined in the input, setting it to "induced"'
                rams%type_pol = 'induced'
            END IF
            !check for type_dipole
            IF (TRIM(dips%type_dipole)=='') THEN
                WRITE (error_unit, '(4X,"[WARN]  ",A)') 'type_dipole not defined in the input setting it to "berry"'
                dips%type_dipole = 'berry'
            END IF
            !check for RT-TDDFT time step
            IF (rams%RR%dt_rtp<0) THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'RTP time step not defined!'
                STOP
            END IF
            !check for the electric field strength
            IF (rams%type_pol.NE.'analytical' .AND. dips%e_field<0) THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'Electric field strength not defined!'
                STOP
            END IF
            !check for the RT-TDDFT frame count
            IF (rams%RR%framecount_rtp<0) THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'RTP framecount not defined!'
                STOP
            END IF
            !check for the damping constant
            IF (rams%RR%damping_constant<0) THEN
                WRITE (error_unit, '(4X,"[WARN]  ",A)') 'Damping constant not defined, setting it to 0.1 eV!'
                rams%RR%damping_constant = 0.1_dp
            END IF
            !check if Pade interpolation is requested
            IF (TRIM(rams%RR%check_pade)=='') THEN
                WRITE (error_unit, '(4X,"[WARN]  ",A)') 'The calculation will continue without Pade approximants!'
                rams%RR%check_pade = 'n'
            END IF
            !check for the number of frames after the Pade interpolation
            IF (TRIM(rams%RR%check_pade)=='y' .AND. rams%RR%framecount_rtp_pade<0) THEN !this can also be adjusted
                WRITE (error_unit, '(4X,"[WARN]  ",A)') 'Pade framecount is set to 80000!'
                rams%RR%framecount_rtp_pade = 80000
            END IF
            !****************************************************************************************************!
            !****************************************************************************************************!
            !Check for static resonance Raman
        ELSEIF (gs%spectral_type%read_function=='RR') THEN
            !check for filename
            IF (TRIM(sys%filename)=='') THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'Filename not defined in the input'
                STOP
            END IF
            !check for force_file
            IF (stats%diag_hessian=='y') THEN
                IF (TRIM(stats%force_file)=='') THEN
                    WRITE (error_unit, '(4X,"[ERROR] ",A)') 'File name of the forces not defined in the input'
                    STOP
                END IF
            ELSEIF (stats%diag_hessian=='n') THEN
                !check for the external file containing the normal mode frequencies
                IF (TRIM(stats%normal_freq_file)=='') THEN
                    WRITE (error_unit, '(4X,"[ERROR] ",A)') 'File name of the normal mode frequencies not defined in the input'
                    STOP
                END IF
                !check for the external file containing the normal mode displacements
                IF (TRIM(stats%normal_displ_file)=='') THEN
                    WRITE (error_unit, '(4X,"[ERROR] ",A)') 'File name of the normal mode displacements not defined in the input'
                    STOP
                END IF
            END IF
            !check for normal mode displacement
            IF (stats%dx<0) THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'Displacement not defined in the input'
                STOP
            END IF
            !check for type_pol
            IF (TRIM(rams%type_pol)=='') THEN
                WRITE (error_unit, '(4X,"[WARN]  ",A)') 'type_pol not defined in the input, setting it to "induced"'
                rams%type_pol = 'induced'
            END IF
            !check for type_dipole
            IF (TRIM(dips%type_dipole)=='') THEN
                WRITE (error_unit, '(4X,"[WARN]  ",A)') 'type_dipole not defined in the input setting it to "berry"'
                dips%type_dipole = 'berry'
            END IF
            !check for the file containing perturbed dipole moments under x_field
            IF (TRIM(dips%dip_x_file)=='') THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'X-field dipole file name not defined in the input'
                STOP
            END IF
            !check for the file containing perturbed dipole moments under y_field
            IF (TRIM(dips%dip_y_file)=='') THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'Y-field dipole file name not defined in the input'
                STOP
            END IF
            !check for the file containing perturbed dipole moments under z_field
            IF (TRIM(dips%dip_z_file)=='') THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'Z-field dipole file name not defined in the input'
                STOP
            END IF
            !check for RT-TDDFT time step
            IF (rams%RR%dt_rtp<0) THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'RTP time step not defined!'
                STOP
            END IF
            !check for the electric field strength
            IF (rams%type_pol.NE.'analytical' .AND. dips%e_field<0) THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'Electric field strength not defined!'
                STOP
            END IF
            !check for RT-TDDFT frame count
            IF (rams%RR%framecount_rtp<0) THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'RTP framecount not defined!'
                STOP
            END IF
            !check if Pade interpolation is requested
            IF (TRIM(rams%RR%check_pade)=='') THEN
                WRITE (error_unit, '(4X,"[WARN]  ",A)') 'The calculation will continue without Pade approximants!'
                rams%RR%check_pade = 'n'
            END IF
            !check for the number of frames after the Pade interpolation
            IF (TRIM(rams%RR%check_pade)=='y' .AND. rams%RR%framecount_rtp_pade<0) THEN !this can also be adjusted
                WRITE (error_unit, '(4X,"[WARN]  ",A)') 'Pade framecount is set to 80000!'
                rams%RR%framecount_rtp_pade = 80000
            END IF
            !check for the damping constant
            IF (rams%RR%damping_constant<0) THEN
                WRITE (error_unit, '(4X,"[WARN]  ",A)') 'Damping constant not defined, setting it to 0.1 eV!'
                rams%RR%damping_constant = 0.1_dp
            END IF
            !check for the incident laser frequency
            IF (.NOT. ALLOCATED(rams%laser_in)) THEN
                WRITE (error_unit, '(4X,"[WARN]  ",A)') 'Incident laser frequency not defined, setting it to 0.5 eV'
                ALLOCATE (rams%laser_in(1))
                IF (.NOT. ALLOCATED(rams%laser_in)) THEN
                    WRITE (error_unit, '(4X,"[WARN]  ",A)') 'Incident laser frequency not defined, setting it to 0.5 eV'
                    ALLOCATE (rams%laser_in(1))
                    rams%laser_in(1) = 0.5
                END IF
            END IF
            !check for the temperature
            IF (gs%temp<0) THEN
                WRITE (error_unit, '(4X,"[WARN]  ",A)') 'Temperature is not defined, setting it to 300 K'
                gs%temp = 300
            END IF
            !check for the full width at half-maximum
            IF (gs%fwhm<0) THEN
                WRITE (error_unit, '(4X,"[WARN]  ",A)') 'Full width at half-maximum is not defined, setting it to 10 cm^{-1}'
                gs%fwhm = 10
            END IF
            !****************************************************************************************************!
            !****************************************************************************************************!
            !Check for MD-based IR
        ELSEIF (gs%spectral_type%read_function=='MD-IR') THEN
            !check for type_dipole
            IF (TRIM(dips%type_dipole)=='') THEN
                WRITE (error_unit, '(4X,"[WARN]  ",A)') 'type_dipole not defined in the input setting it to "berry"'
                dips%type_dipole = 'berry'
            END IF
            !check for dipole file
            IF (TRIM(dips%dip_file)=='') THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'Dipole filename not defined in the input'
                STOP
            END IF
            !check MD time step
            IF (md%dt<0) THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'time_step not defined in the input'
                STOP
            END IF
            !check correlation depth
            IF (md%t_cor<0) THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'correlation depth not defined in the input, we will continue with an estimate' !can be worded differently
                STOP
            END IF

            CALL check_lattice_parameters(sys)
            !****************************************************************************************************!
            !****************************************************************************************************!
            !Check for MD-based Raman
        ELSEIF (gs%spectral_type%read_function=='MD-R') THEN
            !check for type_pol
            IF (TRIM(rams%type_pol)=='') THEN
                WRITE (error_unit, '(4X,"[WARN]  ",A)') 'type_pol not defined in the input, please set it to "induced" or "analytical"'
                STOP
            END IF
            !check for type_dipole
            IF (TRIM(dips%type_dipole)=='') THEN
                WRITE (error_unit, '(4X,"[WARN]  ",A)') 'type_dipole not defined in the input setting it to "berry"'
                dips%type_dipole = 'berry'
            END IF
            !check for electric field strength
            IF (rams%type_pol.NE.'analytical' .AND. dips%e_field<0) THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'Electric field strength not defined!'
                STOP
            END IF
            !check for dipole file
            IF (TRIM(dips%dip_file)=='') THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'Dipole filename not defined in the input'
                STOP
            END IF
            !check for the file containing perturbed dipole moments under x_field
            IF (TRIM(dips%dip_x_file)=='') THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'X-field dipole file name not defined in the input'
                STOP
            END IF
            !check for the file containing perturbed dipole moments under y_field
            IF (TRIM(dips%dip_y_file)=='') THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'Y-field dipole file name not defined in the input'
                STOP
            END IF
            !check for the file containing perturbed dipole moments under z_field
            IF (TRIM(dips%dip_z_file)=='') THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'Z-field dipole file name not defined in the input'
                STOP
            END IF
            !check MD time step
            IF (md%dt<0) THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'time_step not defined in the input'
                STOP
            END IF
            !check correlation depth
            IF (md%t_cor<0) THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'correlation depth not defined in the input, we will continue with an estimate' !can be worded differently
                STOP
            END IF
            !check for the temperature
            IF (gs%temp<0) THEN
                WRITE (error_unit, '(4X,"[WARN]  ",A)') 'Temperature is not defined, setting it to 300 K'
                gs%temp = 300
            END IF
            !check for incident laser wavelength
            IF (.NOT. ALLOCATED(rams%laser_in)) THEN
                WRITE (error_unit, '(4X,"[WARN]  ",A)') 'Incident laser frequency not defined, setting it to 0.5 eV'
                ALLOCATE (rams%laser_in(1))
                rams%laser_in(1) = 0.5
            END IF
            IF (TRIM(dips%type_dipole)=='wannier') THEN
                CALL check_lattice_parameters(sys)
            END IF
            !****************************************************************************************************!
            !****************************************************************************************************!
            !Check for MD-based absorption spectrum
        ELSEIF (gs%spectral_type%read_function=='MD-ABS') THEN
            !check for type_pol
            IF (TRIM(rams%type_pol)=='') THEN
                WRITE (error_unit, '(4X,"[WARN]  ",A)') 'type_pol not defined in the input, setting it to "induced"'
                rams%type_pol = 'induced'
            END IF
            !check for type_dipole
            IF (TRIM(dips%type_dipole)=='') THEN
                WRITE (error_unit, '(4X,"[WARN]  ",A)') 'type_dipole not defined in the input setting it to berry'
                dips%type_dipole = 'berry'
            END IF
            !check for electric field strength
            IF (rams%type_pol.NE.'analytical' .AND. dips%e_field<0) THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'Electric field strength not defined!'
                STOP
            END IF
            !check for the file containing perturbed dipole moments under x_field
            IF (TRIM(dips%dip_x_file)=='') THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'X-field dipole file name not defined in the input'
                STOP
            END IF
            !check for the file containing perturbed dipole moments under y_field
            IF (TRIM(dips%dip_y_file)=='') THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'Y-field dipole file name not defined in the input'
                STOP
            END IF
            !check for the file containing perturbed dipole moments under z_field
            IF (TRIM(dips%dip_z_file)=='') THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'Z-field dipole file name not defined in the input'
                STOP
            END IF
            !check for RT-TDDFT time step
            IF (rams%RR%dt_rtp<0) THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'RTP time step not defined!'
                STOP
            END IF
            !check for RT-TDDFT frame count
            IF (rams%RR%framecount_rtp<0) THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'RTP framecount not defined!'
                STOP
            END IF
            !check for the damping constant
            IF (rams%RR%damping_constant<0) THEN
                WRITE (error_unit, '(4X,"[WARN]  ",A)') 'Damping constant not defined, setting it to 0.1 eV!'
                rams%RR%damping_constant = 0.1_dp
            END IF
            !check if Pade interpolation is requested
            IF (TRIM(rams%RR%check_pade)=='') THEN
                WRITE (error_unit, '(4X,"[WARN]  ",A)') 'The calculation will continue without Pade approximants!'
                rams%RR%check_pade = 'n'
            END IF
            !check for the number of frames after the Pade interpolation
            IF (TRIM(rams%RR%check_pade)=='y' .AND. rams%RR%framecount_rtp_pade<0) THEN !this can also be adjusted
                WRITE (error_unit, '(4X,"[WARN]  ",A)') 'Pade framecount is set to 80000!'
                rams%RR%framecount_rtp_pade = 20000
            END IF
            !****************************************************************************************************!
            !****************************************************************************************************!
            !Check for MD-based resonance Raman spectrum
        ELSEIF (gs%spectral_type%read_function=='MD-RR') THEN
            !check for type_pol
            IF (TRIM(rams%type_pol)=='') THEN
                WRITE (error_unit, '(4X,"[WARN]  ",A)') 'type_pol not defined in the input, setting it to "induced"'
                rams%type_pol = 'induced'
            END IF
            !check for type_dipole
            IF (TRIM(dips%type_dipole)=='') THEN
                WRITE (error_unit, '(4X,"[WARN]  ",A)') 'type_dipole not defined in the input setting it to berry'
                dips%type_dipole = 'berry'
            END IF
            !check for electric field strength
            IF (rams%type_pol.NE.'analytical' .AND. dips%e_field<0) THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'Electric field strength not defined!'
                STOP
            END IF
            IF (rams%type_pol=='induced') THEN
                !check for the file containing perturbed dipole moments under x_field
                IF (TRIM(dips%dip_x_file)=='') THEN
                    WRITE (error_unit, '(4X,"[ERROR] ",A)') 'X-field dipole file name not defined in the input'
                    STOP
                END IF
                !check for the file containing perturbed dipole moments under y_field
                IF (TRIM(dips%dip_y_file)=='') THEN
                    WRITE (error_unit, '(4X,"[ERROR] ",A)') 'Y-field dipole file name not defined in the input'
                    STOP
                END IF
                !check for the file containing perturbed dipole moments under z_field
                IF (TRIM(dips%dip_z_file)=='') THEN
                    WRITE (error_unit, '(4X,"[ERROR] ",A)') 'Z-field dipole file name not defined in the input'
                    STOP
                END IF
            ELSEIF (rams%type_pol=='analytical') THEN
                !check for the file containing x-components of the polarizabilities
                IF (TRIM(rams%pol_x_file)=='') THEN
                    WRITE (error_unit, '(4X,"[ERROR] ",A)') 'file name for the x-component of polarizabilities not defined in the input '
                    STOP
                END IF
                !check for the file containing y-components of the polarizabilities
                IF (TRIM(rams%pol_y_file)=='') THEN
                    WRITE (error_unit, '(4X,"[ERROR] ",A)') 'file name for the y-component of polarizabilities not defined in the input '
                    STOP
                END IF
                !check for the file containing z-components of the polarizabilities
                IF (TRIM(rams%pol_z_file)=='') THEN
                    WRITE (error_unit, '(4X,"[ERROR] ",A)') 'file name for the z-component of polarizabilities not defined in the input '
                    STOP
                END IF

            END IF
            !check MD time step
            IF (md%dt<0) THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'time_step not defined in the input'
                STOP
            END IF
            !check correlation depth
            IF (md%t_cor<0) THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'correlation depth not defined in the input, we will continue with an estimate' !can be worded differently
                STOP
            END IF
            !check for the temperature
            IF (gs%temp<0) THEN
                WRITE (error_unit, '(4X,"[WARN]  ",A)') 'Temperature is not defined, setting it to 300 K'
                gs%temp = 300
            END IF
            !check for incident laser wavelength
            IF (.NOT. ALLOCATED(rams%laser_in)) THEN
                WRITE (error_unit, '(4X,"[WARN]  ",A)') 'Incident laser frequency not defined, setting it to 0.5 eV'
                ALLOCATE (rams%laser_in(1))
                rams%laser_in(1) = 0.5
            END IF
            !check for RT-TDDFT time step
            IF (rams%RR%dt_rtp<0) THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'RTP time step not defined!'
                STOP
            END IF
            !check for RT-TDDFT frame count
            IF (rams%RR%framecount_rtp<0) THEN
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'RTP framecount not defined!'
                STOP
            END IF
            !check for the damping constant
            IF (rams%RR%damping_constant<0) THEN
                WRITE (error_unit, '(4X,"[WARN]  ",A)') 'Damping constant not defined, setting it to 0.1 eV!'
                rams%RR%damping_constant = 0.1_dp
            END IF
            !check if Pade interpolation is requested
            IF (TRIM(rams%RR%check_pade)=='') THEN
                WRITE (error_unit, '(4X,"[WARN]  ",A)') 'The calculation will continue without Pade approximants!'
                rams%RR%check_pade = 'n'
            END IF
            !check for the number of frames after the Pade interpolation
            IF (TRIM(rams%RR%check_pade)=='y' .AND. rams%RR%framecount_rtp_pade<0) THEN !this can also be adjusted
                WRITE (error_unit, '(4X,"[WARN]  ",A)') 'Pade framecount is set to 80000!'
                rams%RR%framecount_rtp_pade = 20000
            END IF
        END IF

    END SUBROUTINE check_input

END MODULE read_input
