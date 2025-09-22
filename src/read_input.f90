MODULE read_input

    USE kinds, ONLY: sp, dp, str_len
    USE vib_types, ONLY: global_settings, systems, static, dipoles, raman, molecular_dynamics

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: parse_command_line, parse_input, check_input

CONTAINS

    !****************************************************************
    ! doxygen doc, to be added
    !****************************************************************
    SUBROUTINE parse_command_line(input_file_name)

        CHARACTER(LEN=str_len), INTENT(OUT) :: input_file_name

        INTEGER(KIND=sp)                                     :: narg

        narg = iargc()

        IF (narg/=1) THEN
            WRITE (*, '(A37)') "Usage: vibrant_input.x your_input.inp"
            STOP
        END IF

        CALL getarg(1, input_file_name)

    END SUBROUTINE parse_command_line

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

    !****************************************************************
    ! doxygen doc, to be added
    !****************************************************************
    SUBROUTINE parse_input(gs, sys, md, stats, dips, rams, input_file_name)

        TYPE(global_settings), INTENT(INOUT) :: gs
        TYPE(systems), INTENT(INOUT)        :: sys
        TYPE(molecular_dynamics), INTENT(INOUT)     :: md
        TYPE(static), INTENT(INOUT)           :: stats
        TYPE(dipoles), INTENT(INOUT)            :: dips
        TYPE(raman), INTENT(INOUT)           :: rams
        CHARACTER(LEN=str_len), INTENT(IN) :: input_file_name
        !INTEGER :: ios
        !CHARACTER(len=256) :: iomsg
        !** intermal variables
        CHARACTER(LEN=str_len) :: line
        CHARACTER(LEN=str_len) :: dummy
        LOGICAL :: in_global = .FALSE.
        LOGICAL :: in_system = .FALSE.
        LOGICAL :: in_cell = .FALSE.
        LOGICAL :: in_coordinates = .FALSE.
        LOGICAL :: in_fragments = .FALSE.
        LOGICAL :: in_md = .FALSE.
        LOGICAL :: in_static = .FALSE.
        LOGICAL :: in_hessian = .FALSE.
        LOGICAL :: in_dipoles = .FALSE.
        LOGICAL :: in_raman = .FALSE.
        LOGICAL :: in_rtp = .FALSE.

        OPEN (unit=999, file=TRIM(input_file_name), status="old")

        DO
            READ (999, '(A)', END=100) line
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
                CYCLE
            END IF

            IF (INDEX(line, '&end cell')>0) THEN
                in_cell = .FALSE.
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

            IF (INDEX(line, '&fragments')>0) THEN
                in_fragments = .TRUE.
                CYCLE
            END IF

            IF (INDEX(line, '&end fragments')>0) THEN
                in_fragments = .FALSE.
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
                    WRITE (*, *) "Temperature (K): ", gs%temp
                    !input%system%cell%present = .TRUE. ! add is present later
                ELSEIF (INDEX(to_lower(line), 'fwhm')>0) THEN
                    READ (line, *) dummy, gs%fwhm
                    WRITE (*, *) "Full width at half-maximum for Gaussian broadening (cm^{-1}): ", gs%fwhm
                ELSEIF (INDEX(to_lower(line), 'spectra')>0) THEN
                    READ (line, *) dummy, gs%spectral_type%read_function
                    WRITE (*, *) "spectra: ", gs%spectral_type%read_function
                END IF
            END IF

            IF (in_system) THEN
                IF (INDEX(to_lower(line), 'filename')>0) THEN
                    READ (line, *) dummy, sys%filename
                    WRITE (*, *) "Coordinates will be read from: ", sys%filename
                ELSEIF (INDEX(to_lower(line), 'type_traj')>0) THEN
                    READ (line, *) dummy, sys%type_traj !'Enter the type of the trajectory (type pos for positions, vel for velocities)'
                    WRITE (*, *) "Type of the trajectory for the power spectrum: ", sys%type_traj
                ELSEIF (INDEX(to_lower(line), 'mass_weighting')>0) THEN !'Do you want to apply mass weighting (y/n)?
                    READ (line, *) dummy, sys%input_mass
                    WRITE (*, *) "Mass weighting: ", sys%input_mass
                ELSEIF (in_cell) THEN
                    IF (INDEX(to_lower(line), 'cell_type')>0) THEN !'Is it the k-point trajectory (1), supercell trajectory (2) or solvent trajectory (3)? '
                        READ (line, *) dummy, sys%cell%cell_type
                        !    ELSEIF (INDEX(to_lower(line), 'abc')>0) THEN
                        !        !READ (line, *) dummy, sys%cell%abc(1:3)
                        !        !input%system%cell%present = .TRUE.
                        !    ELSEIF (INDEX(to_lower(line), 'alpha_beta_gamma')>0) THEN
                        !        !READ (line, *) dummy, input%system%cell%alpha_beta_gamma(1:3)
                    END IF
                ELSEIF (INDEX(to_lower(line), 'frag_type ')>0) THEN !'Does the system contain more than one molecule? (y/n)'
                    READ (line, *) dummy, sys%frag_type
                    WRITE (*, *) "frag_type: ", sys%frag_type
                END IF
            END IF

            IF (in_static) THEN
                IF (in_hessian) THEN
                    IF (INDEX(to_lower(line), 'force_file')>0) THEN
                        READ (line, *) dummy, stats%force_file
                        WRITE (*, *) "Forces will be read from: ", stats%force_file
                    END IF
                END IF
                IF (INDEX(to_lower(ADJUSTL(line)), 'displacement ')==1) THEN  ! only match if first token
                    READ (line, *) dummy, stats%dx
                    WRITE (*, *) "Finite difference displacement (Angstrom): ", stats%dx
                END IF
                !READ (line, *, IOSTAT=ios, IOMSG=iomsg) dummy, stats%dx
                !IF (ios /= 0) THEN
                !   PRINT *, "Error reading displacement line: ", TRIM(iomsg)
                !   PRINT *, "Line content was: ", TRIM(line)
                !   STOP
                !END IF

                IF (INDEX(to_lower(line), 'diag_hessian')>0) THEN !Diagonalize hessian or read the normal mode freqs/disps from a file
                    READ (line, *) dummy, stats%diag_hessian
                    WRITE (*, *) "Hessian diagonalization: ", stats%diag_hessian
                END IF
                IF (INDEX(to_lower(line), 'normal_freq_file')>0) THEN !Read normal mode frequencies
                    READ (line, *) dummy, stats%normal_freq_file
                    WRITE (*, *) "Normal mode frequencies will be read from: ", stats%normal_freq_file
                END IF
                IF (INDEX(to_lower(line), 'normal_displ_file')>0) THEN !Read normal mode displacements
                    READ (line, *) dummy, stats%normal_displ_file
                    WRITE (*, *) "Normal mode displacements will be read from: ", stats%normal_displ_file
                END IF

            END IF

            IF (in_dipoles) THEN
                IF (INDEX(to_lower(line), 'type_dipole')>0) THEN !Type of the dipole moment
                    READ (line, *) dummy, dips%type_dipole
                    WRITE (*, *) "Type of the dipole moments: ", dips%type_dipole
                END IF
                IF (INDEX(to_lower(line), 'dip_file')>0) THEN !Type of the dipole moment
                    READ (line, *) dummy, dips%dip_file
                    WRITE (*, *) "Dipole file: ", dips%dip_file
                END IF
                IF (INDEX(to_lower(line), 'dip_x_file')>0) THEN !Type of the dipole moment
                    READ (line, *) dummy, dips%dip_x_file
                    WRITE (*, *) "Dipole file under x-field: ", dips%dip_x_file
                END IF
                IF (INDEX(to_lower(line), 'dip_y_file')>0) THEN !Type of the dipole moment
                    READ (line, *) dummy, dips%dip_y_file
                    WRITE (*, *) "Dipole file under y-field: ", dips%dip_y_file
                END IF
                IF (INDEX(to_lower(line), 'dip_z_file')>0) THEN !Type of the dipole moment
                    READ (line, *) dummy, dips%dip_z_file
                    WRITE (*, *) "Dipole file under z-field: ", dips%dip_z_file
                END IF
                IF (INDEX(to_lower(line), 'static_pol_file')>0) THEN !Type of the dipole moment
                    READ (line, *) dummy, rams%static_pol_file
                    WRITE (*, *) "Polarizability file: ", rams%static_pol_file
                END IF
                IF (INDEX(to_lower(line), 'field_strength')>0) THEN !Field strength
                    READ (line, *) dummy, dips%e_field
                    WRITE (*, *) "Electric field strength (a.u.): ", dips%e_field
                END IF
            END IF

            IF (in_raman) THEN
                IF (INDEX(to_lower(line), 'laser_in')>0) THEN !Type of the dipole moment
                    READ (line, *) dummy, rams%laser_in
                    WRITE (*, *) "Incident laser frequency (eV): ", rams%laser_in
                END IF
            END IF

            IF (in_rtp) THEN
                IF (INDEX(to_lower(line), 'rtp_time_step')>0) THEN !RTP time step
                    READ (line, *) dummy, rams%RR%dt_rtp
                    WRITE (*, *) "RTP time step (fs): ", rams%RR%dt_rtp
                END IF
                IF (INDEX(to_lower(line), 'rtp_framecount')>0) THEN !RTP time step
                    READ (line, *) dummy, rams%RR%framecount_rtp
                    WRITE (*, *) "RTP framecount: ", rams%RR%framecount_rtp
                END IF
                IF (INDEX(to_lower(line), 'check_pade')>0) THEN !RTP time step
                    READ (line, *) dummy, rams%RR%check_pade
                    WRITE (*, *) "Apply Pade: ", rams%RR%check_pade
                END IF
                IF (INDEX(to_lower(line), 'pade_framecount')>0) THEN !RTP time step
                    READ (line, *) dummy, rams%RR%framecount_rtp_pade
                    WRITE (*, *) "Requested framecount after Pade: ", rams%RR%framecount_rtp_pade
                END IF
                IF (INDEX(to_lower(line), 'damping_constant')>0) THEN !RTP time step
                    READ (line, *) dummy, rams%RR%damping_constant
                    WRITE (*, *) "Damping constant (eV): ", rams%RR%damping_constant
                END IF
            END IF
            !IF (in_coordinates) THEN
            !IF (in_coordinates) THEN
            !    IF (INDEX(line, 'xyz_filename')>0) THEN
            !        READ (line, *) dummy, input%system%coordinates%xyz_filename
            !        input%system%coordinates%present = .TRUE.
            !    END IF
            !END IF
            !IF (in_fragments) THEN
            !    IF (INDEX(line, 'pdb_filename')>0) THEN
            !        READ (line, *) dummy, input%system%fragments%pdb_filename
            !        input%system%fragments%present = .TRUE.
            !    ELSEIF (INDEX(line, 'atom_list')>0) THEN
            !        ! parse atom lists, e.g.
            !        ! atom_list 4 1 2 3 4
            !        ! store in an allocatable array
            !    END IF
            !END IF
            IF (in_md) THEN
                IF (INDEX(line, 'time_step')>0) THEN
                    READ (line, *) dummy, md%dt
                    WRITE (*, *) "MD time step (fs): ", md%dt
                END IF
                IF (INDEX(line, 'correlation_depth')>0) THEN
                    READ (line, *) dummy, md%t_cor
                    WRITE (*, *) "Correlation depth: ", md%t_cor
                END IF
            END IF
        END DO
100     CONTINUE
        CLOSE (999)

    END SUBROUTINE parse_input

    SUBROUTINE check_input(gs, sys, md, stats, dips, rams)

        TYPE(global_settings)       :: gs
        TYPE(systems)               :: sys
        TYPE(molecular_dynamics)    :: md
        TYPE(static)                :: stats
        TYPE(dipoles)               :: dips
        TYPE(raman)                 :: rams

        IF (TRIM(gs%spectral_type%read_function)=='') THEN
            PRINT *, 'Error: Spectra not defined in the input'
            STOP
            !check for power spectrum
        ELSEIF (gs%spectral_type%read_function=='P') THEN
            !check for input_type
            IF (TRIM(sys%type_traj)=='') THEN
                PRINT *, 'Error: type_traj not defined in the input - provide "type_traj pos" for positions, &
        &                  "type_traj vel" for velocities'
                STOP
            END IF
            !check for filename
            IF (TRIM(sys%filename)=='') THEN
                PRINT *, 'Error: Filename not defined in the input'
                STOP
            END IF
            !check for mass_weighting
            IF (TRIM(sys%input_mass)=='') THEN
                PRINT *, 'Error: mass_weighting not defined in the input'
                STOP
            END IF
            !check time step
            IF (md%dt<0) THEN
                PRINT *, 'Error: time_step not defined in the input'
                STOP
            END IF
            !check t_cor
            IF (md%t_cor<0) THEN
                PRINT *, 'Error: correlation depth not defined in the input, we will continue with an estimate' !can be worded differently
                STOP
            END IF
            !check for the temperature
            IF (gs%temp<0) THEN
                PRINT *, 'Warning: Temperature is not defined, setting it to 300 K'
                gs%temp = 300
            END IF
            !check for incident laser wavelength
            !check for normal mode analysis
        ELSEIF (gs%spectral_type%read_function=='NMA') THEN
            !check for input_dipole not needed for P but set to a default value
            IF (TRIM(dips%type_dipole)=='') THEN
                PRINT *, 'Warning: type_dipole not defined in the input setting it to 1'
                dips%type_dipole = '1'
            END IF
            !check for filename
            IF (TRIM(sys%filename)=='') THEN
                PRINT *, 'Error: Filename not defined in the input'
                STOP
            END IF
            !check for force_file
            IF (TRIM(stats%force_file)=='') THEN
                PRINT *, 'Error: Force filename not defined in the input'
                STOP
            END IF
            !check for displacement
            IF (stats%dx<0) THEN
                PRINT *, 'Error: Displacement not defined in the input'
                STOP
            END IF
            !check for static IR
        ELSEIF (gs%spectral_type%read_function=='IR') THEN
            !check for filename
            IF (TRIM(sys%filename)=='') THEN
                PRINT *, 'Error: Filename not defined in the input'
                STOP
            END IF
            !check for force_file
            IF (stats%diag_hessian=='y') THEN
                IF (TRIM(stats%force_file)=='') THEN
                    PRINT *, 'Error: File name of the forces not defined in the input'
                    STOP
                END IF
            ELSEIF (stats%diag_hessian=='n') THEN
                IF (TRIM(stats%normal_freq_file)=='') THEN
                    PRINT *, 'Error: File name of the normal mode frequencies not defined in the input'
                    STOP
                END IF
                IF (TRIM(stats%normal_displ_file)=='') THEN
                    PRINT *, 'Error: File name of the normal mode displacements not defined in the input'
                    STOP
                END IF
            END IF
            !check for displacement
            IF (stats%dx<0) THEN
                PRINT *, 'Error: Displacement not defined in the input'
                STOP
            END IF
            !check for dipole file
            IF (TRIM(dips%dip_file)=='') THEN
                PRINT *, 'Error: Dipole filename not defined in the input'
                STOP
            END IF
            !check for type_dipole
            IF (TRIM(dips%type_dipole)=='') THEN
                PRINT *, 'Warning: type_dipole not defined in the input setting it to 1'
                dips%type_dipole = 'berry'
            END IF
            !check for the temperature
            IF (gs%temp<0) THEN
                PRINT *, 'Warning: Temperature is not defined, setting it to 300 K'
                gs%temp = 300
            END IF
            IF (gs%fwhm<0) THEN
                PRINT *, 'Warning: Full width at half-maximum is not defined, setting it to 10 cm^{-1}'
                gs%temp = 10
            END IF
            !check for static raman
        ELSEIF (gs%spectral_type%read_function=='R') THEN
            !check for filename
            IF (TRIM(sys%filename)=='') THEN
                PRINT *, 'Error: Filename not defined in the input'
                STOP
            END IF
            !check for force_file
            IF (stats%diag_hessian=='y') THEN
                IF (TRIM(stats%force_file)=='') THEN
                    PRINT *, 'Error: File name of the forces not defined in the input'
                    STOP
                END IF
            ELSEIF (stats%diag_hessian=='n') THEN
                IF (TRIM(stats%normal_freq_file)=='') THEN
                    PRINT *, 'Error: File name of the normal mode frequencies not defined in the input'
                    STOP
                END IF
                IF (TRIM(stats%normal_displ_file)=='') THEN
                    PRINT *, 'Error: File name of the normal mode displacements not defined in the input'
                    STOP
                END IF
            END IF
            !check for displacement
            IF (stats%dx<0) THEN
                PRINT *, 'Error: Displacement not defined in the input'
                STOP
            END IF
            !check for dipole file
            IF (TRIM(rams%static_pol_file)=='') THEN
                PRINT *, 'Error: Polarizability filename not defined in the input'
                STOP
            END IF
            !check for type_dipole
            IF (TRIM(dips%type_dipole)=='') THEN
                PRINT *, 'Warning: type_dipole not defined in the input setting it to 1'
                dips%type_dipole = 'dfpt'
            END IF
            IF (dips%type_dipole.NE.'dfpt' .AND. dips%e_field<0) THEN
                PRINT *, 'Error: Electric field strength not defined!'
                STOP
            END IF
            IF (rams%laser_in<0) THEN
                PRINT *, 'Warning: Incident laser frequency not defined, setting it to 1 0.5 cm⁻1'
                rams%laser_in = 0.5
            END IF
            !check for the temperature
            IF (gs%temp<0) THEN
                PRINT *, 'Warning: Temperature is not defined, setting it to 300 K'
                gs%temp = 300
            END IF
            IF (gs%fwhm<0) THEN
                PRINT *, 'Warning: Full width at half-maximum is not defined, setting it to 10 cm^{-1}'
                gs%temp = 10
            END IF

        ELSEIF (gs%spectral_type%read_function=='ABS') THEN
            !check for filename
            IF (TRIM(sys%filename)=='') THEN
                PRINT *, 'Error: Filename not defined in the input'
                STOP
            END IF
            !check for dipole file
            IF (TRIM(dips%dip_x_file)=='') THEN
                PRINT *, 'Error: X-field dipole file name not defined in the input'
                STOP
            END IF
            IF (TRIM(dips%dip_y_file)=='') THEN
                PRINT *, 'Error: Y-field dipole file name not defined in the input'
                STOP
            END IF
            IF (TRIM(dips%dip_z_file)=='') THEN
                PRINT *, 'Error: Z-field dipole file name not defined in the input'
                STOP
            END IF
            !check for type_dipole
            IF (TRIM(dips%type_dipole)=='') THEN
                PRINT *, 'Warning: type_dipole not defined in the input setting it to 1'
                dips%type_dipole = 'berry'
            END IF
            IF (rams%RR%dt_rtp<0) THEN
                PRINT *, 'Error: RTP time step not defined!'
                STOP
            END IF
            IF (dips%type_dipole.NE.'dfpt' .AND. dips%e_field<0) THEN
                PRINT *, 'Error: Electric field strength not defined!'
                STOP
            END IF
            IF (rams%RR%framecount_rtp<0) THEN
                PRINT *, 'Error: RTP framecount not defined!'
                STOP
            END IF
            IF (rams%RR%damping_constant<0) THEN
                PRINT *, 'Warning: Damping constant not defined, setting it to 0.1 eV!'
                rams%RR%damping_constant = 0.1_dp
            END IF
            IF (TRIM(rams%RR%check_pade)=='') THEN
                PRINT *, 'Warning: The calculation will continue without Pade approximants!'
                rams%RR%check_pade = 'y'
            END IF
            IF (TRIM(rams%RR%check_pade)=='y' .AND. rams%RR%framecount_rtp_pade<0) THEN !this can also be adjusted
                PRINT *, 'Warning: Pade framecount is set to 80000!'
                rams%RR%framecount_rtp_pade = 80000
            END IF

        !!Check for RR
        ELSEIF (gs%spectral_type%read_function=='RR') THEN
            !check for filename
            IF (TRIM(sys%filename)=='') THEN
                PRINT *, 'Error: Filename not defined in the input'
                STOP
            END IF
            !check for force_file
            IF (stats%diag_hessian=='y') THEN
                IF (TRIM(stats%force_file)=='') THEN
                    PRINT *, 'Error: File name of the forces not defined in the input'
                    STOP
                END IF
            ELSEIF (stats%diag_hessian=='n') THEN
                IF (TRIM(stats%normal_freq_file)=='') THEN
                    PRINT *, 'Error: File name of the normal mode frequencies not defined in the input'
                    STOP
                END IF
                IF (TRIM(stats%normal_displ_file)=='') THEN
                    PRINT *, 'Error: File name of the normal mode displacements not defined in the input'
                    STOP
                END IF
            END IF
            !check for displacement
            IF (stats%dx<0) THEN
                PRINT *, 'Error: Displacement not defined in the input'
                STOP
            END IF
            !check for type_dipole
            IF (TRIM(dips%type_dipole)=='') THEN
                PRINT *, 'Warning: type_dipole not defined in the input setting it to 1'
                dips%type_dipole = 'berry'
            END IF
            !check for dipole files
            IF (TRIM(dips%dip_x_file)=='') THEN
                PRINT *, 'Error: X-field dipole file name not defined in the input'
                STOP
            END IF
            IF (TRIM(dips%dip_y_file)=='') THEN
                PRINT *, 'Error: Y-field dipole file name not defined in the input'
                STOP
            END IF
            IF (TRIM(dips%dip_z_file)=='') THEN
                PRINT *, 'Error: Z-field dipole file name not defined in the input'
                STOP
            END IF
            IF (rams%RR%dt_rtp<0) THEN
                PRINT *, 'Error: RTP time step not defined!'
                STOP
            END IF
            IF (dips%type_dipole.NE.'dfpt' .AND. dips%e_field<0) THEN
                PRINT *, 'Error: Electric field strength not defined!'
                STOP
            END IF
            IF (rams%RR%framecount_rtp<0) THEN
                PRINT *, 'Error: RTP framecount not defined!'
                STOP
            END IF
            IF (TRIM(rams%RR%check_pade)=='') THEN
                PRINT *, 'Warning: The calculation will continue without Pade approximants!'
                rams%RR%check_pade = 'y'
            END IF
            IF (TRIM(rams%RR%check_pade)=='y' .AND. rams%RR%framecount_rtp_pade<0) THEN !this can also be adjusted
                PRINT *, 'Warning: Pade framecount is set to 80000!'
                rams%RR%framecount_rtp_pade = 80000
            END IF
            IF (rams%RR%damping_constant<0) THEN
                PRINT *, 'Warning: Damping constant not defined, setting it to 0.1 eV!'
                rams%RR%damping_constant = 0.1_dp
            END IF
            IF (rams%laser_in<0) THEN
                PRINT *, 'Warning: Incident laser frequency not defined, setting it to 1 0.5 cm⁻1'
                rams%laser_in = 0.5
            END IF
            !check for the temperature
            IF (gs%temp<0) THEN
                PRINT *, 'Warning: Temperature is not defined, setting it to 300 K'
                gs%temp = 300
            END IF
            IF (gs%fwhm<0) THEN
                PRINT *, 'Warning: Full width at half-maximum is not defined, setting it to 10 cm^{-1}'
                gs%temp = 10
            END IF

            !check for MD-IR
        ELSEIF (gs%spectral_type%read_function=='MD-IR') THEN
            !check for dipole file
            IF (TRIM(dips%dip_file)=='') THEN
                PRINT *, 'Error: Dipole filename not defined in the input'
                STOP
            END IF
            !check for type_dipole
            IF (TRIM(dips%type_dipole)=='') THEN
                PRINT *, 'Warning: type_dipole not defined in the input setting it to 1'
                dips%type_dipole = 'berry'
            END IF
            !check for dipole file
            IF (TRIM(dips%dip_file)=='') THEN
                PRINT *, 'Error: Dipole filename not defined in the input'
                STOP
            END IF
            !check time step
            IF (md%dt<0) THEN
                PRINT *, 'Error: time_step not defined in the input'
                STOP
            END IF
            !check t_cor
            IF (md%t_cor<0) THEN
                PRINT *, 'Error: correlation depth not defined in the input, we will continue with an estimate' !can be worded differently
                STOP
            END IF
            !check for the temperature
            IF (gs%temp<0) THEN
                PRINT *, 'Warning: Temperature is not defined, setting it to 300 K'
                gs%temp = 300
            END IF
            !check for incident laser wavelength
            !check for MD-Raman
        ELSEIF (gs%spectral_type%read_function=='MD-R') THEN
            !check for dipole file
            IF (TRIM(dips%dip_file)=='') THEN
                PRINT *, 'Error: Dipole filename not defined in the input'
                STOP
            END IF
            !check for type_dipole
            IF (TRIM(dips%type_dipole)=='') THEN
                PRINT *, 'Warning: type_dipole not defined in the input setting it to berry'
                dips%type_dipole = 'berry'
            END IF
            !check for electric field strength
            IF (dips%type_dipole.NE.'dfpt' .AND. dips%e_field<0) THEN
                PRINT *, 'Error: Electric field strength not defined!'
                STOP
            END IF
            !check for dipole file
            IF (TRIM(dips%dip_file)=='') THEN
                PRINT *, 'Error: Dipole filename not defined in the input'
                STOP
            END IF
            IF (TRIM(dips%dip_x_file)=='') THEN
                PRINT *, 'Error: X-field dipole file name not defined in the input'
                STOP
            END IF
            IF (TRIM(dips%dip_y_file)=='') THEN
                PRINT *, 'Error: Y-field dipole file name not defined in the input'
                STOP
            END IF
            IF (TRIM(dips%dip_z_file)=='') THEN
                PRINT *, 'Error: Z-field dipole file name not defined in the input'
                STOP
            END IF
            IF (rams%RR%dt_rtp<0) THEN
                PRINT *, 'Error: RTP time step not defined!'
                STOP
            END IF
            !check time step
            IF (md%dt<0) THEN
                PRINT *, 'Error: time_step not defined in the input'
                STOP
            END IF
            !check t_cor
            IF (md%t_cor<0) THEN
                PRINT *, 'Error: correlation depth not defined in the input, we will continue with an estimate' !can be worded differently
                STOP
            END IF
            !check for the temperature
            IF (gs%temp<0) THEN
                PRINT *, 'Warning: Temperature is not defined, setting it to 300 K'
                gs%temp = 300
            END IF
            !check for incident laser wavelength
            IF (rams%laser_in<0) THEN
                PRINT *, 'Warning: Incident laser frequency not defined, setting it to 1 0.5 cm⁻1'
                rams%laser_in = 0.5
            END IF

        END IF

    END SUBROUTINE check_input
END MODULE read_input
