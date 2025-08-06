MODULE read_input

    USE kinds, ONLY: dp, str_len
    USE vib_types, ONLY: global_settings, systems, static, dipoles, raman, molecular_dynamics

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: parse_command_line,  parse_input, check_input

CONTAINS

    !****************************************************************
    ! doxygen doc, to be added
    !****************************************************************
    SUBROUTINE parse_command_line(input_file_name)

        CHARACTER(LEN=str_len), INTENT(OUT) :: input_file_name

        CHARACTER(LEN=str_len)                :: arg
        INTEGER(KIND=4)                                     :: narg
        INTEGER                                             :: i, stat

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

        !** intermal variables
        CHARACTER(LEN=str_len) :: line
        CHARACTER(LEN=str_len) :: dummy
        LOGICAL :: in_global= .FALSE.
        LOGICAL :: in_system = .FALSE.
        LOGICAL :: in_cell = .FALSE.
        LOGICAL :: in_coordinates = .FALSE.
        LOGICAL :: in_fragments = .FALSE.
        LOGICAL :: in_md = .FALSE.
        LOGICAL :: in_static = .FALSE.
        LOGICAL :: in_hessian = .FALSE.
        
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
            ENDIF
  
            IF (INDEX(line, '&end hessian')>0) THEN
                in_hessian = .FALSE.
                CYCLE
            END IF            

            ! Parse within active sections
            IF (in_global) THEN
                IF (INDEX(to_lower(line), 'temperature')>0) THEN
                    READ (line, *) dummy, gs%temp
                    write(*,*) "temperature: ", gs%temp
                    !input%system%cell%present = .TRUE. ! add is present later
                ELSEIF (INDEX(to_lower(line), 'laser_in')>0) THEN
                    READ (line, *) dummy, gs%laser_in
                    write(*,*) "laser_in: ", gs%laser_in
                ELSEIF (INDEX(to_lower(line), 'spectra')>0) THEN
                        READ (line, *) dummy, gs%spectral_type%read_function
                        write(*,*) "spectra: ", gs%spectral_type%read_function
                ELSEIF (INDEX(to_lower(line), 'type_static')>0) THEN
                    READ (line, *) dummy, gs%spectral_type%type_static !'Do you want the normal modes to be calculated (type "1") or read from an external file (type "2")?'
                    write(*,*) "type_static: ", gs%spectral_type%type_static
                ELSEIF (INDEX(to_lower(line), 'type_dipole')>0) THEN  !Which one do you want to use: Wannier centers (1) or Berry phase dipole moments (2)?'
                    READ (line, *) dummy, gs%spectral_type%type_dipole
                    write(*,*) "type_dipole: ", gs%spectral_type%type_dipole
                END IF
            END IF
            
            IF (in_system) THEN
                IF (INDEX(to_lower(line), 'filename')>0) THEN
                    READ (line, *) dummy, sys%filename
                    write(*,*) "filename: ", sys%filename
                ELSEIF (INDEX(to_lower(line), 'type_traj')>0) THEN
                    READ (line, *) dummy, sys%type_traj !'Enter the type of the input file (type 1 for positions, 2 for velocities)'
                    write(*,*) "type_traj: ", sys%type_traj
                ELSEIF (INDEX(to_lower(line), 'mass_weighting')>0) THEN !'Do you want to apply mass weighting (y/n)?
                    READ (line, *) dummy, sys%input_mass
                    write(*,*) "mass_weighting: ",  sys%input_mass
                ELSEIF (in_cell) THEN
                    IF (INDEX(to_lower(line), 'cell_type')>0) THEN !'Is it the k-point trajectory (1), supercell trajectory (2) or solvent trajectory (3)? '
                      READ (line, *) dummy, sys%cell%cell_type
                !    ELSEIF (INDEX(to_lower(line), 'abc')>0) THEN
                !        !READ (line, *) dummy, sys%cell%abc(1:3)
                !        !input%system%cell%present = .TRUE.
                !    ELSEIF (INDEX(to_lower(line), 'alpha_beta_gamma')>0) THEN
                !        !READ (line, *) dummy, input%system%cell%alpha_beta_gamma(1:3)
                    END IF
                ELSEIF (INDEX(to_lower(line), 'periodic ')>0) THEN !'Does the system contain more than one molecule? (y/n)'
                    READ (line, *) dummy, sys%periodic
                    write(*,*) "periodic: ",  sys%periodic
                ELSEIF (INDEX(to_lower(line), 'frag_type ')>0) THEN !'Does the system contain more than one molecule? (y/n)'
                    READ (line, *) dummy, sys%frag_type
                    write(*,*) "frag_type: ",  sys%frag_type
                END IF
            END IF

            IF (in_static) THEN
                IF (in_hessian) THEN                
                    IF (INDEX(to_lower(line), 'force_file')>0) THEN
                        READ (line, *) dummy, stats%force_file
                        write(*,*) "force_file: ", stats%force_file
                    ENDIF
                ENDIF
                IF (INDEX(to_lower(line), 'displacement')>0) THEN !'Do you want to apply mass weighting (y/n)?
                    READ (line, *) dummy, stats%dx
                    write(*,*) "displacement: ",  stats%dx
                END iF

            END IF

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
                    write(*,*) "time_step: ",  md%dt
                END IF
                IF (INDEX(line, 'correlation_depth')>0) THEN
                    READ (line, *) dummy, md%t_cor
                    write(*,*) "correlation depth: ",  md%t_cor
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

        IF (trim(gs%spectral_type%read_function) == '') THEN
            print *, 'Error: Spectra not defined in the input'
            stop
        
        ELSEIF (gs%spectral_type%read_function=='P') THEN
            !check for input_type
            IF (trim(sys%type_traj) == '') THEN
                print *, 'Error: type_traj not defined in the input - provide "type_traj pos" for positions, "type_traj vel" for velocities'
                stop
            END IF
            !check for input_dipole, not needed for P but set to a default value
            IF (trim(gs%spectral_type%type_dipole) == '') THEN !<----- NEEDED IN READ_COORD FUNCTION ...
                print *, 'Warning: type_dipole not defined in the input setting it to 1'
                gs%spectral_type%type_dipole = '1'
            !check for filename
            END IF
            IF (trim(sys%filename) == '') THEN
                print *, 'Error: Filename not defined in the input'
                stop
            END IF
            !check for mass_weighting 
            IF (trim(sys%input_mass ) == '') THEN
                print *, 'Error: mass_weighting not defined in the input'
                stop
            END IF
            !check time step
            IF (md%dt < 0) THEN
                print *, 'Error: time_step not defined in the input'
                stop
            END IF
            !check t_cor
            IF (md%t_cor < 0) THEN
                print *, 'Error: correlation depth not defined in the input, we will continue with an estimate' !can be worded differently
                stop
            ENDIF
 
        ELSEIF (gs%spectral_type%read_function=='NMA') THEN
            !check for input_dipole not needed for P but set to a default value
            IF (trim(gs%spectral_type%type_dipole) == '') THEN
                print *, 'Warning: type_dipole not defined in the input setting it to 1'
                gs%spectral_type%type_dipole = '1'
            END IF
            !check for filename
            IF (trim(sys%filename) == '') THEN
                print *, 'Error: Filename not defined in the input'
                stop
            END IF
            !check for force_file
            IF (trim(stats%force_file) == '') THEN
                print *, 'Error: Force filename not defined in the input'
                stop
            END IF
            !check for displacement
            IF (stats%dx  < 0) THEN
                print *, 'Error: Displacement not defined in the input'
                stop
            END IF

        END IF


   END SUBROUTINE check_input
END MODULE read_input
