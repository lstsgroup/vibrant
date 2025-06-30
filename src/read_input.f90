MODULE read_input

 USE types,    ONLY: section_type
 USE kinds,    ONLY: dp, default_string_length
 
 IMPLICIT NONE

 PRIVATE
 
 PUBLIC :: parse_command_line, parse_input

CONTAINS

!****************************************************************
! doxygen doc, to be added
!****************************************************************
 SUBROUTINE parse_command_line(input, input_file_name)

 TYPE(section_type) :: input
 CHARACTER(LEN=default_string_length), INTENT(OUT) :: input_file_name

 CHARACTER(LEN=default_string_length)                :: arg
 INTEGER(KIND=4)                                     :: narg
 INTEGER                                             :: i, stat

 narg = iargc()

 IF (narg /= 1) THEN
   WRITE(*,'(A37)') "Usage: vibrant_input.x your_input.inp" 
   STOP
 ENDIF

 CALL getarg(1, input_file_name)

 END SUBROUTINE parse_command_line
!****************************************************************
! doxygen doc, to be added
!****************************************************************
 SUBROUTINE parse_input(input,input_file_name)

 TYPE(section_type) :: input
 CHARACTER(LEN=default_string_length), INTENT(IN) :: input_file_name

 !** intermal variables
 CHARACTER(LEN=default_string_length) :: line
 CHARACTER(LEN=default_string_length) :: dummy
 LOGICAL :: in_system = .FALSE.
 LOGICAL :: in_cell = .FALSE.
 LOGICAL :: in_coordinates = .FALSE.
 LOGICAL :: in_fragments = .FALSE.
 LOGICAL :: in_md = .FALSE.
 

 OPEN(unit=999, file=TRIM(input_file_name), status="old")
  
 DO
    READ(999,'(A)',end=100) line
    line = ADJUSTL(line)

    ! Identify section starts
    IF (INDEX(line,'&system') > 0) THEN
        in_system = .TRUE.
        CYCLE
    ENDIF

    IF (INDEX(line,'&end system') > 0) THEN
        in_system = .FALSE.
        CYCLE
    ENDIF

    IF (INDEX(line,'&cell') > 0) THEN
        in_cell = .TRUE.
        CYCLE
    ENDIF

    IF (INDEX(line,'&end cell') > 0) THEN
        in_cell = .FALSE.
        CYCLE
    ENDIF

    IF (INDEX(line,'&coordinates') > 0) THEN
        in_coordinates = .TRUE.
        CYCLE
    ENDIF

    IF (INDEX(line,'&end coordinates') > 0) THEN
        in_coordinates = .FALSE.
        CYCLE
    ENDIF

    IF (INDEX(line,'&fragments') > 0) THEN
        in_fragments = .TRUE.
        CYCLE
    ENDIF

    IF (INDEX(line,'&end fragments') > 0) THEN
        in_fragments = .FALSE.
        CYCLE
    ENDIF

    IF (INDEX(line,'&md') > 0) THEN
        in_md = .TRUE.
        CYCLE
    ENDIF

    IF (INDEX(line,'&end md') > 0) THEN
        in_md = .FALSE.
        CYCLE
    ENDIF

    ! Parse within active sections
    IF (in_cell) THEN
        IF (INDEX(line,'ABC') > 0) THEN
            READ(line,*) dummy, input%system%cell%abc(1:3)
            input%system%cell%present = .TRUE.
        ELSEIF (INDEX(line,'alpha_beta_gamma') > 0) THEN
            READ(line,*) dummy, input%system%cell%alpha_beta_gamma(1:3)
        ENDIF
    ENDIF

    IF (in_coordinates) THEN
        IF (INDEX(line,'xyz_filename') > 0) THEN
            READ(line,*) dummy, input%system%coordinates%xyz_filename
            input%system%coordinates%present = .TRUE.
        ENDIF
    ENDIF

    IF (in_fragments) THEN
        IF (INDEX(line,'pdb_filename') > 0) THEN
            READ(line,*) dummy, input%system%fragments%pdb_filename
            input%system%fragments%present = .TRUE.
        ELSEIF (INDEX(line,'atom_list') > 0) THEN
            ! parse atom lists, e.g.
            ! atom_list 4 1 2 3 4
            ! store in an allocatable array
        ENDIF
    ENDIF

    IF (in_md) THEN
        IF (INDEX(line,'trajectory_file') > 0) THEN
            READ(line,*) dummy, input%md%trajectory_file
        ELSEIF (INDEX(line,'velocity_file') > 0) THEN
            READ(line,*) dummy, input%md%velocity_file
        ELSEIF (INDEX(line,'snapshot_time_step') > 0) THEN
            READ(line,*) dummy, input%md%snapshot_time_step
        ELSEIF (INDEX(line,'correlation_depth') > 0) THEN
            READ(line,*) dummy, input%md%correlation_depth
        ENDIF
    ENDIF 
 ENDDO
100 CONTINUE
 CLOSE(999)
 
 END SUBROUTINE parse_input
END MODULE read_input
