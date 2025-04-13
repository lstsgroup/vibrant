MODULE read_traj

    USE kinds, ONLY: dp

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: read_coord, read_coord_frame, read_normal_modes, read_static, read_static_resraman

CONTAINS
    SUBROUTINE read_coord(natom, framecount, element, coord, filename, periodic, mol_num, system, read_function, &
                          framecount_rtp, type_dipole)

        CHARACTER(LEN=40), INTENT(IN)                               :: filename, periodic, system, read_function, type_dipole
        INTEGER, INTENT(IN)                                         :: framecount_rtp
        INTEGER, INTENT(OUT)                                        :: natom, framecount, mol_num
        CHARACTER(LEN=2), DIMENSION(:), ALLOCATABLE, INTENT(OUT)      :: element
        REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE, INTENT(OUT)        :: coord

        INTEGER                                                    :: i, j, stat

        framecount = 0

        IF (read_function.NE.'MD-RR') THEN
            OPEN (UNIT=50, FILE=filename, STATUS='old', IOSTAT=stat)
            READ (50, *) natom
            CLOSE (50)
        ELSEIF (read_function=='MD-RR') THEN
            natom = framecount_rtp
        END IF

        ALLOCATE (element(natom), coord(natom, 3))

        OPEN (UNIT=51, FILE=filename, STATUS='old', IOSTAT=stat)
        DO
            READ (51, *, END=998)
            READ (51, *)
            framecount = framecount + 1
            DO i = 1, natom
                READ (51, *) element(i), coord(i, 1), coord(i, 2), coord(i, 3)
            END DO
        END DO
998     CONTINUE
        CLOSE (51)

        IF (type_dipole=='2' .OR. type_dipole=='3') THEN !!gas phase
            mol_num = 1
        ELSEIF ((periodic=='n' .AND. system=='1') .OR. type_dipole=='1') THEN !!fragment approach
            mol_num = 44 !20 !! fix later to 20
        END IF
        PRINT *, mol_num, 'mol num'
    END SUBROUTINE read_coord

!********************************************************************************************
!********************************************************************************************
    SUBROUTINE read_coord_frame(natom, framecount, element, filename, coord_v)
        CHARACTER(LEN=40), INTENT(IN)                               :: filename
        INTEGER, INTENT(INOUT)                                      :: natom, framecount
        CHARACTER(LEN=2), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)    :: element
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE, INTENT(OUT)      :: coord_v

        INTEGER                                                    :: i, j, stat

        ALLOCATE (coord_v(framecount, natom, 3))
        OPEN (UNIT=52, FILE=filename, STATUS='old', IOSTAT=stat)
        DO
            DO j = 1, framecount
                READ (52, *, END=999)
                READ (52, *)
                DO i = 1, natom
                    READ (52, *) element(i), coord_v(j, i, 1), coord_v(j, i, 2), coord_v(j, i, 3)
                END DO
            END DO
        END DO
999     CONTINUE
        CLOSE (52)

    END SUBROUTINE read_coord_frame

!********************************************************************************************
!********************************************************************************************

    SUBROUTINE read_normal_modes(natom, element, normal_freq_file, normal_displ_file, freq, disp, nmodes, &
                                 read_function, type_static, force_file, force)

        CHARACTER(LEN=40), INTENT(IN)                               :: read_function, type_static
        CHARACTER(LEN=40), INTENT(IN)                               :: normal_freq_file, normal_displ_file, force_file
        INTEGER, INTENT(INOUT)                                      :: natom
        INTEGER, INTENT(OUT)                                        :: nmodes
        CHARACTER(LEN=2), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)    :: element
        REAL(kind=dp), DIMENSION(:, :, :, :, :), ALLOCATABLE, INTENT(OUT)  :: force
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(OUT)          :: freq
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE, INTENT(OUT)      :: disp

        CHARACTER(LEN=40)                                          :: chara
        INTEGER                                                    :: i, j, k, m, n
        INTEGER                                                    :: stat

        IF (read_function=='NMA' .OR. type_static=='1') THEN
            ALLOCATE (force(2, natom, 3, natom, 3))
            OPEN (UNIT=49, FILE=force_file, STATUS='old', IOSTAT=stat) !Reading forces
            DO i = 1, 2
                DO j = 1, natom
                    DO m = 1, 3
                        DO n = 1, natom
                            READ (49, *) force(i, j, m, n, 1), force(i, j, m, n, 2), force(i, j, m, n, 3)
                        END DO
                    END DO
                END DO
            END DO
            CLOSE (49)

        ELSEIF (type_static=='2') THEN
            nmodes = 0
            OPEN (UNIT=50, FILE=normal_freq_file, STATUS='old', IOSTAT=stat) !Reading normal freqs/coords
            DO
                READ (50, *, END=998) chara
                nmodes = nmodes + 1
            END DO
998         CONTINUE
            CLOSE (50)

            ALLOCATE (freq(nmodes), disp(nmodes, natom, 3))
            
            OPEN (UNIT=51, FILE=normal_freq_file, STATUS='old', IOSTAT=stat) !Reading normal freqs/coords
            DO i = 1, nmodes
                READ (51, *, END=997) freq(i)
            END DO
997         CONTINUE
            CLOSE (51)

            OPEN (UNIT=51, FILE=normal_displ_file, STATUS='old', IOSTAT=stat) !Reading normal freqs/coords
            DO i = 1, nmodes
                DO j = 1, natom
                    READ (51, *, END=996) disp(i, j, 1), disp(i, j, 2), disp(i, j, 3)
                END DO
            END DO
996         CONTINUE
            CLOSE (51)
        END IF


    END SUBROUTINE read_normal_modes

!********************************************************************************************
!********************************************************************************************

    SUBROUTINE read_static(natom, element, static_pol_file, pol, static_dip_file, type_dipole, static_dip, &
                           type_static)

        CHARACTER(LEN=40), INTENT(IN)                               :: static_pol_file, static_dip_file
        CHARACTER(LEN=40), INTENT(IN)                               :: type_static, type_dipole
        INTEGER, INTENT(INOUT)                                      :: natom
        CHARACTER(LEN=2), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)    :: element
        REAL(kind=dp), DIMENSION(:, :, :, :, :), ALLOCATABLE, INTENT(OUT)  :: pol
        REAL(kind=dp), DIMENSION(:, :, :, :), ALLOCATABLE, INTENT(OUT)    :: static_dip

        CHARACTER(LEN=40)                                          :: chara
        INTEGER                                                    :: i, j, k, m, n
        INTEGER                                                    :: stat

        ALLOCATE (pol(natom, 3, 2, 3, 3), static_dip(natom, 3, 2, 3))

        IF (type_dipole=='3') THEN
            OPEN (UNIT=52, FILE=static_pol_file, STATUS='old', IOSTAT=stat) !Reading polarizabilties
            DO
                    DO k = 1, 2
                        DO i = 1, natom
                            DO j = 1, 3
                                READ (52, *, END=995)
                                READ (52, *)
                                READ (52, *)
                                READ (52, *)
                                READ (52, *)
                                READ (52, *)
                                READ (52, *) chara, chara, chara, pol(i, j, k, 1, 1), pol(i, j, k, 2, 2), pol(i, j, k, 3, 3)
                                READ (52, *) chara, chara, chara, pol(i, j, k, 1, 2), pol(i, j, k, 1, 3), pol(i, j, k, 2, 3)
                                READ (52, *) chara, chara, chara, pol(i, j, k, 2, 1), pol(i, j, k, 3, 1), pol(i, j, k, 3, 2)
                            END DO
                        END DO
                    END DO
                END DO
995             CONTINUE
                CLOSE (52)

        ELSEIF (type_dipole=='2') THEN
                OPEN (UNIT=53, FILE=static_dip_file, STATUS='old', IOSTAT=stat) !Reading dipoles
                DO
                    DO k = 1, 2
                        DO i = 1, natom
                            DO j = 1, 3
                                READ (53, *, END=994)
                                READ (53, *)
                                READ (53, *) chara, static_dip(i, j, k, 1), static_dip(i, j, k, 2), static_dip(i, j, k, 3)
                            END DO
                        END DO
                    END DO
                END DO
994             CONTINUE
                CLOSE (53)
            
        END IF

    END SUBROUTINE read_static

!********************************************************************************************
!********************************************************************************************

    SUBROUTINE read_static_resraman(natom, element, static_dip_file, framecount_rtp, static_dip_rtp)

        CHARACTER(LEN=40), INTENT(IN)                               :: static_dip_file
        CHARACTER(LEN=2), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)    :: element
        INTEGER, INTENT(INOUT)                                      :: natom, framecount_rtp
        REAL(kind=dp), DIMENSION(:, :, :, :, :), ALLOCATABLE, INTENT(OUT)  :: static_dip_rtp

        CHARACTER(LEN=40)                                          :: chara
        INTEGER                                                    :: i, j, k, m, stat

        ALLOCATE (static_dip_rtp(natom, 3, 2, 3, framecount_rtp + 1))

        OPEN (UNIT=53, FILE=static_dip_file, STATUS='old', IOSTAT=stat) !Reading polarizabilties
        DO
            DO k = 1, 2
                DO i = 1, natom
                    DO j = 1, 3
                        READ (53, *, END=994)
                        READ (53, *)
                        DO m = 1, framecount_rtp + 1
                            READ (53, *) chara, static_dip_rtp(i, j, k, 1, m), static_dip_rtp(i, j, k, 2, m), &
                                static_dip_rtp(i, j, k, 3, m)
                        END DO
                    END DO
                END DO
            END DO
        END DO
994     CONTINUE
        CLOSE (53)

    END SUBROUTINE read_static_resraman
END MODULE read_traj
