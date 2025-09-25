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

MODULE read_traj

    USE kinds, ONLY: dp, str_len
    USE iso_fortran_env, ONLY: output_unit, error_unit
    USE vib_types, ONLY: global_settings, systems, static, dipoles, raman, molecular_dynamics, static_property

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: read_coord, read_coord_frame, read_normal_modes, read_static, read_static_resraman, check_file_open

CONTAINS

    SUBROUTINE check_file_open(stat, msg, filename)
        INTEGER,      INTENT(IN) :: stat
        CHARACTER(*), INTENT(IN) :: msg
        CHARACTER(*), INTENT(IN) :: filename

        IF (stat /= 0) THEN
            WRITE(error_unit,'(4X,"[ERROR] could not open file ",A)') TRIM(filename)
            WRITE(error_unit,'(4X,"I/O error message: ",A)')         TRIM(msg)
        STOP 
        END IF
    END SUBROUTINE check_file_open

    SUBROUTINE read_coord(filename, gs, sys, dips, rams)

        TYPE(global_settings), INTENT(INOUT)   :: gs
        TYPE(systems), INTENT(INOUT)        :: sys
        TYPE(dipoles), OPTIONAL        :: dips
        TYPE(raman), OPTIONAL        :: rams
        CHARACTER(LEN=40), INTENT(IN)                               :: filename

        CHARACTER(len=str_len)                                     :: msg  ! store error message
        INTEGER                                                   :: i, j, stat, runit

        sys%framecount = 0

        IF (gs%spectral_type%read_function/='MD-RR') THEN
            OPEN (FILE=filename, STATUS='old', ACTION='read', IOSTAT=stat, IOMSG=msg, NEWUNIT=runit)
            READ (runit, *) sys%natom
            CLOSE (runit)
        ELSEIF (gs%spectral_type%read_function=='MD-RR') THEN
            sys%natom = rams%RR%framecount_rtp
        END IF

        ALLOCATE (sys%element(sys%natom), sys%coord(sys%natom, 3))

        OPEN (FILE=filename, STATUS='old', ACTION='read', IOSTAT=stat, IOMSG=msg, NEWUNIT=runit)
        !Check if file exists
        CALL check_file_open(stat, msg, filename)
        DO
            READ (runit, *, END=998)
            READ (runit, *)
            sys%framecount = sys%framecount + 1
            DO i = 1, sys%natom
                READ (runit, *) sys%element(i), sys%coord(i, 1), sys%coord(i, 2), sys%coord(i, 3)
            END DO
        END DO
998     CONTINUE
        CLOSE (runit)

        IF (gs%spectral_type%read_function/='P') THEN
            IF (dips%type_dipole=='berry' .OR. dips%type_dipole=='dfpt' .OR. dips%type_dipole=='wannier') THEN !!gas phase
                sys%mol_num = 1
                ! ELSEIF ((sys%periodic=='n' .AND. sys%system=='1') .OR. dips%type_dipole=='wannier') THEN !!fragment approach
                !     sys%mol_num = 44 !20 !! fix later to 20
            END IF
        END IF
       
        WRITE (*, '(4X,"Number of mols", T60, I0)') sys%mol_num
        WRITE (*, '(4X,"Number of frames", T60, I0)')   sys%framecount

    END SUBROUTINE read_coord

!********************************************************************************************
!********************************************************************************************
    SUBROUTINE read_coord_frame(natom, filename, coord_v, sys)

        TYPE(systems), INTENT(INOUT)        :: sys
        CHARACTER(LEN=40), INTENT(IN)                               :: filename
        INTEGER, INTENT(INOUT)                               :: natom
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE, INTENT(OUT)      :: coord_v

        CHARACTER(len=str_len)                                     :: msg  ! store error message
        INTEGER                                                    :: i, j, stat, runit

        ALLOCATE (coord_v(sys%framecount, natom, 3))
        OPEN (FILE=filename, STATUS='old', ACTION='read', IOSTAT=stat, IOMSG=msg, NEWUNIT=runit)
        !Check if file exists
        CALL check_file_open(stat, msg, filename)
        !Start reading if file found
        DO
            DO j = 1, sys%framecount
                READ (runit, *, END=999)
                READ (runit, *)
                DO i = 1, natom
                    READ (runit, *) sys%element(i), coord_v(j, i, 1), coord_v(j, i, 2), coord_v(j, i, 3)
                END DO
            END DO
        END DO
999     CONTINUE
        CLOSE (runit)

    END SUBROUTINE read_coord_frame
!********************************************************************************************
!********************************************************************************************

    SUBROUTINE read_normal_modes(gs, sys, stats)

        ! Variables of your derived types:
        TYPE(global_settings), INTENT(INOUT)   :: gs
        TYPE(systems), INTENT(INOUT)        :: sys
        TYPE(static), INTENT(INOUT)        :: stats

        CHARACTER(LEN=str_len)                                          :: chara, msg
        INTEGER                                                    :: i, j, k, m, n, d, xyz
        INTEGER                                                    :: stat, runit

        IF (gs%spectral_type%read_function=='NMA' .OR. stats%diag_hessian=='y') THEN

            CALL stats%init_force(sys%natom, 1)

            OPEN (FILE=stats%force_file, STATUS='old', ACTION='read', IOSTAT=stat, IOMSG=msg, NEWUNIT=runit)
            !Check if file exists
            CALL check_file_open(stat, msg, stats%force_file)
            DO i = 1, 2
                DO j = 1, sys%natom
                    DO m = 1, 3
                        DO n = 1, sys%natom
                            READ (runit, *) stats%force(n, 1)%atom(j)%displacement(i)%XYZ(m)%frame(1), &
                                stats%force(n, 2)%atom(j)%displacement(i)%XYZ(m)%frame(1), &
                                stats%force(n, 3)%atom(j)%displacement(i)%XYZ(m)%frame(1)
                        END DO
                    END DO
                END DO
            END DO
            CLOSE (runit)

        ELSEIF (stats%diag_hessian=='n') THEN
            stats%nmodes = 0
            OPEN (FILE=stats%normal_freq_file, STATUS='old', ACTION='read', IOSTAT=stat, IOMSG=msg, NEWUNIT=runit) !Reading normal freqs/coords
            !Check if file exists
            CALL check_file_open(stat, msg, stats%normal_freq_file)
            DO
                READ (runit, *, END=998) chara
                stats%nmodes = stats%nmodes + 1
            END DO
998         CONTINUE
            CLOSE (runit)

            ALLOCATE (stats%freq(stats%nmodes), stats%disp(stats%nmodes, sys%natom, 3))

            OPEN (FILE=stats%normal_freq_file, STATUS='old', ACTION='read', IOSTAT=stat, IOMSG=msg, NEWUNIT=runit) !Reading normal freqs/coords
            !Check if file exists
            CALL check_file_open(stat, msg, stats%normal_freq_file)
            DO i = 1, stats%nmodes
                READ (runit, *, END=997) stats%freq(i)
            END DO
997         CONTINUE
            CLOSE (runit)

            OPEN (FILE=stats%normal_displ_file, STATUS='old', ACTION='read', IOSTAT=stat, IOMSG=msg, NEWUNIT=runit) !Reading normal freqs/coords
            !Check if file exists
            CALL check_file_open(stat, msg, stats%normal_displ_file)
            DO i = 1, stats%nmodes
                DO j = 1, sys%natom
                    READ (runit, *, END=996) stats%disp(i, j, 1), stats%disp(i, j, 2), stats%disp(i, j, 3)
                END DO
            END DO
996         CONTINUE
            CLOSE (runit)
        END IF

    END SUBROUTINE read_normal_modes

!********************************************************************************************
!********************************************************************************************

    SUBROUTINE read_static(gs, sys, dips, rams)
        TYPE(global_settings), INTENT(INOUT)   :: gs
        TYPE(systems), INTENT(INOUT)        :: sys
        TYPE(dipoles), INTENT(INOUT), OPTIONAL        :: dips
        TYPE(raman), INTENT(INOUT), OPTIONAL        :: rams

        CHARACTER(LEN=str_len)                                          :: chara, msg

        INTEGER                                                    :: i, j, k, m, n, d, runit
        INTEGER                                                    :: stat, i_pol, j_pol, xyz

        !IF (PRESENT(static_dip)) ALLOCATE (static_dip(sys%natom, 3, 2, 3))
        IF (PRESENT(dips)) CALL dips%init_dip(sys%natom, 1)
        IF (PRESENT(rams)) CALL rams%init_pol(sys%natom, 1)

        IF (dips%type_dipole=='dfpt') THEN
            OPEN (FILE=rams%static_pol_file, STATUS='old', ACTION='read', IOSTAT=stat, IOMSG=msg, NEWUNIT=runit) !Reading polarizabilties
            !Check if file exists
            CALL check_file_open(stat, msg, rams%static_pol_file)
            DO
                DO k = 1, 2
                    DO i = 1, sys%natom
                        DO j = 1, 3
                            READ (runit, *, END=995)
                            READ (runit, *)
                            READ (runit, *)
                            READ (runit, *)
                            READ (runit, *)
                            READ (runit, *)
                            READ (runit, *) chara, chara, chara, rams%pol(1, 1)%atom(i)%displacement(k)%XYZ(j)%frame(1), &
                                rams%pol(2, 2)%atom(i)%displacement(k)%XYZ(j)%frame(1), rams%pol(3, 3)%atom(i)%displacement(k)%XYZ(j)%frame(1)
                            READ (runit, *) chara, chara, chara, rams%pol(1, 2)%atom(i)%displacement(k)%XYZ(j)%frame(1), &
                                rams%pol(1, 3)%atom(i)%displacement(k)%XYZ(j)%frame(1), rams%pol(2, 3)%atom(i)%displacement(k)%XYZ(j)%frame(1)
                            READ (runit, *) chara, chara, chara, rams%pol(2, 1)%atom(i)%displacement(k)%XYZ(j)%frame(1), &
                                rams%pol(3, 1)%atom(i)%displacement(k)%XYZ(j)%frame(1), rams%pol(3, 2)%atom(i)%displacement(k)%XYZ(j)%frame(1)
                        END DO
                    END DO
                END DO
            END DO
995         CONTINUE
            CLOSE (runit)

        ELSEIF (dips%type_dipole=='berry') THEN
            OPEN (FILE=dips%dip_file, STATUS='old', ACTION='read', IOSTAT=stat, IOMSG=msg, NEWUNIT=runit)!Reading dipoles
            !Check if file exists
            CALL check_file_open(stat, msg, dips%dip_file)
            DO
                DO k = 1, 2
                    DO i = 1, sys%natom
                        DO j = 1, 3
                            READ (runit, *, END=994)
                            READ (runit, *)
                            READ (runit, *) chara, dips%static_dip(1)%atom(i)%displacement(k)%XYZ(j)%frame(1), &
                                dips%static_dip(2)%atom(i)%displacement(k)%XYZ(j)%frame(1), &
                                dips%static_dip(3)%atom(i)%displacement(k)%XYZ(j)%frame(1)
                        END DO
                    END DO
                END DO
            END DO
994         CONTINUE
            CLOSE (runit)

        END IF

    END SUBROUTINE read_static

!********************************************************************************************
!********************************************************************************************

    SUBROUTINE read_static_resraman(static_dip_file, static_dip_rtp, sys, rams)

        TYPE(systems), INTENT(INOUT)        :: sys
        TYPE(raman), INTENT(INOUT)        :: rams
        CHARACTER(LEN=40), INTENT(IN)                               :: static_dip_file
        TYPE(static_property), DIMENSION(3), INTENT(OUT)  :: static_dip_rtp

        CHARACTER(LEN=str_len)                                          :: chara, msg
        INTEGER                                                    :: x, y, i, j, k, m, stat, stat2, runit

        !! Allocate
        CALL rams%RR%init_rr_static_dip(static_dip_rtp, sys%natom, rams%RR%framecount_rtp + 1)

        OPEN (FILE=static_dip_file, STATUS='old', ACTION='read', IOSTAT=stat, IOMSG=msg, NEWUNIT=runit) !Reading polarizabilties
        !Check if file exists
        CALL check_file_open(stat, msg, static_dip_file)
        DO
            DO k = 1, 2
                DO i = 1, sys%natom
                    DO j = 1, 3
                        READ (runit, *, END=994)
                        READ (runit, *)
                        DO m = 1, rams%RR%framecount_rtp + 1
                            READ (runit, *) chara, static_dip_rtp(1)%atom(i)%displacement(k)%XYZ(j)%frame(m), &
                                static_dip_rtp(2)%atom(i)%displacement(k)%XYZ(j)%frame(m), &
                                static_dip_rtp(3)%atom(i)%displacement(k)%XYZ(j)%frame(m)

                        END DO
                    END DO
                END DO
            END DO
        END DO
994     CONTINUE
        CLOSE (runit)
    END SUBROUTINE read_static_resraman
END MODULE read_traj
