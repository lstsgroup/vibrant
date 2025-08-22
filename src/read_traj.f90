MODULE read_traj

   USE kinds, ONLY: dp
   USE vib_types, ONLY: global_settings, systems, static, dipoles, raman, molecular_dynamics

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: read_coord, read_coord_frame, read_normal_modes, read_static, read_static_resraman

CONTAINS
   SUBROUTINE read_coord(filename, gs, sys, dips, rams)

      TYPE(global_settings), INTENT(INOUT)   :: gs
      TYPE(systems), INTENT(INOUT)        :: sys
      TYPE(dipoles), optional        :: dips
      TYPE(raman), optional        :: rams
      CHARACTER(LEN=40), INTENT(IN)                               :: filename

      INTEGER                                                   :: i, j, stat

      sys%framecount = 0

      IF (gs%spectral_type%read_function /= 'MD-RR') THEN
         OPEN (UNIT=50, FILE=filename, STATUS='old', IOSTAT=stat)
         READ (50, *) sys%natom
         CLOSE (50)
      ELSEIF (gs%spectral_type%read_function == 'MD-RR') THEN
         sys%natom = rams%RR%framecount_rtp
      END IF

      ALLOCATE (sys%element(sys%natom), sys%coord(sys%natom, 3))

      OPEN (UNIT=51, FILE=filename, STATUS='old', IOSTAT=stat)
      DO
         READ (51, *, END=998)
         READ (51, *)
         sys%framecount = sys%framecount + 1
         DO i = 1, sys%natom
            READ (51, *) sys%element(i), sys%coord(i, 1), sys%coord(i, 2), sys%coord(i, 3)
         END DO
      END DO
998   CONTINUE
      CLOSE (51)

      IF (gs%spectral_type%read_function /= 'P') THEN
    !  IF (dips%type_dipole == 'berry' .OR. dips%type_dipole == 'dfpt') THEN !!gas phase
         sys%mol_num = 1
         ! ELSEIF ((sys%periodic=='n' .AND. sys%system=='1') .OR. dips%type_dipole=='wannier') THEN !!fragment approach
         !     sys%mol_num = 44 !20 !! fix later to 20
    !  END IF
      END IF
      PRINT *, sys%mol_num, 'mol num', sys%framecount

   END SUBROUTINE read_coord

!********************************************************************************************
!********************************************************************************************
   SUBROUTINE read_coord_frame(natom, filename, coord_v, sys)

      TYPE(systems), INTENT(INOUT)        :: sys
      CHARACTER(LEN=40), INTENT(IN)                               :: filename
      INTEGER, INTENT(INOUT)                               :: natom
      REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE, INTENT(OUT)      :: coord_v

      INTEGER                                                    :: i, j, stat

      ALLOCATE (coord_v(sys%framecount, natom, 3))
      OPEN (UNIT=52, FILE=filename, STATUS='old', IOSTAT=stat)
      DO
         DO j = 1, sys%framecount
            READ (52, *, END=999)
            READ (52, *)
            DO i = 1, natom
               READ (52, *) sys%element(i), coord_v(j, i, 1), coord_v(j, i, 2), coord_v(j, i, 3)
            END DO
         END DO
      END DO
999   CONTINUE
      CLOSE (52)

   END SUBROUTINE read_coord_frame
!********************************************************************************************
!********************************************************************************************

   SUBROUTINE read_normal_modes(gs, sys, stats)

      ! Variables of your derived types:
      TYPE(global_settings), INTENT(INOUT)   :: gs
      TYPE(systems), INTENT(INOUT)        :: sys
      TYPE(static), INTENT(INOUT)        :: stats

      CHARACTER(LEN=40)                                          :: chara
      INTEGER                                                    :: i, j, k, m, n
      INTEGER                                                    :: stat

      IF (gs%spectral_type%read_function == 'NMA' .OR. stats%diag_hessian == 'y') THEN
         ALLOCATE (stats%force(2, sys%natom, 3, sys%natom, 3))
         OPEN (UNIT=49, FILE=stats%force_file, STATUS='old', IOSTAT=stat) !Reading forces
         DO i = 1, 2
            DO j = 1, sys%natom
               DO m = 1, 3
                  DO n = 1, sys%natom
                     READ (49, *) stats%force(i, j, m, n, 1), stats%force(i, j, m, n, 2), stats%force(i, j, m, n, 3)
                  END DO
               END DO
            END DO
         END DO
         CLOSE (49)

      ELSEIF (stats%diag_hessian == 'n') THEN
         stats%nmodes = 0
         OPEN (UNIT=50, FILE=stats%normal_freq_file, STATUS='old', IOSTAT=stat) !Reading normal freqs/coords
         DO
            READ (50, *, END=998) chara
            stats%nmodes = stats%nmodes + 1
         END DO
998      CONTINUE
         CLOSE (50)

         ALLOCATE (stats%freq(stats%nmodes), stats%disp(stats%nmodes, sys%natom, 3))

         OPEN (UNIT=51, FILE=stats%normal_freq_file, STATUS='old', IOSTAT=stat) !Reading normal freqs/coords
         DO i = 1, stats%nmodes
            READ (51, *, END=997) stats%freq(i)
         END DO
997      CONTINUE
         CLOSE (51)

         OPEN (UNIT=51, FILE=stats%normal_displ_file, STATUS='old', IOSTAT=stat) !Reading normal freqs/coords
         DO i = 1, stats%nmodes
            DO j = 1, sys%natom
               READ (51, *, END=996) stats%disp(i, j, 1), stats%disp(i, j, 2), stats%disp(i, j, 3)
            END DO
         END DO
996      CONTINUE
         CLOSE (51)
      END IF

   END SUBROUTINE read_normal_modes

!********************************************************************************************
!********************************************************************************************

   SUBROUTINE read_static(dip_file, static_dip, gs, sys, dips, rams)
      TYPE(global_settings), INTENT(INOUT)   :: gs
      TYPE(systems), INTENT(INOUT)        :: sys
      TYPE(dipoles), INTENT(INOUT)        :: dips
      TYPE(raman), INTENT(INOUT)        :: rams
      CHARACTER(LEN=40), INTENT(IN)                               :: dip_file
      REAL(kind=dp), DIMENSION(:, :, :, :), ALLOCATABLE, INTENT(OUT)    :: static_dip

      CHARACTER(LEN=40)                                          :: chara
      INTEGER                                                    :: i, j, k, m, n
      INTEGER                                                    :: stat

      ALLOCATE (rams%pol(sys%natom, 3, 2, 3, 3), static_dip(sys%natom, 3, 2, 3))

      IF (dips%type_dipole == 'dfpt') THEN
         OPEN (UNIT=52, FILE=rams%static_pol_file, STATUS='old', IOSTAT=stat) !Reading polarizabilties
         DO
            DO k = 1, 2
               DO i = 1, sys%natom
                  DO j = 1, 3
                     READ (52, *, END=995)
                     READ (52, *)
                     READ (52, *)
                     READ (52, *)
                     READ (52, *)
                     READ (52, *)
                     READ (52, *) chara, chara, chara, rams%pol(i, j, k, 1, 1), &
                        rams%pol(i, j, k, 2, 2), rams%pol(i, j, k, 3, 3)
                     READ (52, *) chara, chara, chara, rams%pol(i, j, k, 1, 2), &
                        rams%pol(i, j, k, 1, 3), rams%pol(i, j, k, 2, 3)
                     READ (52, *) chara, chara, chara, rams%pol(i, j, k, 2, 1), &
                        rams%pol(i, j, k, 3, 1), rams%pol(i, j, k, 3, 2)
                  END DO
               END DO
            END DO
         END DO
995      CONTINUE
         CLOSE (52)

      ELSEIF (dips%type_dipole == 'berry') THEN
         OPEN (UNIT=53, FILE=dip_file, STATUS='old', IOSTAT=stat) !Reading dipoles
         DO
            DO k = 1, 2
               DO i = 1, sys%natom
                  DO j = 1, 3
                     READ (53, *, END=994)
                     READ (53, *)
                     READ (53, *) chara, static_dip(i, j, k, 1), static_dip(i, j, k, 2), static_dip(i, j, k, 3)
                  END DO
               END DO
            END DO
         END DO
994      CONTINUE
         CLOSE (53)

      END IF

   END SUBROUTINE read_static

!********************************************************************************************
!********************************************************************************************

   SUBROUTINE read_static_resraman(dip_file, static_dip_rtp, sys, rams)

      TYPE(systems), INTENT(INOUT)        :: sys
      TYPE(raman), INTENT(INOUT)        :: rams
      CHARACTER(LEN=40), INTENT(IN)                               :: dip_file
      REAL(kind=dp), DIMENSION(:, :, :, :, :), ALLOCATABLE, INTENT(OUT)  :: static_dip_rtp

      CHARACTER(LEN=40)                                          :: chara
      INTEGER                                                    :: i, j, k, m, stat

      ALLOCATE (static_dip_rtp(sys%natom, 3, 2, 3, rams%RR%framecount_rtp + 1))

      OPEN (UNIT=53, FILE=dip_file, STATUS='old', IOSTAT=stat) !Reading polarizabilties
      DO
         DO k = 1, 2
            DO i = 1, sys%natom
               DO j = 1, 3
                  READ (53, *, END=994)
                  READ (53, *)
                  DO m = 1, rams%RR%framecount_rtp + 1
                     READ (53, *) chara, static_dip_rtp(i, j, k, 1, m), static_dip_rtp(i, j, k, 2, m), &
                        static_dip_rtp(i, j, k, 3, m)
                  END DO
               END DO
            END DO
         END DO
      END DO
994   CONTINUE
      CLOSE (53)

   END SUBROUTINE read_static_resraman
END MODULE read_traj
