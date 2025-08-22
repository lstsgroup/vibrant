MODULE dipole_calc

   USE kinds, ONLY: dp
   USE constants, ONLY: debye
   USE vib_types, ONLY: global_settings, systems, molecular_dynamics, static, dipoles
   USE setup, ONLY: build_hexagonal_hmat, pbc_orthorombic, pbc_hexagonal, pbc_hexagonal, pbc_orthorombic_old, invert3x3, &
                      build_oblique_hmat

   IMPLICIT NONE

   PUBLIC :: wannier, center_mass, wannier_frag!, solv_frag_index ! wannier_frag_old,center_mass_old
   PRIVATE

CONTAINS

SUBROUTINE center_mass(filename, fragment, gs, sys, md, dips)

   USE kinds, ONLY : dp
   IMPLICIT NONE

   TYPE(global_settings),    INTENT(INOUT) :: gs
   TYPE(systems),            INTENT(INOUT) :: sys
   TYPE(molecular_dynamics), INTENT(INOUT) :: md
   TYPE(dipoles),            INTENT(INOUT) :: dips

   CHARACTER(LEN=*), INTENT(INOUT) :: filename
   INTEGER, DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT) :: fragment

   ! locals
   INTEGER :: m, j, stat
   REAL(dp) :: mass_tot
   REAL(dp) :: COM_cart(3)
REAL(dp) :: COC_frac(3), COC_cart(3), xfrac(3)
REAL(dp) :: charge_tot
   REAL(dp) :: hmat(3,3), h_inv(3,3)

   !-----------------------------------------------------------
   ! Allocate arrays
   !-----------------------------------------------------------
   ALLOCATE (sys%fragments%mass_tot_frag(sys%framecount, sys%mol_num))
   ALLOCATE (sys%fragments%refpoint(sys%framecount, 3))
   ALLOCATE (fragment(sys%framecount, sys%mol_num, 32))
   ALLOCATE (sys%fragments%natom_frag(sys%mol_num))

   sys%fragments%refpoint      = 0.0_dp
   sys%fragments%mass_tot_frag = 0.0_dp
   fragment                    = 0
   sys%fragments%natom_frag    = 0

   !-----------------------------------------------------------
!    Compute COM for each frame (simple Cartesian average)
   !-----------------------------------------------------------
   DO m = 1, sys%framecount
      COM_cart = 0.0_dp
      mass_tot = 0.0_dp

      DO j = 1, sys%natom
         COM_cart = COM_cart + md%coord_v(m,j,:) * sys%mass_atom(j)
         mass_tot = mass_tot + sys%mass_atom(j)
      END DO

      sys%fragments%refpoint(m,:) = COM_cart / mass_tot
   END DO

   !-----------------------------------------------------------
   ! Write COM for visualization
   !-----------------------------------------------------------
   OPEN (UNIT=12, FILE=TRIM(filename)//'-COM.xyz', STATUS='unknown', IOSTAT=stat)
   IF (stat /= 0) STOP "Error opening COM output file"

   DO m = 1, sys%framecount
      WRITE (12,*) sys%natom+1
      WRITE (12,*)
      DO j = 1, sys%natom
         WRITE (12,*) sys%element(j), md%coord_v(m,j,:)
      END DO
      WRITE (12,*) 'N', sys%fragments%refpoint(m,:)   ! COM marker
   END DO

   CLOSE (12)

END SUBROUTINE center_mass

!!***************************************************************************************************
SUBROUTINE wannier_frag(natom_frag, filename, dipole, fragment, gs, sys, md, dips)

   USE kinds, ONLY : dp
   IMPLICIT NONE

   TYPE(global_settings),    INTENT(INOUT) :: gs
   TYPE(systems),            INTENT(INOUT) :: sys
   TYPE(molecular_dynamics), INTENT(INOUT) :: md
   TYPE(dipoles),            INTENT(INOUT) :: dips

   CHARACTER(LEN=*), INTENT(IN) :: filename
   INTEGER, DIMENSION(:),      ALLOCATABLE, INTENT(INOUT) :: natom_frag
   INTEGER, DIMENSION(:,:,:),  ALLOCATABLE, INTENT(INOUT) :: fragment
   REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT)   :: dipole

   ! locals
   REAL(dp) :: hmat(3,3), h_inv(3,3)
   REAL(dp) :: COM_frac(3), COM_cart(3), xfrac(3)
REAL(dp) :: COC_frac(3), COC_cart(3)
   REAL(dp) :: mass_tot, charge_tot
   REAL(dp) :: r_atom(3), r_wc(3)
   REAL(dp) :: dr(3), dr_wc(3), dr_atom(3), dist, dmin, dist_prev, dist_nearest
   INTEGER  :: m, i, iwc, parent, stat, prev_parent, nearest
   REAL(dp), PARAMETER :: HYST = 0.20_dp   ! Å, hysteresis threshold
    REAL(kind=dp)                                                :: pox_x, pox_y, pox_z, pox_all
   INTEGER,  ALLOCATABLE :: wc_parent(:,:)   ! [frame, iWC]
   INTEGER,  ALLOCATABLE :: n_wc_per_atom(:),counter(:)   ! [frame, iWC]
INTEGER :: nat_tot
CHARACTER(LEN=20) :: frame_label

nat_tot = sys%natom   ! includes both atoms and Wannier centers
   !-----------------------------------------------------------
   ! Allocate outputs
   !-----------------------------------------------------------
   ALLOCATE(dipole(sys%framecount, sys%mol_num, 3))
   dipole = 0.0_dp
   !-----------------------------------------------------------
   ! Reassign Wannier center indices (per frame)
   !-----------------------------------------------------------
  
  ! DO m = 1, 1
     ! DO i = 1, sys%natom
      !counter = 0
        ! IF (TRIM(sys%element(i)) == 'X') CYCLE
        ! IF (TRIM(sys%element(i)) == 'H') CYCLE
          !  DO iwc = 1, sys%natom
          !   IF (TRIM(sys%element(iwc)) .NE. 'X') CYCLE
          !   CALL pbc_hexagonal(md%coord_v(m,i,:), md%coord_v(m,iwc,:), sys, dr)
           !       dist = SQRT(SUM(dr**2))
          !        IF (dist < 1.00) THEN
         !         counter = counter + 1
                !  IF (counter >= sys%wc_number(i)) CYCLE
        !             IF (ANY(iwc == fragment(m, :, :))) CYCLE
       !              fragment(m, i, n) = k
                   
                  
!-----------------------------------------------------------
! Allocate outputs
!-----------------------------------------------------------

! lattice matrices
CALL build_hexagonal_hmat(sys, hmat)
CALL invert3x3(hmat, h_inv)

!-----------------------------------------------------------
! Assign parent atom for each Wannier center (per frame)
!-----------------------------------------------------------
ALLOCATE(wc_parent(sys%framecount, sys%natom),counter(sys%natom))
wc_parent = -1

! fixed quota per element
ALLOCATE(n_wc_per_atom(sys%natom))
n_wc_per_atom = 0
DO i = 1, sys%natom
   SELECT CASE (TRIM(sys%element(i)))
   CASE ("H")
      n_wc_per_atom(i) = 0
   CASE ("B")
      n_wc_per_atom(i) = 1
   CASE ("C")
      n_wc_per_atom(i) = 2
   CASE ("N")
      n_wc_per_atom(i) = 3
   CASE ("O")
      n_wc_per_atom(i) = 4
   CASE DEFAULT
      n_wc_per_atom(i) = 0   ! safe fallback
   END SELECT
END DO

! loop over frames
DO m = 1, sys%framecount
   counter = 0   ! reset WC count per atom for this frame

   DO iwc = 1, sys%natom
      IF (TRIM(sys%element(iwc)) == 'X') THEN
         dmin   = HUGE(1.0_dp)
         parent = -1
         DO i = 1, sys%natom
            IF (TRIM(sys%element(i)) /= 'X') THEN
               ! skip if this parent already got its quota
               IF (counter(i) >= n_wc_per_atom(i)) CYCLE
               CALL pbc_hexagonal(md%coord_v(m,iwc,:), md%coord_v(m,i,:), sys, dr)
               dist = SQRT(SUM(dr**2))
               IF (dist < dmin) THEN
                  dmin   = dist
                  parent = i
               END IF
            END IF
         END DO
         wc_parent(m, iwc) = parent
         IF (parent > 0) counter(parent) = counter(parent) + 1
      END IF
   END DO
!WRITE(*,*) "Frame ", m, " WC assignment summary:"
!DO i = 1, sys%natom
 !  IF (n_wc_per_atom(i) > 0) THEN
 !     WRITE(*,'(A3,I6,2I6)') TRIM(sys%element(i)), i, counter(i), n_wc_per_atom(i)
 !  END IF
!END DO
END DO

! Print how many WCs were assigned in the first frame
!DO i = 1, sys%natom
  ! IF (n_wc_per_atom(i) > 0) THEN
  !    WRITE(*,'(A3,I6,2I6)') TRIM(sys%element(i)), i, wc_parent(1, i), wc_parent(5000,i), counter(i), "warning"
      !WRITE(*,'(A3,I6,2I6)') TRIM(sys%element(i)), i, n_wc_per_atom(i), counter(i), "warning"
 !  END IF
!END DO

   !-----------------------------------------------------------
   ! Assign parent atom for each Wannier center (per frame)
   !-----------------------------------------------------------
  ! wc_parent = -1

   !DO m = 1, sys%framecount
    !  DO iwc = 1, sys%natom
     !    IF (TRIM(sys%element(iwc)) == 'X') THEN

            ! ---- find nearest atom (this frame) ----
      !      nearest      = -1
      !      dist_nearest = HUGE(1.0_dp)
       !     DO i = 1, sys%natom
        !       IF (TRIM(sys%element(i)) /= 'X') THEN
         !         CALL pbc_hexagonal(md%coord_v(m,iwc,:), md%coord_v(m,i,:), sys, dr)
          !        dist = SQRT(SUM(dr**2))
           !       IF (dist < dist_nearest) THEN
            !         dist_nearest = dist
             !        nearest      = i
             !     END IF
             !  END IF
           ! END DO

            ! ---- apply hysteresis vs previous parent ----
           ! IF (m > 1) THEN
              ! prev_parent = wc_parent(m-1, iwc)
              ! IF (prev_parent > 0) THEN
                !  CALL pbc_hexagonal(md%coord_v(m,iwc,:), md%coord_v(m,prev_parent,:), sys, dr)
                !  dist_prev = SQRT(SUM(dr**2))

               !   IF (dist_nearest < dist_prev - HYST) THEN
               !      parent = nearest
              !    ELSE
             !        parent = prev_parent
            !      END IF
           !    ELSE
          !        parent = nearest
         !      END IF
        !    ELSE
       !        parent = nearest
       !     END IF

    !        wc_parent(m, iwc) = parent
   !      END IF
  !    END DO
 !  END DO


!-----------------------------------------------------------
! Physically shift Wannier centers next to their parent atoms
!-----------------------------------------------------------
DO m = 1, sys%framecount
   DO iwc = 1, sys%natom
      IF (TRIM(sys%element(iwc)) == 'X') THEN
         parent = wc_parent(m, iwc)
         IF (parent > 0) THEN
            ! get PBC minimum-image vector WC−parent
            CALL pbc_hexagonal(md%coord_v(m,iwc,:), md%coord_v(m,parent,:), sys, dr)
           ! print*, SQRT(SUM(sys%cell%vec(:)**2)), SQRT(SUM(dr**2)), md%coord_v(m,parent,:), md%coord_v(m,iwc,:)
            ! overwrite WC coordinate: parent position + minimum-image vector
           !  dr(3) = md%coord_v(m,iwc,3) - md%coord_v(m,parent,3)
            print*, md%coord_v(m,parent,:), dr, md%coord_v(m,iwc,:)
             md%coord_v(m,iwc,:) = md%coord_v(m,parent,:) + dr
           ! CALL pbc_hexagonal(md%coord_v(m,iwc,:), md%coord_v(m,parent,:), sys, dr)
            print*, md%coord_v(m,iwc,:)
         END IF
      END IF
   END DO
END DO

   OPEN(UNIT=69, FILE='vmd_dipole.xyz', STATUS='unknown', IOSTAT=stat)
   IF (stat /= 0) STOP "Error opening dipole output file"
   DO m = 1, sys%framecount
      WRITE(69,*) sys%natom
      WRITE(69,*) 
   DO i = 1, sys%natom
      WRITE(69,*) TRIM(sys%element(i)), md%coord_v(m,i,1:3)
   END DO
   END DO
   CLOSE(69)

   !-----------------------------------------------------------
   ! Loop over frames: COM and dipole
   !-----------------------------------------------------------

      pox_x = 40.0_dp
      pox_y = 40.0_dp
      pox_z = 40.0_dp
      pox_all = 40.0_dp

   DO m = 1, sys%framecount

      ! --- mass-weighted COM under PBC (fractional) ---
      COM_frac = 0.0_dp; mass_tot = 0.0_dp
      DO i = 1, sys%natom
         xfrac = MATMUL(h_inv, md%coord_v(m,i,:))
         xfrac = xfrac - FLOOR(xfrac)
         COM_frac = COM_frac + xfrac * REAL(sys%mass_atom(i), dp)
         mass_tot = mass_tot + REAL(sys%mass_atom(i), dp)
      END DO
      COM_frac = COM_frac / MAX(mass_tot, 1.0_dp)
      COM_frac = COM_frac - FLOOR(COM_frac)
      COM_cart = MATMUL(hmat, COM_frac)
      sys%fragments%refpoint(m,:) = COM_cart
      sys%fragments%refpoint(m,:) = 0

      ! --- dipole for this frame ---
      dipole(m,1,:) = 0.0_dp

      ! (1) atoms: q_i * (r_i - COM)
      DO i = 1, sys%natom
         IF (TRIM(sys%element(i)) /= 'X') THEN
    CALL pbc_orthorombic_old(md%coord_v(m,i,:), sys%fragments%refpoint(m, :), sys%cell%vec, sys%cell%vec_pbc, &
                                           pox_all, pox_x, pox_y, pox_z) !! IF COM is for whole supercel!!!
   !         CALL pbc_hexagonal(md%coord_v(m,i,:), sys%fragments%refpoint(m,:), sys, dr)
   !         dipole(m,1,:) = dipole(m,1,:) + sys%charge(i) * dr
             dipole(m, 1, :) = dipole(m, 1, :) + sys%cell%vec_pbc*1.889725989_dp*sys%charge(i)
         END IF
      END DO

      ! (2) Wannier centers: q_wc * [(WC - parent) + (parent - COM)]
      DO iwc = 1, sys%natom
         IF (TRIM(sys%element(iwc)) == 'X') THEN
            parent = wc_parent(m, iwc)
            IF (parent > 0) THEN
               r_atom = md%coord_v(m,parent,:)
               r_wc   = md%coord_v(m,iwc,:)
    CALL pbc_orthorombic_old(r_wc, r_atom, sys%cell%vec, dr_wc, &
                                           pox_all, pox_x, pox_y, pox_z) !! IF COM is for whole supercel!!!
    !           CALL pbc_hexagonal(r_wc,   r_atom,   sys, dr_wc)   ! WC−parent
    CALL pbc_orthorombic_old(r_atom, sys%fragments%refpoint(m, :), sys%cell%vec, dr_atom, &
                                           pox_all, pox_x, pox_y, pox_z) !! IF COM is for whole supercel!!!
     !          CALL pbc_hexagonal(r_atom, sys%fragments%refpoint(m,:), sys, dr_atom) ! parent−COM
               dipole(m,1,:) = dipole(m,1,:) + sys%charge(iwc) * (dr_wc + dr_atom)*1.889725989_dp
            END IF
         END IF
      END DO

      ! Debye
      dipole(m,1,:) = dipole(m,1,:) * 4.80320427_dp

   END DO

   !-----------------------------------------------------------
   ! Write result
   !-----------------------------------------------------------
   OPEN(UNIT=68, FILE=TRIM(filename)//'_dipole.xyz', STATUS='unknown', IOSTAT=stat)
   IF (stat /= 0) STOP "Error opening dipole output file"
   DO m = 1, sys%framecount
      WRITE(68,'(I8,3F20.10)') m, dipole(m,1,:)
   END DO
   CLOSE(68)

   !-----------------------------------------------------------
   ! Clean up
   !-----------------------------------------------------------
   DEALLOCATE(wc_parent)

END SUBROUTINE wannier_frag



!!***************************************************************************************************
   SUBROUTINE wannier_frag2(natom_frag, filename, dipole, fragment, gs, sys, md, dips)

      TYPE(global_settings), INTENT(INOUT)        :: gs
      TYPE(systems), INTENT(INOUT)                :: sys
      TYPE(molecular_dynamics), INTENT(INOUT)     :: md
      TYPE(dipoles), INTENT(INOUT)     :: dips

      CHARACTER(LEN=40), INTENT(IN)                                :: filename
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT)              :: natom_frag
      INTEGER, DIMENSION(:, :, :), ALLOCATABLE, INTENT(INOUT)          :: fragment
      REAL(kind=dp), DIMENSION(5000, 3)                              :: bec, bec_pbc
      REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE, INTENT(OUT)       :: dipole

      CHARACTER(LEN=40)                                           :: length, filename_dip
      INTEGER                                                     :: stat   ! error status of OPEN statements
      INTEGER                                                     :: i, j, l, n, m
      REAL(kind=dp)                                                :: dist!,mass_tot(20)
      REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE                     :: mass
      REAL(kind=dp), DIMENSION(:), ALLOCATABLE                       :: coord2
      REAL(kind=dp)                                                :: pox_x, pox_y, pox_z, pox_all
      REAL(kind=dp)                                                :: hmat(3, 3), sqrt3, acosa, asina, a

      ALLOCATE (dipole(sys%framecount, sys%mol_num, 3), coord2(1))

      dist = 1.2_dp
      dipole = 0.0_dp

      pox_x = 40.0_dp
      pox_y = 40.0_dp
      pox_z = 40.0_dp
      pox_all = 40.0_dp

      sqrt3 = 1.73205080756887729352744634_dp
      a = 0.5_dp*(sys%cell%box_x + sys%cell%box_y)
      acosa = 0.5_dp*a
      asina = sqrt3*acosa
      hmat(1, 1) = a; hmat(1, 2) = acosa; hmat(1, 3) = 0.0_dp
      hmat(2, 1) = 0.0_dp; hmat(2, 2) = asina; hmat(2, 3) = 0.0_dp
      hmat(3, 1) = 0.0_dp; hmat(3, 2) = 0.0_dp; hmat(3, 3) = sys%cell%box_z
      PRINT *, sys%mol_num, 'check mol num'


dipole = 0.0_dp

DO m = 1, sys%framecount


   !--- Loop over atoms ---
   DO i = 1, sys%natom
!      CALL pbc_hexagonal(md%coord_v(m,i,:), sys%fragments%refpoint(m,:), sys)

      ! Accumulate dipole (in e·Å)
      dipole(m,1,:) = dipole(m,1,:) + sys%cell%vec_pbc(:) * sys%charge(i)
   END DO

   !--- Convert to Debye (1 e·Å = 4.80320427 Debye) ---
   dipole(m,1,:) = dipole(m,1,:) * 4.80320427_dp

END DO


!      DO m = 1, sys%framecount
 !        DO i = 1, sys%natom
  !          CALL pbc_hexagonal(md%coord_v(m, i, :), sys%fragments%refpoint(m, :), sys)
   !           dipole(m, 1, :) = dipole(m,1,:) + sys%cell%vec_pbc(:)*1.889725989_dp*sys%charge(i)
            !   END IF
   !         END DO
   !   END DO
!print*,sys%charge(1)
    !  PRINT *, sys%fragments%mass_tot_cell, 'mass tot cell'

      DO m = 1, sys%framecount
            IF (dipole(m, 1, 3) > 10) THEN
 !            dipole(m, 1, 3) = dipole(m, 1, 3) -  (hmat(3, 3)*2)
        ! IF (sys%system == '2' .AND. dips%type_dipole == '1') THEN
         !   IF (dipole(m, 1, 1) > 120) THEN
        !       dipole(m, 1, 1) = dipole(m, 1, 1) - (hmat(1, 1)*3*1.889725989)
       !     END IF
!  dipole(m,1,1)=dipole(m,1,1)-(hmat(1,1)*3*1.889725989)
      !      dipole(m, 1, 2) = dipole(m, 1, 2) + (hmat(2, 2)*4*1.889725989)!+(hmat(1,2)*3*1.889725989)
     !       dipole(m, 1, :) = REAL((dipole(m, 1, :)*2.54174622741_dp)/sys%fragments%mass_tot_cell, kind=dp) !!Dividing by total mass, change it later
    !     ELSEIF (sys%system == '1') THEN
       !     DO i = 1, sys%mol_num !! define something for this
      !         dipole(m, 1, :) = REAL((dipole(m, 1, :)*2.54174622741_dp)/sys%mass_tot(m, i), kind=dp) !!Dividing by total mass, change it later
     !       END DO
         END IF
      END DO

!DO j=1,sys%framecount !! define something for this
! WRITE(filename_dip, '(i0)') j
! OPEN(UNIT=68,FILE='COF-1_dipoles-'//trim(filename_dip)//'.txt',STATUS='unknown',IOSTAT=stat)
! DO m=1,sys%framecount
!  WRITE(68,*) m,dipole(m,j,1:3)!,vec_pbc(3)
! ENDDO
!ENDDO
!CLOSE(68)

!dipole = dipole/sys%mass_tot

OPEN(UNIT=68,FILE='dipole_result_final.xyz',STATUS='unknown',IOSTAT=stat)
DO m=1,sys%framecount
!DO j=1,20
!dipole(m,1,1) = dipole(m,1,3)
!dipole(m,1,2) = dipole(m,1,3)
        WRITE(68,*) m,dipole(m, 1, :)
        !WRITE(68,*) "net dipole",SQRT(DOT_PRODUCT(dipole(m,1,:),dipole(m,1,:)))*2.54174622741_dp

 ENDDO
! ENDDO
CLOSE(68)

!OPEN(UNIT=61,FILE='dipole_result_vmd.xyz',STATUS='unknown',IOSTAT=stat)
!DO m=1,sys%framecount
! DO j=1,20
      ! WRITE(61,'(2X,A15,4X,I2,6X,3F20.12)') "center of mass", j, sys%fragments%refpoint(m,j,:)

      !WRITE(61,'(2X,A15,4X,I2,6X,3F20.12)') "com+net dipole", j,  sys%fragments%refpoint(m,j,:)+(dipole(m,j,:)*2.541746227414447_dp)

! ENDDO
! ENDDO
!CLOSE(61)

!OPEN(UNIT=51,FILE='dipole_result_vmd_final.xyz',STATUS='unknown',IOSTAT=stat)
!DO m=1,sys%framecount
!WRITE(51,*) sys%framecount
!WRITE(51,*)
! DO j=1,20
      ! WRITE(51,*) "draw arrow     ", "{",sys%fragments%refpoint(m,j,:)&
      !        ,"}      ","{",sys%fragments%refpoint(m,j,:)+(dipole(m,j,:)*2.541746227414447_dp),"}"
! ENDDO
!ENDDO
!CLOSE(51)
      DEALLOCATE (sys%fragments%mass_tot_frag)
      DEALLOCATE (sys%fragments%refpoint)
   END SUBROUTINE wannier_frag2
!
!!***************************************************************************************************
!!***************************************************************************************************

   SUBROUTINE wannier(filename, dip, sys, md)

      TYPE(systems), INTENT(INOUT)                :: sys
      TYPE(molecular_dynamics), INTENT(INOUT)     :: md
      CHARACTER(LEN=40), INTENT(INOUT)                                  :: filename
      REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE, INTENT(OUT)            :: dip

      INTEGER                                                          :: stat                            ! error status of OPEN statements
      INTEGER                                                          :: k, i, j, l, n
      REAL(kind=dp)                                                     :: dist, coord1(3)
      REAL, DIMENSION(:), ALLOCATABLE                                   :: dip_tot !--> in your case of 1 molecule the dimension is one, i.e., it is a simple real
      REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE                          :: refpoint2
      REAL, DIMENSION(:, :, :), ALLOCATABLE                               :: com

      !        dist = 1.2_dp
      !        coord1(:) = 0.00_dp
      !
      !        ALLOCATE (dip(sys%framecount, sys%mol_num, 3))
      !        ALLOCATE (dip_tot(sys%mol_num))
      !        ALLOCATE (refpoint2(sys%framecount, 3))
      !        ALLOCATE (com(sys%framecount, sys%natom, 3))
      !
      !        dip = 0.0_dp
      !        l = 0
      !
      !        IF (sys%periodic=='y' .OR. sys%periodic=='yes') THEN
      !
      !            DO i = 1, sys%natom
      !                IF (sys%element(i).NE."O") CYCLE
      !                l = l + 1
      !                DO j = 1, sys%natom
      !                    IF (sys%element(j)=="O") CYCLE
      !                    DO k = 1, sys%framecount
      !                        CALL pbc_orthorombic(md%coord_v(k, j, :), md%coord_v(k, i, :), sys%cell_type%vec, sys%cell_type%vec_pbc, sys%cell_type%box_all, sys%cell_type%box_x, sys%cell_type%box_y, sys%cell_type%box_z)
      !
      !                        IF (SQRT(DOT_PRODUCT(sys%cell%sys%cell_type%vec_pbc, sys%cell%sys%cell_type%vec_pbc))>dist) CYCLE
      !                        IF (sys%element(j)=="X") THEN
      !                            dip(k, l, :) = REAL(dip(k, l, :) + sys%cell_type%vec_pbc*1.889725989_dp*(-2.0_dp)/debye, kind=dp)
      !                        ELSEIF (sys%element(j)=="H") THEN
      !                            dip(k, l, :) = REAL(dip(k, l, :) + sys%cell_type%vec_pbc*1.889725989_dp*(1.0_dp)/debye, kind=dp)
      !                        END IF
      !                    END DO
      !                    !dip_tot(l)=SQRT(dipole2(l,1)**2+dipole2(l,2)**2+dipole2(l,3)**2)
      !                END DO
      !            END DO
      !
      !        ELSE IF (sys%periodic=='n' .OR. sys%periodic=='no') THEN
      !
      !            refpoint2 = 0.0_dp
      !
      !            DO j = 1, sys%natom
      !                DO k = 1, sys%framecount
      !
      !                    CALL pbc_orthorombic(md%coord_v(k, j, :), coord1, sys%cell_type%vec, sys%cell_type%vec_pbc, sys%cell_type%box_all, sys%cell_type%box_x, sys%cell_type%box_y, sys%cell_type%box_z)
      !
      !                    com(k, j, :) = sys%cell_type%vec_pbc*sys%mass_atom(j)
      !                    refpoint2(k, :) = refpoint2(k, :) + com(k, j, :)
      !                END DO
      !            END DO
      !
      !            refpoint2 = REAL(refpoint2/sys%mass_tot, kind=dp)
      !
      !            DO j = 1, sys%natom
      !                l = 1
      !                DO k = 1, sys%framecount
      !                    CALL pbc_orthorombic(md%coord_v(k, j, :), refpoint2(k, :), sys%cell_type%vec, sys%cell_type%vec_pbc, sys%cell_type%box_all, sys%cell_type%box_x, sys%cell_type%box_y, sys%cell_type%box_z)
      !
      !                    IF (sys%element(j)=="X") THEN
      !                        dip(k, l, :) = REAL(dip(k, l, :) + sys%cell_type%vec_pbc*1.889725989_dp*(-2.0_dp)/debye, kind=dp)
      !                    ELSEIF (sys%element(j)=="H") THEN
      !                        dip(k, l, :) = REAL(dip(k, l, :) + sys%cell_type%vec_pbc*1.889725989_dp*(1.0_dp)/debye, kind=dp)
      !                    ELSEIF (sys%element(j)=="O") THEN
      !                        dip(k, l, :) = REAL(dip(k, l, :) + sys%cell_type%vec_pbc*1.889725989_dp*(6.0_dp)/debye, kind=dp)
      !                    ELSEIF (sys%element(j)=="B") THEN
      !                        dip(k, l, :) = REAL(dip(k, l, :) + sys%cell_type%vec_pbc*1.889725989_dp*(3.0_dp)/debye, kind=dp)
      !                    END IF
      !                    !dip_tot(l)=SQRT(dipole2(l,1)**2+dipole2(l,2)**2+dipole2(l,3)**2)
      !                END DO
      !            END DO
      !        END IF
      !
      !        OPEN (UNIT=60, FILE='dipole_result.xyz', STATUS='unknown', IOSTAT=stat)
      !
      !        DO i = 1, sys%framecount
      !            WRITE (60, *) sys%framecount
      !            WRITE (60, *) sys%mol_num
      !
      !            DO j = 1, sys%mol_num
      !                WRITE (60, *) j, dip(i, j, :)
      !            END DO
      !        END DO
      !        CLOSE (60)
      !
      !        DEALLOCATE (dip_tot, com, refpoint2)
      !
   END SUBROUTINE wannier

END MODULE dipole_calc

