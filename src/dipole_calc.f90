MODULE dipole_calc

   USE kinds, ONLY: dp
   USE constants, ONLY: debye
   USE vib_types, ONLY: global_settings, systems, molecular_dynamics, static, dipoles
   USE setup, ONLY: build_hexagonal_hmat, pbc_orthorombic, pbc_hexagonal, pbc_hexagonal, pbc_orthorombic_old, invert3x3, &
                      build_oblique_hmat, build_triclinic_hmat, determinant3x3

   IMPLICIT NONE

   PUBLIC :: wannier, wannier_frag!, solv_frag_index ! wannier_frag_old,center_mass_old
   PRIVATE

CONTAINS

!!***************************************************************************************************
SUBROUTINE wannier_frag(natom_frag, filename, dipole, fragment, gs, sys, md, dips)

   USE kinds, ONLY : dp
   IMPLICIT NONE

   TYPE(global_settings),    INTENT(INOUT) :: gs
   TYPE(systems),            INTENT(INOUT) :: sys
   TYPE(molecular_dynamics), INTENT(INOUT) :: md
   TYPE(dipoles),            INTENT(INOUT) :: dips

   CHARACTER(LEN=40), INTENT(IN) :: filename
   INTEGER, DIMENSION(:),      ALLOCATABLE, INTENT(INOUT) :: natom_frag
   INTEGER, DIMENSION(:,:,:),  ALLOCATABLE, INTENT(INOUT) :: fragment
   REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT)   :: dipole

   ! locals
   REAL(dp) :: hmat(3,3), h_inv(3,3), V
   REAL(dp) :: COM_frac(3), COM_cart(3), xfrac(3)
   REAL(dp) :: COC_frac(3), COC_cart(3)
   REAL(dp) :: mass_tot, charge_tot
   REAL(dp) :: dr(3), dr_wc(3), dr_atom(3), dist, dmin, dist_prev, dist_nearest
   INTEGER  :: m, i, iwc, j, parent, stat, prev_parent, nearest
   INTEGER :: nshift, n_wc
   REAL(dp), DIMENSION(:,:,:), ALLOCATABLE    :: dipole1, dipole2
   REAL(dp) :: pol_quantum(3), delta

   !-----------------------------------------------------------
   ! Allocate outputs
   !-----------------------------------------------------------
   ALLOCATE(dipole(sys%framecount, 1, 3))
   ALLOCATE(dipole1(sys%framecount, 1, 3))
   ALLOCATE(dipole2(sys%framecount, 1, 3))
   ALLOCATE (sys%fragments%refpoint1(sys%framecount, 3), sys%fragments%refpoint2(sys%framecount, 3))

   dipole = 0.0_dp

 !  CALL build_oblique_hmat(sys, hmat)   ! or build_hexagonal_hmat
 !  CALL build_hexagonal_hmat(sys, hmat)   ! or build_hexagonal_hmat
   CALL build_triclinic_hmat (sys, hmat)
   CALL invert3x3(hmat, h_inv)

   DO m = 1, sys%framecount

      ! --- mass-weighted COM under PBC (fractional) ---
      COM_frac = 0.0_dp; mass_tot = 0.0_dp
      COC_frac = 0.0_dp; charge_tot = 0.0_dp
      DO i = 1, sys%natom
         xfrac = MATMUL(h_inv, md%coord_v(m,i,:))
         xfrac = xfrac - FLOOR(xfrac)
         IF (TRIM(sys%element(i)) /= 'X') THEN
             COM_frac = COM_frac + xfrac * REAL(sys%mass_atom(i), dp)
             mass_tot = mass_tot + REAL(sys%mass_atom(i), dp)
         ELSEIF (TRIM(sys%element(i)) == 'X') THEN
             COC_frac = COC_frac + xfrac * REAL(sys%charge(i), dp)
             charge_tot = charge_tot + REAL(sys%charge(i), dp)
         ENDIF
      END DO
      COM_frac = COM_frac / MAX(mass_tot, 1.0_dp)
      COC_frac = COC_frac / MAX(charge_tot, 1.0_dp)
      COM_frac = COM_frac - FLOOR(COM_frac)
      COC_frac = COC_frac - FLOOR(COC_frac)
      COM_cart = MATMUL(hmat, COM_frac)
      COC_cart = MATMUL(hmat, COC_frac)
      sys%fragments%refpoint1(m,:) = COM_cart
      sys%fragments%refpoint2(m,:) = COC_cart

      ! --- dipole for this frame ---
      dipole1(m,1,:) = 0.0_dp
      dipole2(m,1,:) = 0.0_dp

      DO i = 1, sys%natom
         IF (TRIM(sys%element(i)) /= 'X') THEN
            CALL pbc_hexagonal(md%coord_v(m,i,:), sys%fragments%refpoint1(m,:), sys, dr)
            dipole1(m,1,:) = dipole1(m,1,:) + sys%charge(i) * dr * 1.889725989_dp 

         ELSEIF (TRIM(sys%element(i)) == 'X') THEN
            CALL pbc_hexagonal(md%coord_v(m,i,:), sys%fragments%refpoint2(m,:), sys, dr)
            dipole2(m,1,:) = dipole2(m,1,:) + sys%charge(i) * dr * 1.889725989_dp
         END IF
      END DO

  !    dipole(m,1,:) = dipole1(m,1,:) + dipole2(m,1,:)
      ! Debye
  !    dipole(m,1,:) = dipole(m,1,:) * 2.54174622741_dp/mass_tot

   END DO
   !-----------------------------------------------------------
   ! Polarization quantum (scalar per axis)
   !-----------------------------------------------------------
   V = ABS(determinant3x3(hmat)) * (1.889725989_dp**3)    ! Å³ → Bohr³
   DO i = 1, 3
     pol_quantum(i) = (-2.0_dp / V) * (NORM2(hmat(:,i)) * 1.889725989_dp) !* 2.54174622741_dp
   END DO

   !-----------------------------------------------------------
   ! Apply continuity correction only to WCs
   !-----------------------------------------------------------
   DO m = 2, sys%framecount
      DO i = 1, 3
         delta = dipole1(m,1,i) - dipole1(m-1,1,i)
         nshift = ANINT(delta / pol_quantum(i))
         IF (nshift /= 0) THEN
            dipole1(m,1,i) = dipole1(m,1,i) + nshift * pol_quantum(i)
         END IF
      END DO
   END DO

   DO m = 2, sys%framecount
      DO i = 1, 3
         delta = dipole2(m,1,i) - dipole2(m-1,1,i)
         nshift = ANINT(delta / pol_quantum(i))
         IF (nshift /= 0) THEN
            dipole2(m,1,i) = dipole2(m,1,i) + nshift * pol_quantum(i)
         END IF
      END DO
   END DO

   !-----------------------------------------------------------
   !-----------------------------------------------------------
   ! Combine and convert to Debye
   !-----------------------------------------------------------
   DO m = 1, sys%framecount
   !   dipole(m,1,:) = (dipole1(m,1,:)) * 2.54174622741_dp
      dipole(m,1,:) = ((dipole1(m,1,:)/mass_tot) + (dipole2(m,1,:)/mass_tot)) * 2.54174622741_dp
   END DO

   !-----------------------------------------------------------
   ! Write result
   !-----------------------------------------------------------
   OPEN(UNIT=68, FILE=TRIM(filename)//'_dipole2.xyz', STATUS='unknown', IOSTAT=stat)
   IF (stat /= 0) STOP "Error opening dipole output file"
   DO m = 1, sys%framecount
      WRITE(68,'(I8,3F20.10)') m, dipole(m,1,:)
   END DO
   CLOSE(68)

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
     ! DEALLOCATE (sys%fragments%mass_tot_frag)
     ! DEALLOCATE (sys%fragments%refpoint)
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

