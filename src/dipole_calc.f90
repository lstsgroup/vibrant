MODULE dipole_calc

   USE kinds, ONLY: dp
   USE constants, ONLY: bohr2ang, debye
   USE vib_types, ONLY: global_settings, systems, molecular_dynamics, static, dipoles
   USE setup, ONLY: build_hexagonal_hmat, pbc_orthorombic, pbc_hexagonal, pbc_hexagonal, pbc_orthorombic_old, invert3x3, &
                      build_oblique_hmat, build_triclinic_hmat, determinant3x3

   IMPLICIT NONE

   PUBLIC :: compute_dipole_unwrapped, wannier_frag!, solv_frag_index ! wannier_frag_old,center_mass_old
   PRIVATE

CONTAINS


SUBROUTINE compute_dipole_unwrapped(dipole, sys, md)
   TYPE(systems),            INTENT(INOUT) :: sys
   TYPE(molecular_dynamics), INTENT(IN) :: md
   REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT) :: dipole  ! (nframes,1,3)

   INTEGER :: nframes, natoms, nwc
   INTEGER :: m, i, iwc, k, stat
   REAL(dp) :: hmat(3,3), h_inv(3,3)
   REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: wc_frac, wc_unwrapped, wc_cart
   REAL(dp) :: delta, coord3(3), pol_quantum(3)
   REAL(dp) :: mass_tot, charge_tot
   REAL(dp) :: COM_frac(3), COM_cart(3), xfrac(3)
   REAL(dp) :: dr(3), dr_wc(3), dr_atom(3), dist, dmin, dist_prev, dist_nearest

   ! --- system sizes
   nframes = sys%framecount
   natoms  = sys%natom
   nwc     = COUNT(sys%element(:) == 'X')
   coord3 = 0.0
   ! --- allocate arrays
   ALLOCATE(dipole(nframes,1,3))
   ALLOCATE(wc_frac(nframes,natoms,3))
   ALLOCATE(wc_unwrapped(nframes,natoms,3))
   ALLOCATE(wc_cart(nframes,natoms,3))
   ALLOCATE (sys%fragments%refpoint1(sys%framecount, 3))
   dipole = 0.0_dp

   ! --- lattice
   CALL build_triclinic_hmat(sys, hmat)
   CALL invert3x3(hmat, h_inv)

   ! --- extract WCs in fractional coords
!   iwc = 0
   DO i = 1, natoms
   !   IF (TRIM(sys%element(i)) == 'X') THEN
  !       iwc = iwc + 1
         DO m = 1, nframes
            wc_frac(m,i,:) = MATMUL(h_inv, md%coord_v(m,i,:))
            wc_frac(m,i,:) = wc_frac(m,i,:) - FLOOR(wc_frac(m,i,:)) ! wrap into [0,1)
         END DO
    !  END IF
   END DO

   ! --- unwrap WCs across frames
  ! wc_unwrapped(1,:,:) = wc_frac(1,:,:)
   !DO m = 2, nframes
     ! DO iwc = 1, nwc
       !  DO k = 1, 3
       !     delta = wc_frac(m,iwc,k) - wc_frac(m-1,iwc,k)
       !     IF (delta >  0.5_dp) THEN
       !        wc_unwrapped(m,iwc,k) = wc_unwrapped(m-1,iwc,k) + (delta - 1.0_dp)
       !     ELSEIF (delta < -0.5_dp) THEN
       !        wc_unwrapped(m,iwc,k) = wc_unwrapped(m-1,iwc,k) + (delta + 1.0_dp)
      !      ELSE
     !          wc_unwrapped(m,iwc,k) = wc_unwrapped(m-1,iwc,k) + delta
    !        END IF
   !      END DO
  !    END DO
 !  END DO

!   wc_unwrapped(:,:,:) = wc_frac(:,:,:)
   ! --- convert WCs to Cartesian
   DO m = 1, nframes
      DO iwc = 1, natoms
         wc_cart(m,iwc,:) = MATMUL(hmat, wc_frac(m,iwc,:))
      END DO
   END DO

   ! --- compute dipole
   DO m = 1, nframes
      dipole(m,1,:) = 0.0_dp
      COM_frac = 0.0_dp; mass_tot = 0.0_dp
      DO i = 1, sys%natom
         IF (TRIM(sys%element(i)) /= 'X') THEN
          !  CALL pbc_hexagonal(wc_cart(m,i,:), coord3(:), sys, dr)
             COM_frac = COM_frac + wc_cart(m,i,:)*REAL(sys%mass_atom(i), dp)
             mass_tot = mass_tot + REAL(sys%mass_atom(i), dp)
         ENDIF
      END DO
      COM_frac = COM_frac / MAX(mass_tot, 1.0_dp)
      sys%fragments%refpoint1(m,:) = COM_frac

      ! nuclei
      DO i = 1, natoms
   !      IF (TRIM(sys%element(i)) /= 'X') THEN
            CALL pbc_hexagonal( wc_cart(m,i,:), sys%fragments%refpoint1(m,:), sys, dr)
           ! dipole(m,1,:) = dipole(m,1,:) + sys%charge(i) * sys%cell%vec(:) / bohr2ang
            dipole(m,1,:) = dipole(m,1,:) + sys%charge(i) * dr / bohr2ang
    !     END IF
      END DO

      ! Wannier centers (–2 e each)
  !    DO iwc = 1, nwc
  !       CALL pbc_hexagonal(wc_cart(m,iwc,:), sys%fragments%refpoint1(m,:), sys, dr)
 !        dipole(m,1,:) = dipole(m,1,:) - 2.0_dp * dr / bohr2ang
!      END DO

      ! convert to Debye
      dipole(m,1,:) = dipole(m,1,:) / debye
   END DO
 DO i = 1, 3
   ! pol_quantum(i) = SUM(hmat(i,:)) * 1.889725989_dp / debye  ! e·Bohr
    pol_quantum(i) = SUM(hmat(:,i)) * 1.889725989_dp / debye  ! e·Bohr
 END DO
print*, pol_quantum(:), NORM2(hmat(:,1)), NORM2(hmat(:,2)), NORM2(hmat(:,3)), NORM2(hmat(:,1))*1.889725989_dp/debye
print*, 1.889725989_dp/debye

   DO m = 1, sys%framecount
      DO i = 1, 3
     !    print*, ANINT(dipole(m,1,i)/pol_quantum(i))
     !    nshift = ANINT(delta / pol_quantum(i))
     !    IF (nshift /= 0) THEN
          !  dipole(m,1,3) = dipole(m,1,3) - ANINT(dipole(m,1,3)/hmat(3,3))*hmat(3,3)
          !  dipole(m,1,1) = dipole(m,1,1) - ANINT(dipole(m,1,1)/(hmat(1,1)+)*hmat(1,1)
          !  dipole(m,1,2) = dipole(m,1,2) - ANINT(dipole(m,1,2)/hmat(2,2))*hmat(2,2)
            dipole(m,1,i) = dipole(m,1,i) - ANINT(dipole(m,1,i)/pol_quantum(i))*pol_quantum(i)
    !     END IF
      END DO
   END DO
   
   ! --- write to file
   OPEN(UNIT=69, FILE='COF-1_refpoint.xyz', STATUS='unknown', IOSTAT=stat)
   DO m = 1, nframes
      WRITE(69,*) sys%natom+1
      WRITE(69,*)
      DO i = 1, sys%natom
      !    IF (sys%element(i) .NE. 'X') THEN
          WRITE(69,*) sys%element(i),  wc_cart(m,i,:)
     !  !   ENDIF
      ENDDO
     ! DO iwc = 1, nwc
      !    WRITE(69,*) 'X', wc_cart(m,iwc,:)
    !  ENDDO
          WRITE(69,*) "N", sys%fragments%refpoint1(m,:)
   END DO
   ! --- write to file
   OPEN(UNIT=68, FILE='dipole.xyz', STATUS='unknown', IOSTAT=stat)
   IF (stat /= 0) STOP "Error opening dipole output file"
   DO m = 1, nframes
      WRITE(68,'(I8,3F20.10)') m, dipole(m,1,:)
   END DO
   CLOSE(68)

   DEALLOCATE(wc_frac, wc_unwrapped, wc_cart)
END SUBROUTINE compute_dipole_unwrapped



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
print*,hmat
   DO m = 1, sys%framecount

      ! --- mass-weighted COM under PBC (fractional) ---
      COM_frac = 0.0_dp; mass_tot = 0.0_dp
      COC_frac = 0.0_dp; charge_tot = 0.0_dp
      DO i = 1, sys%natom
         xfrac = MATMUL(h_inv, md%coord_v(m,i,:))
         xfrac = xfrac - FLOOR(xfrac)
         IF (TRIM(sys%element(i)) /= 'X') THEN
             COM_frac = COM_frac +  REAL(sys%mass_atom(i), dp)
             !COM_frac = COM_frac + xfrac * REAL(sys%mass_atom(i), dp)
             mass_tot = mass_tot + REAL(sys%mass_atom(i), dp)
         ELSEIF (TRIM(sys%element(i)) == 'X') THEN
          !   COC_frac = COC_frac + xfrac * REAL(sys%charge(i), dp)
             COC_frac = COC_frac +  REAL(sys%charge(i), dp)
             charge_tot = charge_tot + REAL(sys%charge(i), dp)
         ENDIF
      END DO
      COM_frac = COM_frac / MAX(mass_tot, 1.0_dp)
    !  COC_frac = COC_frac / MAX(charge_tot, 1.0_dp)
  !    COM_frac = COM_frac - FLOOR(COM_frac)
   !   COC_frac = COC_frac - FLOOR(COC_frac)
   !   COM_cart = MATMUL(hmat, COM_frac)
   !   COC_cart = MATMUL(hmat, COC_frac)
   !   sys%fragments%refpoint1(m,:) = COM_cart
      sys%fragments%refpoint1(m,:) = COM_frac
     ! sys%fragments%refpoint1(m,:) = 0.0_dp
   !   sys%fragments%refpoint2(m,:) = COC_cart
!      sys%fragments%refpoint1(m,:) = 0.0_dp

      ! --- dipole for this frame ---
      dipole1(m,1,:) = 0.0_dp
      dipole2(m,1,:) = 0.0_dp

      DO i = 1, sys%natom
    !     IF (TRIM(sys%element(i)) /= 'X') THEN
            CALL pbc_hexagonal(md%coord_v(m,i,:), sys%fragments%refpoint1(m,:), sys, dr)
        !    dipole1(m,1,:) = dipole1(m,1,:) + sys%charge(i) * md%coord_v(m,i,:) * 1.889725989_dp 
            dipole1(m,1,:) = dipole1(m,1,:) + sys%charge(i) * dr * 1.889725989_dp 
           
            print*, sys%element(i), sys%cell%vec, dr
     !    ELSEIF (TRIM(sys%element(i)) == 'X') THEN
     !       CALL pbc_hexagonal(md%coord_v(m,i,:), sys%fragments%refpoint1(m,:), sys, dr)
     !       dipole2(m,1,:) = dipole2(m,1,:) + sys%charge(i) * dr * 1.889725989_dp
           ! print*, sys%element(i),sys%cell%vec, dr
       !  END IF
      END DO
   END DO
   !-----------------------------------------------------------
   ! Polarization quantum (scalar per axis)
   !-----------------------------------------------------------
   V = ABS(determinant3x3(hmat)) * (1.889725989_dp**3)    ! Å³ → Bohr³

DO i = 1, 3
   pol_quantum(i) = -2.0_dp * NORM2(hmat(:,i)) * 1.889725989_dp  ! e·Bohr
END DO

 !  DO i = 1, 3
 !    pol_quantum(i) = (-2.0_dp / V) * (NORM2(hmat(:,i)) * 1.889725989_dp) !* 2.54174622741_dp
!   END DO

   !-----------------------------------------------------------
   ! Apply continuity correction only to WCs
   !-----------------------------------------------------------
   DO m = 2, sys%framecount
      DO i = 1, 3
         delta = dipole1(m,1,i) - dipole1(m-1,1,i)
         nshift = ANINT(delta / pol_quantum(i))
         IF (nshift /= 0) THEN
          !  dipole1(m:,1,i) = dipole1(m:,1,i) - nshift * pol_quantum(i)
         END IF
      END DO
   END DO

   DO m = 2, sys%framecount
      DO i = 1, 3
         delta = dipole2(m,1,i) - dipole2(m-1,1,i)
         nshift = ANINT(delta / pol_quantum(i))
      !   print*, nshift
         IF (nshift /= 0) THEN
   !         dipole2(m:,1,i) = dipole2(m:,1,i) - nshift * pol_quantum(i)
         END IF
      END DO
   END DO

   !-----------------------------------------------------------
   !-----------------------------------------------------------
   ! Combine and convert to Debye
   !-----------------------------------------------------------
   DO m = 1, sys%framecount
      dipole(m,1,:) = dipole1(m,1,:) * 2.54174622741_dp
      !dipole(m,1,:) = ((dipole1(m,1,:)) + (dipole2(m,1,:))) * 2.54174622741_dp
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

END MODULE dipole_calc
