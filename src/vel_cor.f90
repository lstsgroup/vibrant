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

MODULE vel_cor

   USE dipole_calc, ONLY: center_mass
   USE kinds, ONLY: dp
   USE constants, ONLY: pi, ang, fs2s, at_u, bohr2ang
   USE vib_types, ONLY: global_settings, systems, molecular_dynamics, static, dipoles, raman

   IMPLICIT NONE
   PUBLIC :: cvv, cvv_iso, cvv_aniso, cvv_only_x, cvv_resraman

CONTAINS
   SUBROUTINE cvv(natom, coord_v, sys, gs, md)

      TYPE(systems), INTENT(INOUT)                :: sys
      TYPE(global_settings), INTENT(INOUT)                :: gs
      TYPE(molecular_dynamics), INTENT(INOUT)     :: md
      INTEGER, INTENT(INOUT)                                    :: natom
      REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE, INTENT(INOUT)  :: coord_v

      CHARACTER(LEN=40)                                        :: chara
      INTEGER                                                  :: stat, i, j, k, m, t0, t1, l
      INTEGER, DIMENSION(:), ALLOCATABLE                         :: norm
      REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE                  :: coord

      ALLOCATE (md%z(0:2*md%t_cor - 1), norm(0:2*md%t_cor - 1))

      norm = 0
      md%z = 0.0_dp
      k = 0
      j = 0

      ! IF (md%t_cor < 0) THEN
      !    t_cor
      ! ENDIF

      DO t0 = 1, sys%framecount - 2
         t1 = MIN(sys%framecount - 2, t0 + md%t_cor)
         DO j = 1, natom
            IF (sys%frag_type == '2') THEN
               k = j + 8
            ELSEIF (sys%frag_type == '3') THEN
               k = j + 20
            ELSE
               k = j
            END IF
            IF (sys%input_mass == 'y') THEN
               md%z(0:t1 - t0) = md%z(0:t1 - t0) + (coord_v(t0, k, 1)*coord_v(t0:t1, j, 1) + coord_v(t0, j, 2)* &
                                                    coord_v(t0:t1, j, 2) + coord_v(t0, j, 3)*coord_v(t0:t1, j, 3))*sys%mass_atom(j)
            ELSE
               DO m = 1, 3
                  md%z(0:t1 - t0) = md%z(0:t1 - t0) + coord_v(t0, k, m)*coord_v(t0:t1, k, m)
               END DO
            END IF
         END DO
         norm(0:t1 - t0) = norm(0:t1 - t0) + 1
      END DO

      md%z(:) = md%z(:)/norm(:) !!Normalization

      !! unit conversion of dipole autocorrelation function to debye^2/s^2
      IF (gs%spectral_type%read_function == 'MD-IR') THEN
         md%z(:) = md%z(:)/(fs2s*fs2s)
      !! unit conversion of velocity autocorrelation function to m^2/s^2
      ELSE IF (gs%spectral_type%read_function == 'P') THEN
         IF (sys%type_traj == 'pos') THEN
            md%z(:) = md%z(:)*ang*ang/(fs2s*fs2s)
         ELSE IF (sys%type_traj == 'vel') THEN
            md%z(:) = md%z(:)*ang*ang*bohr2ang*bohr2ang/(at_u*at_u)
         END IF
      END IF

      md%z(md%t_cor) = 0.0_dp
      DO i = 1, md%t_cor - 1
         md%z(md%t_cor + i) = md%z(md%t_cor - i) !!Data mirroring
      END DO

      DO i = 0, 2*md%t_cor - 1
         md%z(i) = md%z(i)*((COS(i/(md%t_cor - 1.0_dp)/2.0_dp*pi))**2) !!Hann Window function
      END DO

      OPEN (UNIT=61, FILE='result_cvv.txt', STATUS='unknown', IOSTAT=stat) !!Write output
      DO i = 0, 2*md%t_cor - 1
         WRITE (61, *) md%z(i)
      END DO
      CLOSE (61)

      DEALLOCATE (norm)

   END SUBROUTINE cvv

!
!!***********************************************************************************
!!***********************************************************************************

   SUBROUTINE cvv_iso(mol_num, z_iso, alpha_diff_x, alpha_diff_y, alpha_diff_z, sys, md)

      TYPE(systems), INTENT(INOUT)                :: sys
      TYPE(molecular_dynamics), INTENT(INOUT)     :: md
      INTEGER, INTENT(INOUT)                                    :: mol_num
      REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(OUT)        :: z_iso
      REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE, INTENT(INOUT)  :: alpha_diff_x, alpha_diff_y, alpha_diff_z

      INTEGER                                                  :: stat, i, j, k, t0, t1
      INTEGER, DIMENSION(:), ALLOCATABLE                         :: norm

      ALLOCATE (z_iso(0:2*md%t_cor - 1), norm(0:2*md%t_cor - 1))
      norm = 0
      z_iso = 0.0_dp
      DO t0 = 1, sys%framecount - 2
         t1 = MIN(sys%framecount - 2, t0 + md%t_cor)
         DO j = 1, mol_num
            IF (sys%frag_type == '2') THEN
               k = j + 8
            ELSEIF (sys%frag_type == '3') THEN
               k = j + 20
            ELSE
               k = j
            END IF
            z_iso(0:t1 - t0) = z_iso(0:t1 - t0) + (alpha_diff_x(t0, k, 1) + alpha_diff_y(t0, k, 2) + alpha_diff_z(t0, k, 3))* &
                               (alpha_diff_x(t0:t1, k, 1) + alpha_diff_y(t0:t1, k, 2) + alpha_diff_z(t0:t1, k, 3))
         END DO
         norm(0:t1 - t0) = norm(0:t1 - t0) + 1
      END DO

      PRINT *, norm(0:3), 'iso norm'
      PRINT *, z_iso(0:3), 'iso z'
      z_iso(:) = z_iso(:)/(norm(:)*9._dp)  !!Normalization
      !   z_iso(:) = z_iso(:)/(2.0_dp*pi)
      !z_iso(:)=z_iso(:)/mol_num

      !!Unit conversion of Debye^2/(E^2*fs^2) into Debye^2/(E^2*s^2)
      z_iso(:) = z_iso(:)/(fs2s*fs2s)

      DO i = 0, md%t_cor - 1
         z_iso(i) = z_iso(i)*((COS(i/(md%t_cor - 1.0_dp)/2.0_dp*pi))**2) !!Hann window function
      END DO

      z_iso(md%t_cor) = 0.0_dp
      DO i = 1, md%t_cor - 1
         z_iso(md%t_cor + i) = z_iso(md%t_cor - i) !!Data mirroring
      END DO

      OPEN (UNIT=61, FILE='result_cvv_iso.txt', STATUS='unknown', IOSTAT=stat)

      DO i = 0, 2*md%t_cor - 1
         WRITE (61, *) z_iso(i)
      END DO
      CLOSE (61)
      DEALLOCATE (norm)

   END SUBROUTINE cvv_iso

!***********************************************************************************
!***********************************************************************************

   SUBROUTINE cvv_aniso(mol_num, z_aniso, alpha_diff_x, alpha_diff_y, alpha_diff_z, sys, md)

      TYPE(systems), INTENT(INOUT)                :: sys
      TYPE(molecular_dynamics), INTENT(INOUT)     :: md
      INTEGER, INTENT(INOUT)                                    :: mol_num
      REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(OUT)        :: z_aniso
      REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE, INTENT(INOUT)  :: alpha_diff_x, alpha_diff_y, alpha_diff_z

      INTEGER                                                  :: stat, i, j, k, t0, t1
      INTEGER, DIMENSION(:), ALLOCATABLE                         :: norm

      ALLOCATE (z_aniso(0:2*md%t_cor - 1), norm(0:2*md%t_cor - 1))

      norm = 0
      z_aniso = 0.0_dp

      DO t0 = 1, sys%framecount - 2
         t1 = MIN(sys%framecount - 2, t0 + md%t_cor)
         DO j = 1, mol_num
         IF (sys%frag_type == '2') THEN
            k = j + 8
         ELSEIF (sys%frag_type == '3') THEN
            k = j + 20
         ELSE
            k = j
         END IF
         z_aniso(0:t1 - t0) = z_aniso(0:t1 - t0) + (alpha_diff_x(t0, k, 1) - alpha_diff_y(t0, k, 2)) &
                              *(alpha_diff_x(t0:t1, k, 1) - alpha_diff_y(t0:t1, k, 2))/2.0_dp
         z_aniso(0:t1 - t0) = z_aniso(0:t1 - t0) + (alpha_diff_y(t0, k, 2) - alpha_diff_z(t0, k, 3)) &
                              *(alpha_diff_y(t0:t1, k, 2) - alpha_diff_z(t0:t1, k, 3))/2.0_dp
         z_aniso(0:t1 - t0) = z_aniso(0:t1 - t0) + (alpha_diff_z(t0, k, 3) - alpha_diff_x(t0, k, 1)) &
                              *(alpha_diff_z(t0:t1, k, 3) - alpha_diff_x(t0:t1, k, 1))/2.0_dp
         z_aniso(0:t1 - t0) = z_aniso(0:t1 - t0) + (alpha_diff_x(t0, k, 2)*0.50_dp + alpha_diff_y(t0, k, 1)*0.50_dp) &
                              *(alpha_diff_x(t0:t1, k, 2)*0.50_dp + alpha_diff_y(t0:t1, k, 1)*0.50_dp)*3.0_dp
         z_aniso(0:t1 - t0) = z_aniso(0:t1 - t0) + (alpha_diff_y(t0, k, 3)*0.50_dp + alpha_diff_z(t0, k, 2)*0.50_dp) &
                              *(alpha_diff_y(t0:t1, k, 3)*0.50_dp + alpha_diff_z(t0:t1, k, 2)*0.50_dp)*3.0_dp
         z_aniso(0:t1 - t0) = z_aniso(0:t1 - t0) + (alpha_diff_z(t0, k, 1)*0.50_dp + alpha_diff_x(t0, k, 3)*0.50_dp) &
                              *(alpha_diff_z(t0:t1, k, 1)*0.50_dp + alpha_diff_x(t0:t1, k, 3)*0.50_dp)*3.0_dp
         END DO
         norm(0:t1 - t0) = norm(0:t1 - t0) + 1
      END DO

      PRINT *, norm(0:3), 'aniso norm'
      PRINT *, z_aniso(0:3), 'aniso z'
      z_aniso(:) = z_aniso(:)/norm(:)
      !  z_aniso(:) = z_aniso(:)/(2.0_dp*pi)
      !z_aniso(:)=REAL(z_aniso(:)/mol_num,kind=dp)

      !!Unit conversion of Debye^2/(E^2*fs^2) into C^4*s^2/kg^2
      z_aniso(:) = z_aniso(:)/(fs2s*fs2s)

      DO i = 0, md%t_cor - 1
         z_aniso(i) = z_aniso(i)*((COS(i/(md%t_cor - 1.0_dp)/2.0_dp*pi))**2) !!Hann window function
      END DO

      z_aniso(md%t_cor) = 0.0_dp
      DO i = 1, md%t_cor - 1
         z_aniso(md%t_cor + i) = z_aniso(md%t_cor - i) !!Data mirroring
      END DO

      OPEN (UNIT=61, FILE='result_cvv_aniso.txt', STATUS='unknown', IOSTAT=stat)

      DO i = 0, 2*md%t_cor - 1
         WRITE (61, *) z_aniso(i)
      END DO
      CLOSE (61)

      DEALLOCATE (norm)

   END SUBROUTINE cvv_aniso

!***********************************************************************************
!***********************************************************************************

   SUBROUTINE cvv_resraman(framecount, natom, dt, alpha_resraman_x_diff_re, alpha_resraman_y_diff_re, &
                           alpha_resraman_z_diff_re, alpha_resraman_x_diff_im, alpha_resraman_y_diff_im, &
                           alpha_resraman_z_diff_im, z_iso_resraman, z_aniso_resraman, md)

      TYPE(molecular_dynamics), INTENT(INOUT)     :: md
      INTEGER, INTENT(INOUT)                                    :: framecount, natom
      REAL(kind=dp), INTENT(IN)                                  :: dt
      REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE, INTENT(INOUT)  :: alpha_resraman_x_diff_re, alpha_resraman_y_diff_re
      REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE, INTENT(INOUT)  :: alpha_resraman_z_diff_re
      REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE, INTENT(INOUT)  :: alpha_resraman_x_diff_im, alpha_resraman_y_diff_im
      REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE, INTENT(INOUT)  :: alpha_resraman_z_diff_im
      COMPLEX(kind=dp), DIMENSION(:, :), ALLOCATABLE, INTENT(OUT)   :: z_iso_resraman, z_aniso_resraman

      INTEGER                                                  :: stat, i, j, k, m, t0, t1
      INTEGER, DIMENSION(:, :), ALLOCATABLE                       :: norm_iso, norm_aniso
      COMPLEX(kind=dp)                                          :: im_unit

!!!ISOTROPIC!!!

      ALLOCATE (z_iso_resraman(0:2*md%t_cor, natom - 1), norm_iso(0:2*md%t_cor, natom - 1))
      ALLOCATE (z_aniso_resraman(0:2*md%t_cor, natom - 1), norm_aniso(0:2*md%t_cor, natom - 1))

      im_unit = (0.0_dp, 1.0_dp)
      z_iso_resraman = (0.0_dp, 0.0_dp)
      z_aniso_resraman = (0.0_dp, 0.0_dp)

      framecount = framecount - 2

      norm_iso = 0.0_dp
      DO t0 = 2, framecount
         t1 = MIN(framecount, t0 + md%t_cor)
         DO k = 1, natom - 1
  !!RE*RE
            z_iso_resraman(0:t1 - t0, k) = z_iso_resraman(0:t1 - t0, k) &
                                           + (alpha_resraman_x_diff_re(t0, k, 1) + alpha_resraman_y_diff_re(t0, k, 2) &
                                              + alpha_resraman_z_diff_re(t0, k, 3)) &
                                           *(alpha_resraman_x_diff_re(t0:t1, k, 1) &
                                             + alpha_resraman_y_diff_re(t0:t1, k, 2) + alpha_resraman_z_diff_re(t0:t1, k, 3))
  !!IM*IM
            z_iso_resraman(0:t1 - t0, k) = z_iso_resraman(0:t1 - t0, k) &
                                           + (alpha_resraman_x_diff_im(t0, k, 1) + alpha_resraman_y_diff_im(t0, k, 2) + &
                                              alpha_resraman_z_diff_im(t0, k, 3)) &
                                           *(alpha_resraman_x_diff_im(t0:t1, k, 1) &
                                             + alpha_resraman_y_diff_im(t0:t1, k, 2) + alpha_resraman_z_diff_im(t0:t1, k, 3))
  !!RE*IM
            z_iso_resraman(0:t1 - t0, k) = z_iso_resraman(0:t1 - t0, k) &
                                           + ((alpha_resraman_x_diff_re(t0, k, 1) + alpha_resraman_y_diff_re(t0, k, 2) &
                                               + alpha_resraman_z_diff_re(t0, k, 3)) &
                                              *(alpha_resraman_x_diff_im(t0:t1, k, 1) &
                                                + alpha_resraman_y_diff_im(t0:t1, k, 2) &
                                                + alpha_resraman_z_diff_im(t0:t1, k, 3)))*im_unit
  !!IM*RE (SUBTRACT)
            z_iso_resraman(0:t1 - t0, k) = z_iso_resraman(0:t1 - t0, k) &
                                           - ((alpha_resraman_x_diff_im(t0, k, 1) + alpha_resraman_y_diff_im(t0, k, 2) &
                                               + alpha_resraman_z_diff_im(t0, k, 3)) &
                                              *(alpha_resraman_x_diff_re(t0:t1, k, 1) + alpha_resraman_y_diff_re(t0:t1, k, 2) &
                                                + alpha_resraman_z_diff_re(t0:t1, k, 3)))*im_unit

            norm_iso(0:t1 - t0, k) = norm_iso(0:t1 - t0, k) + 1.0_dp
         END DO
      END DO

      z_iso_resraman(:, :) = z_iso_resraman(:, :)/norm_iso(:, :)
      z_iso_resraman(:, :) = z_iso_resraman(:, :)/9._dp
      z_iso_resraman(:, :) = z_iso_resraman(:, :)/(2.0_dp*pi)

      DO i = 0, md%t_cor - 1
         DO j = 1, natom - 1
            z_iso_resraman(i, j) = z_iso_resraman(i, j)*((COS(i/(md%t_cor - 1.0_dp)/2.0_dp*3.14_dp))**2)
         END DO
      END DO

      DO i = 1, natom - 1
         z_iso_resraman(md%t_cor, i) = 0.0_dp
      END DO

      DO i = 1, md%t_cor - 1
         DO j = 1, natom - 1
            z_iso_resraman(md%t_cor + i, j) = z_iso_resraman(md%t_cor - i, j)
         END DO
      END DO

      OPEN (UNIT=61, FILE='result_cvv_iso_resraman.txt', STATUS='unknown', IOSTAT=stat)

      DO i = 0, 2*md%t_cor - 1
         DO j = 1, natom - 1
            WRITE (61, *) z_iso_resraman(i, j)
         END DO
      END DO
      CLOSE (61)
      DEALLOCATE (norm_iso)

!!!ANISOTROPIC!!!

      norm_aniso = 0.0_dp
      DO t0 = 2, framecount
         t1 = MIN(framecount, t0 + md%t_cor)
         DO k = 1, natom - 1
  !!RE*RE
            z_aniso_resraman(0:t1 - t0, k) = z_aniso_resraman(0:t1 - t0, k) + (alpha_resraman_x_diff_re(t0, k, 1) &
                                                                               - alpha_resraman_y_diff_re(t0, k, 2)) &
                                             *(alpha_resraman_x_diff_re(t0:t1, k, 1) &
                                               - alpha_resraman_y_diff_re(t0:t1, k, 2))/2.0_dp
            z_aniso_resraman(0:t1 - t0, k) = z_aniso_resraman(0:t1 - t0, k) + (alpha_resraman_y_diff_re(t0, k, 2) &
                                                                               - alpha_resraman_z_diff_re(t0, k, 3)) &
                                             *(alpha_resraman_y_diff_re(t0:t1, k, 2) &
                                               - alpha_resraman_z_diff_re(t0:t1, k, 3))/2.0_dp
            z_aniso_resraman(0:t1 - t0, k) = z_aniso_resraman(0:t1 - t0, k) + (alpha_resraman_z_diff_re(t0, k, 3) &
                                                                               - alpha_resraman_x_diff_re(t0, k, 1)) &
                                             *(alpha_resraman_z_diff_re(t0:t1, k, 3) &
                                               - alpha_resraman_x_diff_re(t0:t1, k, 1))/2.0_dp

            z_aniso_resraman(0:t1 - t0, k) = z_aniso_resraman(0:t1 - t0, k) &
                                             + (alpha_resraman_x_diff_re(t0, k, 2)*0.50_dp + &
                                                alpha_resraman_y_diff_re(t0, k, 1)*0.50_dp) &
                                             *(alpha_resraman_x_diff_re(t0:t1, k, 2)*0.50_dp + &
                                               alpha_resraman_y_diff_re(t0:t1, k, 1)*0.50_dp)*3.0_dp

            z_aniso_resraman(0:t1 - t0, k) = z_aniso_resraman(0:t1 - t0, k) &
                                             + (alpha_resraman_y_diff_re(t0, k, 3)*0.50_dp + &
                                                alpha_resraman_z_diff_re(t0, k, 2)*0.50_dp) &
                                             *(alpha_resraman_y_diff_re(t0:t1, k, 3)*0.50_dp + &
                                               alpha_resraman_z_diff_re(t0:t1, k, 2)*0.50_dp)*3.0_dp

            z_aniso_resraman(0:t1 - t0, k) = z_aniso_resraman(0:t1 - t0, k) &
                                             + (alpha_resraman_z_diff_re(t0, k, 1)*0.50_dp + &
                                                alpha_resraman_x_diff_re(t0, k, 3)*0.50_dp) &
                                             *(alpha_resraman_z_diff_re(t0:t1, k, 1)*0.50_dp + &
                                               alpha_resraman_x_diff_re(t0:t1, k, 3)*0.50_dp)*3.0_dp

  !!IM*IM
            z_aniso_resraman(0:t1 - t0, k) = z_aniso_resraman(0:t1 - t0, k) &
                                             + (alpha_resraman_x_diff_im(t0, k, 1) - alpha_resraman_y_diff_im(t0, k, 2)) &
                                             *(alpha_resraman_x_diff_im(t0:t1, k, 1) &
                                               - alpha_resraman_y_diff_im(t0:t1, k, 2))/2.0_dp
            z_aniso_resraman(0:t1 - t0, k) = z_aniso_resraman(0:t1 - t0, k) &
                                             + (alpha_resraman_y_diff_im(t0, k, 2) - alpha_resraman_z_diff_im(t0, k, 3)) &
                                             *(alpha_resraman_y_diff_im(t0:t1, k, 2) &
                                               - alpha_resraman_z_diff_im(t0:t1, k, 3))/2.0_dp
            z_aniso_resraman(0:t1 - t0, k) = z_aniso_resraman(0:t1 - t0, k) &
                                             + (alpha_resraman_z_diff_im(t0, k, 3) - alpha_resraman_x_diff_im(t0, k, 1)) &
                                             *(alpha_resraman_z_diff_im(t0:t1, k, 3) &
                                               - alpha_resraman_x_diff_im(t0:t1, k, 1))/2.0_dp

            z_aniso_resraman(0:t1 - t0, k) = z_aniso_resraman(0:t1 - t0, k) &
                                             + (alpha_resraman_x_diff_im(t0, k, 2)*0.50_dp + &
                                                alpha_resraman_y_diff_im(t0, k, 1)*0.50_dp) &
                                             *(alpha_resraman_x_diff_im(t0:t1, k, 2)*0.50_dp + &
                                               alpha_resraman_y_diff_im(t0:t1, k, 1)*0.50_dp)*3.0_dp

            z_aniso_resraman(0:t1 - t0, k) = z_aniso_resraman(0:t1 - t0, k) &
                                             + (alpha_resraman_y_diff_im(t0, k, 3)*0.50_dp + &
                                                alpha_resraman_z_diff_im(t0, k, 2)*0.50_dp) &
                                             *(alpha_resraman_y_diff_im(t0:t1, k, 3)*0.50_dp + &
                                               alpha_resraman_z_diff_im(t0:t1, k, 2)*0.50_dp)*3.0_dp

            z_aniso_resraman(0:t1 - t0, k) = z_aniso_resraman(0:t1 - t0, k) &
                                             + (alpha_resraman_z_diff_im(t0, k, 1)*0.50_dp + &
                                                alpha_resraman_x_diff_im(t0, k, 3)*0.50_dp) &
                                             *(alpha_resraman_z_diff_im(t0:t1, k, 1)*0.50_dp + &
                                               alpha_resraman_x_diff_im(t0:t1, k, 3)*0.50_dp)*3.0_dp

  !!RE*IM
            z_aniso_resraman(0:t1 - t0, k) = z_aniso_resraman(0:t1 - t0, k) &
                                             + ((alpha_resraman_x_diff_re(t0, k, 1) - alpha_resraman_y_diff_re(t0, k, 2)) &
                                                *(alpha_resraman_x_diff_im(t0:t1, k, 1) &
                                                  - alpha_resraman_y_diff_im(t0:t1, k, 2)))/2.0_dp*im_unit
            z_aniso_resraman(0:t1 - t0, k) = z_aniso_resraman(0:t1 - t0, k) &
                                             + ((alpha_resraman_y_diff_re(t0, k, 2) - alpha_resraman_z_diff_re(t0, k, 3)) &
                                                *(alpha_resraman_y_diff_im(t0:t1, k, 2) &
                                                  - alpha_resraman_z_diff_im(t0:t1, k, 3)))/2.0_dp*im_unit
            z_aniso_resraman(0:t1 - t0, k) = z_aniso_resraman(0:t1 - t0, k) &
                                             + ((alpha_resraman_z_diff_re(t0, k, 3) - alpha_resraman_x_diff_re(t0, k, 1)) &
                                                *(alpha_resraman_z_diff_im(t0:t1, k, 3) &
                                                  - alpha_resraman_x_diff_im(t0:t1, k, 1)))/2.0_dp*im_unit

            z_aniso_resraman(0:t1 - t0, k) = z_aniso_resraman(0:t1 - t0, k) &
                                             + ((alpha_resraman_x_diff_re(t0, k, 2)* &
                                                 0.50_dp + alpha_resraman_y_diff_re(t0, k, 1) &
                                                 *0.50_dp)*(alpha_resraman_x_diff_im(t0:t1, k, 2)*0.50_dp &
                                                            + alpha_resraman_y_diff_im(t0:t1, k, 1)*0.50_dp))*3.0_dp*im_unit
            z_aniso_resraman(0:t1 - t0, k) = z_aniso_resraman(0:t1 - t0, k) &
                                             + ((alpha_resraman_y_diff_re(t0, k, 3)* &
                                                 0.50_dp + alpha_resraman_z_diff_re(t0, k, 2) &
                                                 *0.50_dp)*(alpha_resraman_y_diff_im(t0:t1, k, 3)*0.50_dp &
                                                            + alpha_resraman_z_diff_im(t0:t1, k, 2)*0.50_dp))*3.0_dp*im_unit
            z_aniso_resraman(0:t1 - t0, k) = z_aniso_resraman(0:t1 - t0, k) &
                                             + ((alpha_resraman_z_diff_re(t0, k, 1)* &
                                                 0.50_dp + alpha_resraman_x_diff_re(t0, k, 3) &
                                                 *0.50_dp)*(alpha_resraman_z_diff_im(t0:t1, k, 1)*0.50_dp &
                                                            + alpha_resraman_x_diff_im(t0:t1, k, 3)*0.50_dp))*3.0_dp*im_unit

  !!IM*RE (SUBTRACT)
            z_aniso_resraman(0:t1 - t0, k) = z_aniso_resraman(0:t1 - t0, k) &
                                             - ((alpha_resraman_x_diff_im(t0, k, 1) - alpha_resraman_y_diff_im(t0, k, 2)) &
                                                *(alpha_resraman_x_diff_re(t0:t1, k, 1) &
                                                  - alpha_resraman_y_diff_re(t0:t1, k, 2)))/2.0_dp*im_unit
            z_aniso_resraman(0:t1 - t0, k) = z_aniso_resraman(0:t1 - t0, k) &
                                             - ((alpha_resraman_y_diff_im(t0, k, 2) - alpha_resraman_z_diff_im(t0, k, 3)) &
                                                *(alpha_resraman_y_diff_re(t0:t1, k, 2) &
                                                  - alpha_resraman_z_diff_re(t0:t1, k, 3)))/2.0_dp*im_unit
            z_aniso_resraman(0:t1 - t0, k) = z_aniso_resraman(0:t1 - t0, k) &
                                             - ((alpha_resraman_z_diff_im(t0, k, 3) - alpha_resraman_x_diff_im(t0, k, 1)) &
                                                *(alpha_resraman_z_diff_re(t0:t1, k, 3) &
                                                  - alpha_resraman_x_diff_re(t0:t1, k, 1)))/2.0_dp*im_unit
            z_aniso_resraman(0:t1 - t0, k) = z_aniso_resraman(0:t1 - t0, k) &
                                             - ((alpha_resraman_x_diff_im(t0, k, 2)*0.50_dp &
                                                 + alpha_resraman_y_diff_im(t0, k, 1) &
                                                 *0.50_dp)*(alpha_resraman_x_diff_re(t0:t1, k, 2)*0.50_dp &
                                                            + alpha_resraman_y_diff_re(t0:t1, k, 1)*0.50_dp))*3.0_dp*im_unit
            z_aniso_resraman(0:t1 - t0, k) = z_aniso_resraman(0:t1 - t0, k) &
                                             - ((alpha_resraman_y_diff_im(t0, k, 3)*0.50_dp &
                                                 + alpha_resraman_z_diff_im(t0, k, 2) &
                                                 *0.50_dp)*(alpha_resraman_y_diff_re(t0:t1, k, 3)*0.50_dp &
                                                            + alpha_resraman_z_diff_re(t0:t1, k, 2)*0.50_dp))*3.0_dp*im_unit
            z_aniso_resraman(0:t1 - t0, k) = z_aniso_resraman(0:t1 - t0, k) &
                                             - ((alpha_resraman_z_diff_im(t0, k, 1)*0.50_dp &
                                                 + alpha_resraman_x_diff_im(t0, k, 3)*0.50_dp) &
                                                *(alpha_resraman_z_diff_re(t0:t1, k, 1)*0.50_dp &
                                                  + alpha_resraman_x_diff_re(t0:t1, k, 3)*0.50_dp))*3.0_dp*im_unit

            norm_aniso(0:t1 - t0, k) = norm_aniso(0:t1 - t0, k) + 1.0_dp
         END DO
      END DO

      z_aniso_resraman(:, :) = REAL(z_aniso_resraman(:, :)/norm_aniso(:, :), kind=dp)
      z_aniso_resraman(:, :) = REAL(z_aniso_resraman(:, :)/(2.0_dp*pi), kind=dp)

      DO i = 0, md%t_cor - 1
         DO j = 1, natom - 1
            z_aniso_resraman(i, j) = z_aniso_resraman(i, j)*0.5_dp*(1 + COS(2.0_dp*3.14_dp*i/(2.0_dp*(md%t_cor - 1))))
         END DO
      END DO

      DO i = 1, natom - 1
         z_aniso_resraman(md%t_cor, i) = 0.0_dp
      END DO

      DO i = 1, md%t_cor - 1
         DO j = 1, natom - 1
            z_aniso_resraman(md%t_cor + i, j) = z_aniso_resraman(md%t_cor - i, j)
         END DO
      END DO

      OPEN (UNIT=62, FILE='result_cvv_aniso_resraman.txt', STATUS='unknown', IOSTAT=stat)

      DO i = 0, 2*md%t_cor - 1
         DO j = 1, natom - 1
            WRITE (62, *) z_aniso_resraman(i, j)
         END DO
      END DO
      CLOSE (62)
      DEALLOCATE (norm_aniso)

   END SUBROUTINE cvv_resraman

!***********************************************************************************
!***********************************************************************************

   SUBROUTINE cvv_only_x(mol_num, framecount, z_para, z_ortho, alpha_diff_x, &
                         alpha_diff_y, alpha_diff_z, direction, md)

      TYPE(molecular_dynamics), INTENT(INOUT)     :: md
      CHARACTER(LEN=40), INTENT(IN)                             :: direction
      INTEGER, INTENT(INOUT)                                    :: framecount, mol_num
      REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(OUT)        :: z_para, z_ortho
      REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE, INTENT(INOUT)  :: alpha_diff_x, alpha_diff_y, alpha_diff_z

      CHARACTER(LEN=40)                                        :: chara
      CHARACTER(LEN=2), DIMENSION(:), ALLOCATABLE                :: element
      INTEGER, DIMENSION(:), ALLOCATABLE                         :: norm
      INTEGER                                                  :: stat, i, j, k, m, t0, t1

      ALLOCATE (z_para(0:md%t_cor*2), norm(0:md%t_cor*2))
      ALLOCATE (z_ortho(0:md%t_cor*2))

      framecount = framecount - 2

      norm = 0.0_dp
      z_para = 0.0_dp
      z_ortho = 0.0_dp
      DO t0 = 1, framecount
         t1 = MIN(framecount, t0 + md%t_cor)
         DO k = 1, mol_num

            IF (direction == '1') THEN
               z_para(0:t1 - t0) = z_para(0:t1 - t0) + alpha_diff_x(t0, k, 1)*alpha_diff_x(t0:t1, k, 1)!*18.01468_dp
               z_ortho(0:t1 - t0) = z_ortho(0:t1 - t0) + alpha_diff_x(t0, k, 2)*alpha_diff_x(t0:t1, k, 2)!*18.01468_dp
            ELSE IF (direction == '2') THEN
               z_para(0:t1 - t0) = z_para(0:t1 - t0) + alpha_diff_y(t0, k, 2)*alpha_diff_y(t0:t1, k, 2)
               z_ortho(0:t1 - t0) = z_ortho(0:t1 - t0) + alpha_diff_y(t0, k, 3)*alpha_diff_y(t0:t1, k, 3)
            ELSEIF (direction == '3') THEN
               z_para(0:t1 - t0) = z_para(0:t1 - t0) + alpha_diff_z(t0, k, 3)*alpha_diff_z(t0:t1, k, 3)
               z_ortho(0:t1 - t0) = z_ortho(0:t1 - t0) + alpha_diff_z(t0, k, 1)*alpha_diff_z(t0:t1, k, 1)
            END IF

         END DO
         norm(0:t1 - t0) = norm(0:t1 - t0) + 1.0_dp
      END DO

      z_para(:) = z_para(:)/norm(:)
      z_para(:) = z_para(:)/mol_num
      z_para(:) = z_para(:)/(2.0_dp*pi)
      z_ortho(:) = z_ortho(:)/norm(:)
      z_ortho(:) = z_ortho(:)/mol_num
      z_ortho(:) = z_ortho(:)/(2.0_dp*pi)

      DO i = 0, md%t_cor - 1
         z_para(i) = z_para(i)*0.5_dp*(1 + COS(2.0_dp*3.14_dp*i/(2.0_dp*(md%t_cor - 1))))
         z_ortho(i) = z_ortho(i)*0.5_dp*(1 + COS(2.0_dp*3.14_dp*i/(2.0_dp*(md%t_cor - 1))))
      END DO

      z_para(md%t_cor) = 0.0_dp
      z_ortho(md%t_cor) = 0.0_dp

      DO i = 1, md%t_cor - 1
         z_para(md%t_cor + i) = z_para(md%t_cor - i)
         z_ortho(md%t_cor + i) = z_ortho(md%t_cor - i)
      END DO

      OPEN (UNIT=60, FILE='result_cvv_para.txt', STATUS='unknown', IOSTAT=stat)

      DO i = 0, 2*md%t_cor - 1
         WRITE (60, *) z_para(i)
      END DO
      CLOSE (60)

      OPEN (UNIT=61, FILE='result_cvv_ortho.txt', STATUS='unknown', IOSTAT=stat)

      DO i = 0, 2*md%t_cor - 1
         WRITE (61, *) z_ortho(i)
      END DO
      CLOSE (61)
      DEALLOCATE (norm)
   END SUBROUTINE cvv_only_x

END MODULE vel_cor
