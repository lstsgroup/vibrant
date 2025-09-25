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

MODULE dipole_calc

    USE kinds, ONLY: dp, str_len
    USE constants, ONLY: debye, bohr2ang
    USE iso_fortran_env, ONLY: output_unit, error_unit
    USE vib_types, ONLY: global_settings, systems, molecular_dynamics, static, dipoles
    USE cell_types, ONLY: build_hmat, pbc, invert3x3, determinant3x3
    USE read_traj, ONLY: check_file_open

    IMPLICIT NONE

    PUBLIC :: compute_dipole
    PRIVATE

CONTAINS

    SUBROUTINE compute_dipole(dipole, sys, md)
        TYPE(systems), INTENT(INOUT) :: sys
        TYPE(molecular_dynamics), INTENT(IN) :: md
        REAL(dp), DIMENSION(:, :, :), ALLOCATABLE, INTENT(OUT) :: dipole  ! (nframes,1,3)

        CHARACTER(LEN=str_len)                                          :: msg
        INTEGER :: m, i, k, stat, runit
        REAL(dp) :: hmat(3, 3), h_inv(3, 3), ratio, n, target_branch
        REAL(dp) :: pol_quantum(3), mass_tot, COM(3), dr(3)

        ! --- allocate arrays
        ALLOCATE (dipole(sys%framecount, 1, 3))
        ALLOCATE (sys%fragments%refpoint(sys%framecount, 1, 3))

        dipole = 0.0_dp

        ! --- lattice
        CALL build_hmat(sys, hmat)
        CALL invert3x3(hmat, h_inv)

        ! --- compute dipole
        DO m = 1, sys%framecount
            dipole(m, 1, :) = 0.0_dp
            COM = 0.0_dp; mass_tot = 0.0_dp
            DO i = 1, sys%natom
                IF (TRIM(sys%element(i))/='X') THEN
                    COM = COM + md%coord_v(m, i, :)*REAL(sys%mass_atom(i), dp)
                    mass_tot = mass_tot + REAL(sys%mass_atom(i), dp)
                END IF
            END DO
            COM = COM/MAX(mass_tot, 1.0_dp)
            sys%fragments%refpoint(m, 1, :) = COM

            ! nuclei
            DO i = 1, sys%natom
                CALL pbc(md%coord_v(m, i, :), sys%fragments%refpoint(m, 1, :), sys, dr)
                dipole(m, 1, :) = dipole(m, 1, :) + sys%charge(i)*dr/bohr2ang
            END DO

            ! convert to Debye
            dipole(m, 1, :) = dipole(m, 1, :)/debye
        END DO

   !!!! Calculate the polarization quantum
        DO i = 1, 3
            pol_quantum(i) = SUM(hmat(i, :))/(bohr2ang*debye)  ! e·Bohr
        END DO

   !!Subtract multiples of polarization quantum
        DO i = 1, 3
            DO m = 2, sys%framecount
                dipole(m, 1, i) = dipole(m, 1, i) - NINT((dipole(m, 1, i) - dipole(m - 1, 1, i))/pol_quantum(i))*pol_quantum(i)
            END DO
        END DO

        ! --- write to file
       ! OPEN (UNIT=69, FILE='COF-1_refpoint.xyz', STATUS='unknown', IOSTAT=stat)
       ! DO m = 1, sys%framecount
         !   WRITE (69, *) sys%natom + 1
         !   WRITE (69, *)
         !   DO i = 1, sys%natom
         !       WRITE (69, *) sys%element(i), md%coord_v(m, i, :)
        !    END DO
       !     WRITE (69, *) "N", sys%fragments%refpoint(m, 1, :)
      !  END DO

        ! --- write to file
     !   OPEN (UNIT=68, FILE='dipole.xyz', STATUS='unknown', IOSTAT=stat)
     !   IF (stat/=0) STOP "Error opening dipole output file"
     !   DO m = 1, sys%framecount
     !       WRITE (68, '(I8,3F20.10)') m, dipole(m, 1, :)
    !    END DO
   !     CLOSE (68)

        DEALLOCATE (sys%fragments%refpoint)
    END SUBROUTINE compute_dipole

END MODULE dipole_calc

