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

MODULE setup

    USE kinds, ONLY: dp
    USE constants, ONLY: speed_light
    USE vib_types, ONLY: global_settings, systems
    USE iso_fortran_env, ONLY: output_unit, error_unit

    IMPLICIT NONE

    PRIVATE

    PUBLIC ::  masses_charges, conversion

CONTAINS


!*********************************************************************************************
!*********************************************************************************************
    SUBROUTINE masses_charges(gs, sys)

        TYPE(global_settings), INTENT(INOUT) :: gs
        TYPE(systems), INTENT(INOUT) :: sys

        INTEGER :: i
        REAL(dp), DIMENSION(:, :), ALLOCATABLE :: mat1, mat2
        CHARACTER(len=:), ALLOCATABLE :: s   ! normalized symbol

        ! Allocate outputs (if not already allocated elsewhere)
        IF (.NOT. ALLOCATED(sys%atom_mass_inv_sqrt)) ALLOCATE (sys%atom_mass_inv_sqrt(sys%natom))
        IF (.NOT. ALLOCATED(sys%mass_mat)) ALLOCATE (sys%mass_mat(sys%natom, sys%natom))
        IF (.NOT. ALLOCATED(sys%mass_atom)) ALLOCATE (sys%mass_atom(sys%natom))
        IF (.NOT. ALLOCATED(sys%charge)) ALLOCATE (sys%charge(sys%natom))

        ALLOCATE (mat1(sys%natom, 1), mat2(1, sys%natom))

        sys%mass_atom = 0.0_dp
        sys%charge = 0.0_dp
        sys%mass_tot = 0.0_dp

        DO i = 1, sys%natom
            ! ---Warning no normalize of element symbol is applied (i.e. correct formatting is expected: 'H', 'C', 'Be' etc.") ---
            s = sys%element(i)

            ! --- Element mapping: atomic mass (IUPACc-style table value)
            !     and "charge" = #e- in outermost (highest n) shell ---
            ! ----- Period 1 -----
            IF (s=='H') THEN
                sys%mass_atom(i) = 1.00784_dp; sys%charge(i) = 1.0_dp
            ELSEIF (s=='He') THEN
                sys%mass_atom(i) = 4.002602_dp; sys%charge(i) = 2.0_dp

                ! ----- Period 2 -----
            ELSEIF (s=='Li') THEN
                sys%mass_atom(i) = 6.94_dp; sys%charge(i) = 1.0_dp
            ELSEIF (s=='Be') THEN
                sys%mass_atom(i) = 9.0121831_dp; sys%charge(i) = 2.0_dp
            ELSEIF (s=='B') THEN
                sys%mass_atom(i) = 10.811_dp; sys%charge(i) = 3.0_dp
            ELSEIF (s=='C') THEN
                sys%mass_atom(i) = 12.011_dp; sys%charge(i) = 4.0_dp
            ELSEIF (s=='N') THEN
                sys%mass_atom(i) = 14.0067_dp; sys%charge(i) = 5.0_dp
            ELSEIF (s=='O') THEN
                sys%mass_atom(i) = 15.999_dp; sys%charge(i) = 6.0_dp
            ELSEIF (s=='F') THEN
                sys%mass_atom(i) = 18.998403_dp; sys%charge(i) = 7.0_dp
            ELSEIF (s=='Ne') THEN
                sys%mass_atom(i) = 20.1797_dp; sys%charge(i) = 8.0_dp

                ! ----- Period 3 -----
            ELSEIF (s=='Na') THEN
                sys%mass_atom(i) = 22.989769_dp; sys%charge(i) = 1.0_dp
            ELSEIF (s=='Mg') THEN
                sys%mass_atom(i) = 24.305_dp; sys%charge(i) = 2.0_dp
            ELSEIF (s=='Al') THEN
                sys%mass_atom(i) = 26.9815385_dp; sys%charge(i) = 3.0_dp
            ELSEIF (s=='Si') THEN
                sys%mass_atom(i) = 28.085_dp; sys%charge(i) = 4.0_dp
            ELSEIF (s=='P') THEN
                sys%mass_atom(i) = 30.973761_dp; sys%charge(i) = 5.0_dp
            ELSEIF (s=='S') THEN
                sys%mass_atom(i) = 32.06_dp; sys%charge(i) = 6.0_dp
            ELSEIF (s=='Cl') THEN
                sys%mass_atom(i) = 35.45_dp; sys%charge(i) = 7.0_dp
            ELSEIF (s=='Ar') THEN
                sys%mass_atom(i) = 39.948_dp; sys%charge(i) = 8.0_dp

                ! ----- Period 4 -----
            ELSEIF (s=='K') THEN
                sys%mass_atom(i) = 39.0983_dp; sys%charge(i) = 1.0_dp
            ELSEIF (s=='Ca') THEN
                sys%mass_atom(i) = 40.078_dp; sys%charge(i) = 2.0_dp
                ! 3d block: "charge" counts electrons only in highest n shell (4s)
            ELSEIF (s=='Sc') THEN
                sys%mass_atom(i) = 44.955908_dp; sys%charge(i) = 2.0_dp
            ELSEIF (s=='Ti') THEN
                sys%mass_atom(i) = 47.867_dp; sys%charge(i) = 2.0_dp
            ELSEIF (s=='V') THEN
                sys%mass_atom(i) = 50.9415_dp; sys%charge(i) = 2.0_dp
            ELSEIF (s=='Cr') THEN
                sys%mass_atom(i) = 51.9961_dp; sys%charge(i) = 1.0_dp   ! 4s1 3d5
            ELSEIF (s=='Mn') THEN
                sys%mass_atom(i) = 54.938044_dp; sys%charge(i) = 2.0_dp
            ELSEIF (s=='Fe') THEN
                sys%mass_atom(i) = 55.845_dp; sys%charge(i) = 2.0_dp
            ELSEIF (s=='Co') THEN
                sys%mass_atom(i) = 58.933194_dp; sys%charge(i) = 2.0_dp
            ELSEIF (s=='Ni') THEN
                sys%mass_atom(i) = 58.6934_dp; sys%charge(i) = 2.0_dp
            ELSEIF (s=='Cu') THEN
                sys%mass_atom(i) = 63.546_dp; sys%charge(i) = 1.0_dp   ! 4s1 3d10
            ELSEIF (s=='Zn') THEN
                sys%mass_atom(i) = 65.38_dp; sys%charge(i) = 2.0_dp
            ELSEIF (s=='Ga') THEN
                sys%mass_atom(i) = 69.723_dp; sys%charge(i) = 3.0_dp
            ELSEIF (s=='Ge') THEN
                sys%mass_atom(i) = 72.630_dp; sys%charge(i) = 4.0_dp
            ELSEIF (s=='As') THEN
                sys%mass_atom(i) = 74.921595_dp; sys%charge(i) = 5.0_dp
            ELSEIF (s=='Se') THEN
                sys%mass_atom(i) = 78.971_dp; sys%charge(i) = 6.0_dp
            ELSEIF (s=='Br') THEN
                sys%mass_atom(i) = 79.904_dp; sys%charge(i) = 7.0_dp
            ELSEIF (s=='Kr') THEN
                sys%mass_atom(i) = 83.798_dp; sys%charge(i) = 8.0_dp

                ! ----- Period 5 -----
            ELSEIF (s=='Rb') THEN
                sys%mass_atom(i) = 85.4678_dp; sys%charge(i) = 1.0_dp      ! 5s1
            ELSEIF (s=='Sr') THEN
                sys%mass_atom(i) = 87.62_dp; sys%charge(i) = 2.0_dp      ! 5s2
            ELSEIF (s=='Y') THEN
                sys%mass_atom(i) = 88.90584_dp; sys%charge(i) = 2.0_dp      ! 5s2
            ELSEIF (s=='Zr') THEN
                sys%mass_atom(i) = 91.224_dp; sys%charge(i) = 2.0_dp
            ELSEIF (s=='Nb') THEN
                sys%mass_atom(i) = 92.90637_dp; sys%charge(i) = 1.0_dp      ! 5s1 (exception)
            ELSEIF (s=='Mo') THEN
                sys%mass_atom(i) = 95.95_dp; sys%charge(i) = 1.0_dp      ! 5s1 (exception)
            ELSEIF (s=='Tc') THEN
                sys%mass_atom(i) = 98.0_dp; sys%charge(i) = 1.0_dp      ! 5s1 (radioactive)
            ELSEIF (s=='Ru') THEN
                sys%mass_atom(i) = 101.07_dp; sys%charge(i) = 1.0_dp      ! 5s1 (common gs)
            ELSEIF (s=='Rh') THEN
                sys%mass_atom(i) = 102.90550_dp; sys%charge(i) = 1.0_dp      ! 5s1
            ELSEIF (s=='Pd') THEN
                sys%mass_atom(i) = 106.42_dp; sys%charge(i) = 0.0_dp      ! 4d10 5s0 (exception)
            ELSEIF (s=='Ag') THEN
                sys%mass_atom(i) = 107.8682_dp; sys%charge(i) = 1.0_dp      ! 5s1
            ELSEIF (s=='Cd') THEN
                sys%mass_atom(i) = 112.414_dp; sys%charge(i) = 2.0_dp      ! 5s2
            ELSEIF (s=='In') THEN
                sys%mass_atom(i) = 114.818_dp; sys%charge(i) = 3.0_dp      ! 5s2 5p1
            ELSEIF (s=='Sn') THEN
                sys%mass_atom(i) = 118.71_dp; sys%charge(i) = 4.0_dp      ! 5s2 5p2
            ELSEIF (s=='Sb') THEN
                sys%mass_atom(i) = 121.760_dp; sys%charge(i) = 5.0_dp      ! 5s2 5p3
            ELSEIF (s=='Te') THEN
                sys%mass_atom(i) = 127.60_dp; sys%charge(i) = 6.0_dp      ! 5s2 5p4
            ELSEIF (s=='I') THEN
                sys%mass_atom(i) = 126.90447_dp; sys%charge(i) = 7.0_dp      ! 5s2 5p5
            ELSEIF (s=='Xe') THEN
                sys%mass_atom(i) = 131.293_dp; sys%charge(i) = 8.0_dp      ! 5s2 5p6

                ! ----- Period 6 -----
            ELSEIF (s=='Cs') THEN
                sys%mass_atom(i) = 132.90545_dp; sys%charge(i) = 1.0_dp      ! 6s1
            ELSEIF (s=='Ba') THEN
                sys%mass_atom(i) = 137.327_dp; sys%charge(i) = 2.0_dp      ! 6s2

                ! Lanthanides (outermost n=6; 6s2 → charge=2 for all La–Lu)
            ELSEIF (s=='La') THEN
                sys%mass_atom(i) = 138.90547_dp; sys%charge(i) = 2.0_dp
            ELSEIF (s=='Ce') THEN
                sys%mass_atom(i) = 140.116_dp; sys%charge(i) = 2.0_dp
            ELSEIF (s=='Pr') THEN
                sys%mass_atom(i) = 140.90766_dp; sys%charge(i) = 2.0_dp
            ELSEIF (s=='Nd') THEN
                sys%mass_atom(i) = 144.242_dp; sys%charge(i) = 2.0_dp
            ELSEIF (s=='Pm') THEN
                sys%mass_atom(i) = 145.0_dp; sys%charge(i) = 2.0_dp      ! radioactive
            ELSEIF (s=='Sm') THEN
                sys%mass_atom(i) = 150.36_dp; sys%charge(i) = 2.0_dp
            ELSEIF (s=='Eu') THEN
                sys%mass_atom(i) = 151.964_dp; sys%charge(i) = 2.0_dp
            ELSEIF (s=='Gd') THEN
                sys%mass_atom(i) = 157.25_dp; sys%charge(i) = 2.0_dp
            ELSEIF (s=='Tb') THEN
                sys%mass_atom(i) = 158.92535_dp; sys%charge(i) = 2.0_dp
            ELSEIF (s=='Dy') THEN
                sys%mass_atom(i) = 162.500_dp; sys%charge(i) = 2.0_dp
            ELSEIF (s=='Ho') THEN
                sys%mass_atom(i) = 164.93033_dp; sys%charge(i) = 2.0_dp
            ELSEIF (s=='Er') THEN
                sys%mass_atom(i) = 167.259_dp; sys%charge(i) = 2.0_dp
            ELSEIF (s=='Tm') THEN
                sys%mass_atom(i) = 168.93422_dp; sys%charge(i) = 2.0_dp
            ELSEIF (s=='Yb') THEN
                sys%mass_atom(i) = 173.045_dp; sys%charge(i) = 2.0_dp
            ELSEIF (s=='Lu') THEN
                sys%mass_atom(i) = 174.9668_dp; sys%charge(i) = 2.0_dp

                ! 6th-period transition/post-transition
            ELSEIF (s=='Hf') THEN
                sys%mass_atom(i) = 178.49_dp; sys%charge(i) = 2.0_dp      ! 6s2
            ELSEIF (s=='Ta') THEN
                sys%mass_atom(i) = 180.94788_dp; sys%charge(i) = 2.0_dp
            ELSEIF (s=='W') THEN
                sys%mass_atom(i) = 183.84_dp; sys%charge(i) = 2.0_dp
            ELSEIF (s=='Re') THEN
                sys%mass_atom(i) = 186.207_dp; sys%charge(i) = 2.0_dp
            ELSEIF (s=='Os') THEN
                sys%mass_atom(i) = 190.23_dp; sys%charge(i) = 2.0_dp
            ELSEIF (s=='Ir') THEN
                sys%mass_atom(i) = 192.217_dp; sys%charge(i) = 2.0_dp
            ELSEIF (s=='Pt') THEN
                sys%mass_atom(i) = 195.084_dp; sys%charge(i) = 1.0_dp      ! 6s1 (exception)
            ELSEIF (s=='Au') THEN
                sys%mass_atom(i) = 196.96657_dp; sys%charge(i) = 1.0_dp      ! 6s1 (exception)
            ELSEIF (s=='Hg') THEN
                sys%mass_atom(i) = 200.592_dp; sys%charge(i) = 2.0_dp      ! 6s2
            ELSEIF (s=='Tl') THEN
                sys%mass_atom(i) = 204.3835_dp; sys%charge(i) = 3.0_dp      ! 6s2 6p1
            ELSEIF (s=='Pb') THEN
                sys%mass_atom(i) = 207.2_dp; sys%charge(i) = 4.0_dp      ! 6s2 6p2
            ELSEIF (s=='Bi') THEN
                sys%mass_atom(i) = 208.98040_dp; sys%charge(i) = 5.0_dp      ! 6s2 6p3
            ELSEIF (s=='Po') THEN
                sys%mass_atom(i) = 209.0_dp; sys%charge(i) = 6.0_dp      ! radioactive
            ELSEIF (s=='At') THEN
                sys%mass_atom(i) = 210.0_dp; sys%charge(i) = 7.0_dp      ! radioactive
            ELSEIF (s=='Rn') THEN
                sys%mass_atom(i) = 222.0_dp; sys%charge(i) = 8.0_dp      ! radioactive
                ! ---Placeholder element ---
            ELSEIF (s=='X') THEN
                sys%mass_atom(i) = 0.0_dp
                sys%charge(i) = -2.0_dp
            ELSE
                ! Unknown symbol: mark and continue
                WRITE(error_unit,'(4X,"[WARN]  ",A,A)') 'Unknown element symbol, got ', s
                WRITE(error_unit,'(4X,"[WARN]  ",A,A)') 'Setting atomc mass and charge do -1.0'
                sys%mass_atom(i) = -1.0_dp
                sys%charge(i) = -1.0_dp
            END IF

            sys%mass_tot = sys%mass_tot + sys%mass_atom(i)
        END DO

        ! Compute inv sqrt masses; guard against zero/negative masses
        DO i = 1, sys%natom
            IF (sys%mass_atom(i)>0.0_dp) THEN
                sys%atom_mass_inv_sqrt(i) = SQRT(1.0_dp/sys%mass_atom(i))
            ELSE
                sys%atom_mass_inv_sqrt(i) = 0.0_dp
            END IF
        END DO

        ! Build outer-product mass matrix (i,j) = invsqrt(m_i) * invsqrt(m_j)
        mat1(:, 1) = sys%atom_mass_inv_sqrt(:)
        mat2(1, :) = sys%atom_mass_inv_sqrt(:)
        sys%mass_mat = MATMUL(mat1, mat2)

        DEALLOCATE (mat1, mat2)
    END SUBROUTINE masses_charges

!!*********************************************************************************************
!!*********************************************************************************************

    SUBROUTINE conversion(dt, freq_range, dt_rtp, freq_range_rtp, freq_res, sinc_const)

        REAL(kind=dp), INTENT(OUT)              :: freq_range, freq_range_rtp, freq_res, sinc_const
        REAL(kind=dp), INTENT(IN)               ::  dt, dt_rtp

        INTEGER                               :: stat   ! error status of OPEN statements
        INTEGER                               :: i, j, k

        freq_range = REAL((1.0_dp/(dt*1e-15))/speed_light, kind=dp)
        freq_range_rtp = REAL((1.0_dp/(dt_rtp*1e-15))/speed_light, kind=dp)

        ! freq_res = REAL(freq_range/(2.0_dp*md%t_cor), kind=dp)
        !sinc_const = freq_res*dt*1.883652d-4 !!for sinc function

    END SUBROUTINE conversion
END MODULE setup
