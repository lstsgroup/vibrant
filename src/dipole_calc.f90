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

    USE kinds, ONLY: dp
    USE constants, ONLY: debye
    USE vib_types, ONLY: global_settings, systems, molecular_dynamics, static, dipoles
    USE setup, ONLY: pbc_orthorombic, pbc_hexagonal, pbc_hexagonal, pbc_orthorombic_old

    IMPLICIT NONE

    PUBLIC :: wannier, center_mass, wannier_frag, solv_frag_index ! wannier_frag_old,center_mass_old
    PRIVATE

CONTAINS

    SUBROUTINE center_mass(filename, fragment, gs, sys, md, dips)

        TYPE(global_settings), INTENT(INOUT)        :: gs
        TYPE(systems), INTENT(INOUT)        :: sys
        TYPE(molecular_dynamics), INTENT(INOUT)        :: md
        TYPE(dipoles), INTENT(INOUT)        :: dips
        CHARACTER(LEN=40), INTENT(INOUT)                             :: filename
        INTEGER, DIMENSION(:, :, :), ALLOCATABLE, INTENT(OUT)            :: fragment

        CHARACTER(LEN=40)                                           :: length
        CHARACTER(LEN=2), DIMENSION(:, :), ALLOCATABLE                 :: element_shifted
        INTEGER                                                     :: stat   ! error status of OPEN statements
        INTEGER                                                     :: r, m, p, q, i, j, k, n, l, o, s, t
        REAL(kind=dp)                                                :: bec(5000, 3), bec_pbc(5000, 3)
        REAL(kind=dp)                                                :: hmat(3, 3), sqrt3, acosa, asina, a
        REAL(kind=dp)                                                :: pox_x, pox_y, pox_z, mass_tot
        REAL(kind=dp), DIMENSION(3)                                   :: coord3
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE                   :: com2, coord_shifted
        REAL(kind=dp), DIMENSION(:, :, :, :), ALLOCATABLE                 :: com

        sys%mol_num = 44
        ALLOCATE (com(sys%framecount, sys%mol_num, 35, 3))
        ALLOCATE (coord_shifted(sys%framecount, sys%natom, 3))
        ALLOCATE (com2(sys%framecount, sys%natom, 3))
        ALLOCATE (sys%fragments%mass_tot_frag(sys%framecount, sys%mol_num))
        ALLOCATE (sys%fragments%refpoint(sys%framecount, sys%mol_num, 3))
        ALLOCATE (fragment(sys%framecount, sys%mol_num, 32))
        ALLOCATE (sys%fragments%natom_frag(sys%mol_num))
        ALLOCATE (element_shifted(sys%framecount, sys%natom))

        coord3(1) = 0.00_dp
        coord3(2) = 0.00_dp
        coord3(3) = 0.00_dp
        sys%fragments%nfrag = 0
        fragment = 0
        sys%fragments%natom_frag = 0
        sys%fragments%refpoint = 0.00_dp
        com = 0.0_dp
        com2 = 0.0_dp
        sys%fragments%mass_tot_frag = 0.0_dp
        n = 0
        l = 0
        q = 0
        sqrt3 = 1.73205080756887729352744634_dp

        a = 0.5_dp*(sys%cell%box_x + sys%cell%box_y)
        acosa = 0.5_dp*a
        asina = sqrt3*acosa
        hmat(1, 1) = a; hmat(1, 2) = acosa; hmat(1, 3) = 0.0_dp
        hmat(2, 1) = 0.0_dp; hmat(2, 2) = asina; hmat(2, 3) = 0.0_dp
        hmat(3, 1) = 0.0_dp; hmat(3, 2) = 0.0_dp; hmat(3, 3) = sys%cell%box_z

        !Assign fragments for the B-O ring!
        n = 0
        p = 0
        DO m = 1, sys%framecount
            j = 0
            DO i = 1, sys%natom
                IF (sys%element(i).NE.'B') CYCLE
                IF (ANY(i==fragment(m, :, :))) CYCLE
                j = j + 1
                n = 0
                fragment(m, j, 1) = i
                DO k = 1, sys%natom
                    IF (i==k) CYCLE
                    IF (ANY(k==fragment(m, :, :))) CYCLE
                    IF (sys%element(k).NE.'O') CYCLE
                    CALL pbc_hexagonal(md%coord_v(m, i, :), md%coord_v(m, k, :), sys)
                    IF (SQRT(DOT_PRODUCT(sys%cell%vec_pbc, sys%cell%vec_pbc))<1.551835_dp) THEN
                        IF (ANY(k==fragment(m, :, :))) CYCLE
                        n = n + 1
                        fragment(m, j, n + 1) = k
                        DO l = 1, sys%natom
                            IF (k==l) CYCLE
                            IF (i==l) CYCLE
                            IF (sys%element(l).NE.'B') CYCLE
                            CALL pbc_hexagonal(md%coord_v(m, k, :), md%coord_v(m, l, :), sys)
                            IF (SQRT(DOT_PRODUCT(sys%cell%vec_pbc, sys%cell%vec_pbc))<1.551835_dp) THEN
                                IF (ANY(l==fragment(m, :, :))) CYCLE
                                n = n + 1
                                fragment(m, j, n + 1) = l
                                DO o = 1, sys%natom
                                    IF (l==o) CYCLE
                                    IF (i==o) CYCLE
                                    IF (k==o) CYCLE
                                    IF (sys%element(o).NE.'O') CYCLE
                                    CALL pbc_hexagonal(md%coord_v(m, l, :), md%coord_v(m, o, :), sys)
                                    IF (SQRT(DOT_PRODUCT(sys%cell%vec_pbc, sys%cell%vec_pbc))<1.551835_dp) THEN
                                        IF (ANY(o==fragment(m, :, :))) CYCLE
                                        n = n + 1
                                        fragment(m, j, n + 1) = o
                                    END IF
                                END DO
                            END IF
                        END DO
                    END IF
                END DO
            END DO
        END DO

        !Assign fragments for the phenyl rings!
        n = 0
        DO m = 1, sys%framecount
            j = 8
            DO i = 1, sys%natom
                IF (sys%element(i).NE.'H') CYCLE
                IF (ANY(i==fragment(m, :, :))) CYCLE
                j = j + 1
                n = 0
                fragment(m, j, 1) = i
                DO k = 1, sys%natom
                    IF (i==k) CYCLE
                    IF (sys%element(k).NE.'C') CYCLE
                    CALL pbc_hexagonal(md%coord_v(m, i, :), md%coord_v(m, k, :), sys)
                    IF (SQRT(DOT_PRODUCT(sys%cell%vec_pbc, sys%cell%vec_pbc))<1.31835_dp) THEN
                        IF (ANY(k==fragment(m, :, :))) CYCLE
                        n = n + 1
                        fragment(m, j, n + 1) = k
                        DO l = 1, sys%natom
                            IF (sys%element(l).NE.'C') CYCLE
                            CALL pbc_hexagonal(md%coord_v(m, k, :), md%coord_v(m, l, :), sys)
                            IF (SQRT(DOT_PRODUCT(sys%cell%vec_pbc, sys%cell%vec_pbc))<1.55835_dp) THEN
                                IF (ANY(l==fragment(m, :, :))) CYCLE
                                n = n + 1
                                fragment(m, j, n + 1) = l
                                DO o = 1, sys%natom
                                    IF (sys%element(o)=='X') CYCLE
                                    CALL pbc_hexagonal(md%coord_v(m, l, :), md%coord_v(m, o, :), sys)
                                    IF (SQRT(DOT_PRODUCT(sys%cell%vec_pbc, sys%cell%vec_pbc))<1.55835_dp) THEN
                                        IF (ANY(o==fragment(m, :, :))) CYCLE
                                        n = n + 1
                                        fragment(m, j, n + 1) = o
                                        DO p = 1, sys%natom
                                            IF (sys%element(p)=='X') CYCLE
                                            CALL pbc_hexagonal(md%coord_v(m, o, :), md%coord_v(m, p, :), sys)
                                            IF (SQRT(DOT_PRODUCT(sys%cell%vec_pbc, sys%cell%vec_pbc))<1.55835_dp) THEN
                                                IF (ANY(p==fragment(m, :, :))) CYCLE
                                                n = n + 1
                                                fragment(m, j, n + 1) = p
                                                DO r = 1, sys%natom
                                                    IF (sys%element(r)=='X') CYCLE
                                                    CALL pbc_hexagonal(md%coord_v(m, p, :), md%coord_v(m, r, :), sys)
                                                    IF (SQRT(DOT_PRODUCT(sys%cell%vec_pbc, sys%cell%vec_pbc))<1.55835_dp) THEN
                                                        IF (ANY(r==fragment(m, :, :))) CYCLE
                                                        n = n + 1
                                                        fragment(m, j, n + 1) = r
                                                        DO s = 1, sys%natom
                                                            IF (sys%element(s)=='X') CYCLE
                                                            CALL pbc_hexagonal(md%coord_v(m, r, :), md%coord_v(m, s, :), sys)
                                                            IF (SQRT(DOT_PRODUCT(sys%cell%vec_pbc, sys%cell%vec_pbc))<1.55835_dp) THEN
                                                                IF (ANY(s==fragment(m, :, :))) CYCLE
                                                                n = n + 1
                                                                fragment(m, j, n + 1) = s
                                                                DO t = 1, sys%natom
                                                                    IF (sys%element(t)=='X') CYCLE
                                                                    CALL pbc_hexagonal(md%coord_v(m, s, :), md%coord_v(m, t, :), sys)
                                                                    IF (SQRT(DOT_PRODUCT(sys%cell%vec_pbc, sys%cell%vec_pbc))<1.55835_dp) THEN
                                                                        IF (ANY(t==fragment(m, :, :))) CYCLE
                                                                        n = n + 1
                                                                        fragment(m, j, n + 1) = t
                                                                    END IF
                                                                END DO
                                                            END IF
                                                        END DO
                                                    END IF
                                                END DO
                                            END IF
                                        END DO
                                    END IF
                                END DO
                            END IF
                        END DO
                    END IF
                END DO
            END DO
        END DO

    !!Assign fragments for the B-C bonds
        DO m = 1, sys%framecount
            j = 20
            DO i = 1, sys%natom
                IF (sys%element(i).NE.'B') CYCLE
                n = 0
                j = j + 1
                fragment(m, j, 1) = i
                DO k = 1, sys%natom
                    IF (sys%element(k).NE.'C') CYCLE
                    CALL pbc_hexagonal(md%coord_v(m, i, :), md%coord_v(m, k, :), sys)
                    IF (SQRT(DOT_PRODUCT(sys%cell%vec_pbc, sys%cell%vec_pbc))<2.0_dp) THEN
                        n = n + 1
                        fragment(m, j, n + 1) = k
                    END IF
                END DO
            END DO
        END DO

    !!! Assign wannier centers!!!

        DO m = 1, sys%framecount
            DO i = 1, 8
                n = 0
                DO j = 1, 6
                    IF (sys%element(fragment(m, i, j)).NE.'O') CYCLE
                    DO k = 1, sys%natom
                        IF (j==k) CYCLE
                        IF (sys%element(k).NE.'X') CYCLE
                        CALL pbc_hexagonal(md%coord_v(m, fragment(m, i, j), :), md%coord_v(m, k, :), sys)
                        IF (SQRT(DOT_PRODUCT(sys%cell%vec_pbc, sys%cell%vec_pbc))<0.7_dp) THEN
                            IF (ANY(k==fragment(m, :, :))) CYCLE
                            n = n + 1
                            fragment(m, i, n + 6) = k
                        END IF
                    END DO
                END DO
            END DO
        END DO

        l = 0

        DO m = 1, sys%framecount
            DO i = 9, 20
                n = 0
                DO j = 1, 10
                    IF (sys%element(fragment(m, i, j)).NE.'C') CYCLE
                    outer: DO k = 1, sys%natom
                        IF (sys%element(k).NE.'X') CYCLE
                        CALL pbc_hexagonal(md%coord_v(m, fragment(m, i, j), :), md%coord_v(m, k, :), sys)
                        IF (SQRT(DOT_PRODUCT(sys%cell%vec_pbc, sys%cell%vec_pbc))<1.0_dp) THEN
                            IF (sys%system=='1') THEN
                                inner: DO l = 1, sys%natom
                                    IF (sys%element(l).NE.'B') CYCLE
                                    CALL pbc_hexagonal(md%coord_v(m, k, :), md%coord_v(m, l, :), sys)
                                    IF (SQRT(DOT_PRODUCT(sys%cell%vec_pbc, sys%cell%vec_pbc))<1.1_dp) CYCLE outer
                                END DO inner
                            END IF
                            IF (ANY(k==fragment(m, :, :))) CYCLE
                            n = n + 1
                            fragment(m, i, n + 10) = k
                        END IF
                    END DO outer
                END DO
            END DO
        END DO

        DO m = 1, sys%framecount
            DO i = 21, sys%mol_num
                DO j = 1, 2
                    IF (sys%element(fragment(m, i, j)).NE.'C') CYCLE
                    DO k = 1, sys%natom
                        IF (sys%element(k).NE.'X') CYCLE
                        IF (ANY(k==fragment(m, :, :))) CYCLE
                        CALL pbc_hexagonal(md%coord_v(m, fragment(m, i, j), :), md%coord_v(m, k, :), sys)
                        IF (SQRT(DOT_PRODUCT(sys%cell%vec_pbc, sys%cell%vec_pbc))<1.0_dp) THEN
                            fragment(m, i, 3) = k
                        END IF
                    END DO
                END DO
            END DO
        END DO

    !!get sys%fragments%natom_frag
        DO i = 1, sys%framecount
            DO j = 1, sys%mol_num
                sys%fragments%natom_frag(j) = COUNT(fragment(i, j, :).NE.0)
            END DO
        END DO

        IF (sys%frag_type=='1') THEN
            sys%fragments%nfrag = 8
        ELSEIF (sys%frag_type=='2') THEN
            sys%fragments%nfrag = 12
        ELSEIF (sys%frag_type=='3') THEN
            sys%fragments%nfrag = 24
        ELSEIF (sys%system=='2' .AND. dips%type_dipole=='1') THEN
            sys%fragments%nfrag = 1
        END IF

    !!get total mass
        DO m = 1, sys%framecount
            DO i = 1, sys%mol_num
                DO j = 1, sys%fragments%natom_frag(i)
                    sys%fragments%mass_tot_frag(m, i) = sys%fragments%mass_tot_frag(m, i) + sys%mass_atom(fragment(m, i, j))
                END DO
            END DO
        END DO

    !!Find center of mass

    !!!For B-O and C-B fragments!!!
        IF (sys%frag_type=='1' .OR. sys%frag_type=='2' .OR. (dips%type_dipole=='1' .AND. sys%system=='2')) THEN
            DO m = 1, sys%framecount
                DO i = 1, 20
                    DO j = 2, sys%fragments%natom_frag(i)
                        CALL pbc_hexagonal(md%coord_v(m, fragment(m, i, j), :), md%coord_v(m, fragment(m, i, 1), :), sys)
                        IF (sys%cell%vec(1)>3.0_dp .AND. sys%cell%vec(2)>5.2_dp) THEN
                            md%coord_v(m, fragment(m, i, j), 1) = md%coord_v(m, fragment(m, i, j), 1) - hmat(1, 2)
                            md%coord_v(m, fragment(m, i, j), 2) = md%coord_v(m, fragment(m, i, j), 2) - hmat(2, 2)
                            CALL pbc_hexagonal(md%coord_v(m, fragment(m, i, j), :), md%coord_v(m, fragment(m, i, 1), :), sys)
                        END IF
                        IF (sys%cell%vec(1)<-3.0_dp .AND. sys%cell%vec(2)>5.2_dp) THEN
                            md%coord_v(m, fragment(m, i, j), 1) = md%coord_v(m, fragment(m, i, j), 1) + hmat(1, 2)
                            md%coord_v(m, fragment(m, i, j), 2) = md%coord_v(m, fragment(m, i, j), 2) - hmat(2, 2)
                            CALL pbc_hexagonal(md%coord_v(m, fragment(m, i, j), :), md%coord_v(m, fragment(m, i, 1), :), sys)
                        END IF
                        IF (sys%cell%vec(1)>3.0_dp .AND. sys%cell%vec(2)<-5.2_dp) THEN
                            md%coord_v(m, fragment(m, i, j), 1) = md%coord_v(m, fragment(m, i, j), 1) - hmat(1, 2)
                            md%coord_v(m, fragment(m, i, j), 2) = md%coord_v(m, fragment(m, i, j), 2) + hmat(2, 2)
                            CALL pbc_hexagonal(md%coord_v(m, fragment(m, i, j), :), md%coord_v(m, fragment(m, i, 1), :), sys)
                        END IF
                        IF (sys%cell%vec(1)<-3.0_dp .AND. sys%cell%vec(2)<-5.2_dp) THEN
                            md%coord_v(m, fragment(m, i, j), 1) = md%coord_v(m, fragment(m, i, j), 1) + hmat(1, 2)
                            md%coord_v(m, fragment(m, i, j), 2) = md%coord_v(m, fragment(m, i, j), 2) + hmat(2, 2)
                            CALL pbc_hexagonal(md%coord_v(m, fragment(m, i, j), :), md%coord_v(m, fragment(m, i, 1), :), sys)
                        END IF
                        IF (sys%cell%vec(1)>4.9_dp .AND. sys%cell%vec(2)<5.2_dp .AND. sys%cell%vec(2)>-5.2_dp) THEN
                            md%coord_v(m, fragment(m, i, j), 1) = md%coord_v(m, fragment(m, i, j), 1) - hmat(1, 1)
                            CALL pbc_hexagonal(md%coord_v(m, fragment(m, i, j), :), md%coord_v(m, fragment(m, i, 1), :), sys)
                        END IF
                        IF (sys%cell%vec(1)<-4.9_dp .AND. sys%cell%vec(2)<5.2_dp .AND. sys%cell%vec(2)>-5.2_dp) THEN
                            md%coord_v(m, fragment(m, i, j), 1) = md%coord_v(m, fragment(m, i, j), 1) + hmat(1, 1)
                            CALL pbc_hexagonal(md%coord_v(m, fragment(m, i, j), :), md%coord_v(m, fragment(m, i, 1), :), sys)
                        END IF
                        IF (sys%cell%vec(3)<-4.9_dp) THEN
                            md%coord_v(m, fragment(m, i, j), 3) = md%coord_v(m, fragment(m, i, j), 3) + hmat(3, 3)
                            CALL pbc_hexagonal(md%coord_v(m, fragment(m, i, j), :), md%coord_v(m, fragment(m, i, 1), :), sys)
                        END IF
                        IF (sys%cell%vec(3)>5.0_dp) THEN
                            md%coord_v(m, fragment(m, i, j), 3) = md%coord_v(m, fragment(m, i, j), 3) - hmat(3, 3)
                            CALL pbc_hexagonal(md%coord_v(m, fragment(m, i, j), :), md%coord_v(m, fragment(m, i, 1), :), sys)
                        END IF
                    END DO
                END DO
            END DO

    !!!WRAPPING UP THE FRAGMENTS TOGETHER
            DO m = 1, sys%framecount
                DO i = 1, 20
                    DO j = 1, sys%fragments%natom_frag(i)
                        IF (md%coord_v(m, fragment(m, i, j), 3)>6.0_dp) THEN
                            md%coord_v(m, fragment(m, i, j), 3) = md%coord_v(m, fragment(m, i, j), 3) - hmat(3, 3)
                        END IF
                        IF (md%coord_v(m, fragment(m, i, j), 3)>0.0_dp .AND. md%coord_v(m, fragment(m, i, j), 3)<2.0_dp) THEN
                            IF (md%coord_v(m, fragment(m, i, j), 2)<-3.0_dp .AND. md%coord_v(m, fragment(m, i, j), 1)<-6.2_dp) THEN
                                md%coord_v(m, fragment(m, i, 1:sys%fragments%natom_frag(i)), 1) = md%coord_v(m, fragment(m, i, 1:sys%fragments%natom_frag(i)), 1) &
                                                                                                  + hmat(1, 2)
                                md%coord_v(m, fragment(m, i, 1:sys%fragments%natom_frag(i)), 2) = md%coord_v(m, fragment(m, i, 1:sys%fragments%natom_frag(i)), 2) &
                                                                                                  + hmat(2, 2)
                            END IF
                        END IF
                        IF (md%coord_v(m, fragment(m, i, j), 1)<-8.0_dp) THEN
                            md%coord_v(m, fragment(m, i, 1:sys%fragments%natom_frag(i)), 1) = md%coord_v(m, fragment(m, i, 1:sys%fragments%natom_frag(i)), 1) &
                                                                                              + hmat(1, 1)
                        END IF
                    END DO
                END DO
            END DO

        ELSEIF (sys%frag_type=='3') THEN
            DO m = 1, sys%framecount
                DO i = 21, sys%mol_num
                    DO j = 2, sys%fragments%natom_frag(i)
                        CALL pbc_hexagonal(md%coord_v(m, fragment(m, i, j), :), md%coord_v(m, fragment(m, i, 1), :), sys)
                        IF (sys%cell%vec(1)>3.0_dp .AND. sys%cell%vec(2)>5.2_dp) THEN
                            md%coord_v(m, fragment(m, i, j), 1) = md%coord_v(m, fragment(m, i, j), 1) - hmat(1, 2)
                            md%coord_v(m, fragment(m, i, j), 2) = md%coord_v(m, fragment(m, i, j), 2) - hmat(2, 2)
                            CALL pbc_hexagonal(md%coord_v(m, fragment(m, i, j), :), md%coord_v(m, fragment(m, i, 1), :), sys)
                        END IF
                        IF (sys%cell%vec(1)<-3.0_dp .AND. sys%cell%vec(2)>5.2_dp) THEN
                            md%coord_v(m, fragment(m, i, j), 1) = md%coord_v(m, fragment(m, i, j), 1) + hmat(1, 2)
                            md%coord_v(m, fragment(m, i, j), 2) = md%coord_v(m, fragment(m, i, j), 2) - hmat(2, 2)
                            CALL pbc_hexagonal(md%coord_v(m, fragment(m, i, j), :), md%coord_v(m, fragment(m, i, 1), :), sys)
                        END IF
                        IF (sys%cell%vec(1)>3.0_dp .AND. sys%cell%vec(2)<-5.2_dp) THEN
                            md%coord_v(m, fragment(m, i, j), 1) = md%coord_v(m, fragment(m, i, j), 1) - hmat(1, 2)
                            md%coord_v(m, fragment(m, i, j), 2) = md%coord_v(m, fragment(m, i, j), 2) + hmat(2, 2)
                            CALL pbc_hexagonal(md%coord_v(m, fragment(m, i, j), :), md%coord_v(m, fragment(m, i, 1), :), sys)
                        END IF
                        IF (sys%cell%vec(1)<-3.0_dp .AND. sys%cell%vec(2)<-5.2_dp) THEN
                            md%coord_v(m, fragment(m, i, j), 1) = md%coord_v(m, fragment(m, i, j), 1) + hmat(1, 2)
                            md%coord_v(m, fragment(m, i, j), 2) = md%coord_v(m, fragment(m, i, j), 2) + hmat(2, 2)
                            CALL pbc_hexagonal(md%coord_v(m, fragment(m, i, j), :), md%coord_v(m, fragment(m, i, 1), :), sys)
                        END IF
                        IF (sys%cell%vec(1)>4.9_dp .AND. sys%cell%vec(2)<5.2_dp .AND. sys%cell%vec(2)>-5.2_dp) THEN
                            md%coord_v(m, fragment(m, i, j), 1) = md%coord_v(m, fragment(m, i, j), 1) - hmat(1, 1)
                            CALL pbc_hexagonal(md%coord_v(m, fragment(m, i, j), :), md%coord_v(m, fragment(m, i, 1), :), sys)
                        END IF
                        IF (sys%cell%vec(1)<-4.9_dp .AND. sys%cell%vec(2)<5.2_dp .AND. sys%cell%vec(2)>-5.2_dp) THEN
                            md%coord_v(m, fragment(m, i, j), 1) = md%coord_v(m, fragment(m, i, j), 1) + hmat(1, 1)
                            CALL pbc_hexagonal(md%coord_v(m, fragment(m, i, j), :), md%coord_v(m, fragment(m, i, 1), :), sys)
                        END IF
                        IF (sys%cell%vec(3)<-4.9_dp) THEN
                            md%coord_v(m, fragment(m, i, j), 3) = md%coord_v(m, fragment(m, i, j), 3) + hmat(3, 3)
                            CALL pbc_hexagonal(md%coord_v(m, fragment(m, i, j), :), md%coord_v(m, fragment(m, i, 1), :), sys)
                        END IF
                        IF (sys%cell%vec(3)>5.0_dp) THEN
                            md%coord_v(m, fragment(m, i, j), 3) = md%coord_v(m, fragment(m, i, j), 3) - hmat(3, 3)
                            CALL pbc_hexagonal(md%coord_v(m, fragment(m, i, j), :), md%coord_v(m, fragment(m, i, 1), :), sys)
                        END IF
                    END DO
                END DO
            END DO

    !!!WRAPPING UP THE FRAGMENTS TOGETHER
            DO m = 1, sys%framecount
                DO i = 21, sys%mol_num
                    DO j = 1, sys%fragments%natom_frag(i)
                        IF (md%coord_v(m, fragment(m, i, j), 3)>6.0_dp) THEN
                            md%coord_v(m, fragment(m, i, j), 3) = md%coord_v(m, fragment(m, i, j), 3) - hmat(3, 3)
                        END IF
                        IF (md%coord_v(m, fragment(m, i, j), 3)>0.0_dp .AND. md%coord_v(m, fragment(m, i, j), 3)<2.0_dp) THEN
                            IF (md%coord_v(m, fragment(m, i, j), 2)<-3.0_dp .AND. md%coord_v(m, fragment(m, i, j), 1)<-6.2_dp) THEN
                                md%coord_v(m, fragment(m, i, 1:sys%fragments%natom_frag(i)), 1) = md%coord_v(m, fragment(m, i, 1:sys%fragments%natom_frag(i)), 1) &
                                                                                                  + hmat(1, 2)
                                md%coord_v(m, fragment(m, i, 1:sys%fragments%natom_frag(i)), 2) = md%coord_v(m, fragment(m, i, 1:sys%fragments%natom_frag(i)), 2) &
                                                                                                  + hmat(2, 2)
                            END IF
                        END IF
                        IF (md%coord_v(m, fragment(m, i, j), 1)<-8.0_dp) THEN
                            md%coord_v(m, fragment(m, i, 1:sys%fragments%natom_frag(i)), 1) = md%coord_v(m, fragment(m, i, 1:sys%fragments%natom_frag(i)), 1) &
                                                                                              + hmat(1, 1)
                        END IF
                    END DO
                END DO
            END DO

        END IF

        bec_pbc = 0.0_dp
        coord_shifted = 0.
        l = 0

        IF (sys%system=='2' .AND. dips%type_dipole=='1') THEN
            sys%mol_num = 20
        END IF

        DO m = 1, sys%framecount
            DO i = 1, sys%mol_num
                DO j = 1, sys%fragments%natom_frag(i)
                    pox_x = 50.0_dp
                    pox_y = 50.0_dp
                    pox_z = 50.0_dp
                    bec(m, :) = md%coord_v(m, fragment(m, i, j), :)
                    bec_pbc(m, 1) = bec(m, 1) - pox_x*ANINT((1./pox_x)*bec(m, 1))
                    bec_pbc(m, 2) = bec(m, 2) - pox_y*ANINT((1./pox_y)*bec(m, 2))
                    bec_pbc(m, 3) = bec(m, 3) - pox_z*ANINT((1./pox_z)*bec(m, 3))
                    com(m, i, j, :) = bec_pbc(m, :)*sys%mass_atom(fragment(m, i, j))
                    IF (sys%system=='1') THEN
                        sys%fragments%refpoint(m, i, :) = sys%fragments%refpoint(m, i, :) + com(m, i, j, :)
                    ELSEIF (dips%type_dipole=='1' .AND. sys%system=='2') THEN
                        sys%fragments%refpoint(m, 1, :) = sys%fragments%refpoint(m, 1, :) + com(m, i, j, :)  !!!For COM of whole sys%system
                    END IF
                END DO
            END DO
        END DO

        sys%fragments%mass_tot_cell = 0.0_dp

    !!get total mass of the cell
        DO j = 1, sys%natom
            sys%fragments%mass_tot_cell = sys%fragments%mass_tot_cell + sys%mass_atom(j)
        END DO

        DO m = 1, sys%framecount
            IF (dips%type_dipole=='1' .AND. sys%system=='2') THEN
                sys%fragments%refpoint(m, 1, :) = sys%fragments%refpoint(m, 1, :)/sys%fragments%mass_tot_cell
            ELSEIF (sys%system=='1') THEN
                DO i = 1, sys%mol_num
                    sys%fragments%refpoint(m, i, :) = sys%fragments%refpoint(m, i, :)/sys%fragments%mass_tot_frag(m, i)
                END DO
            END IF
        END DO

        OPEN (UNIT=60, FILE=TRIM(filename)//'-dipole_result_wholesystem.xyz', STATUS='unknown', IOSTAT=stat)
        DO m = 1, sys%framecount
            WRITE (60, *) sys%natom + 1
            WRITE (60, *)
            DO i = 1, 20
                DO j = 1, sys%fragments%natom_frag(i)
                    WRITE (60, *) sys%element(fragment(m, i, j)), md%coord_v(m, fragment(m, i, j), :)
                END DO
                !WRITE(60,*) 'N', sys%fragments%refpoint(m,i,:)
            END DO
            WRITE (60, *) 'N', sys%fragments%refpoint(m, 1, :)
        END DO
        CLOSE (60)

        OPEN (UNIT=12, FILE=TRIM(filename)//'-dipole_result_Ph.xyz', STATUS='unknown', IOSTAT=stat)
        DO m = 1, sys%framecount
            WRITE (12, *) 288
            WRITE (12, *)
            DO i = 9, 20
                DO j = 1, sys%fragments%natom_frag(i)
                    WRITE (12, *) sys%element(fragment(m, i, j)), md%coord_v(m, fragment(m, i, j), :)
                END DO
                WRITE (12, *) 'N', sys%fragments%refpoint(m, i, :)
            END DO
        END DO
        CLOSE (12)

        DEALLOCATE (com)
        PRINT *, sys%mass_atom(fragment(1, 9, 1)), sys%mass_atom(fragment(1, 20, 2)), sys%element(fragment(1, 9, 1)), &
            sys%element(fragment(1, 20, 2))
        PRINT *, sys%mass_atom(10), sys%element(10), 'mass atoms'
    END SUBROUTINE center_mass

!***************************************************************************************************
!***************************************************************************************************

    SUBROUTINE solv_frag_index(filename, natom_frag, fragment, sys, md, dips)

        TYPE(systems), INTENT(INOUT)                :: sys
        TYPE(molecular_dynamics), INTENT(INOUT)     :: md
        TYPE(dipoles), INTENT(INOUT)     :: dips
        CHARACTER(LEN=40), INTENT(INOUT)                             :: filename
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT)              :: natom_frag
        INTEGER, DIMENSION(:, :, :), ALLOCATABLE, INTENT(INOUT)          :: fragment

        CHARACTER(LEN=40)                                           :: length
        INTEGER                                                     :: stat   ! error status of OPEN statements
        INTEGER                                                     :: r, m, p, q, i, j, k, n, l, o, s, t
        REAL(kind=dp)                                                :: bec(5000, 3), bec_pbc(5000, 3)
        REAL(kind=dp)                                                :: hmat(3, 3), sqrt3, acosa, asina, a
        REAL(kind=dp)                                                :: pox_x, pox_y, pox_z
        REAL(kind=dp), DIMENSION(3)                                   :: coord3
        REAL(kind=dp), DIMENSION(:, :, :, :), ALLOCATABLE                 :: com

        !ALLOCATE(fragment(sys%framecount,37,32))
        ALLOCATE (sys%fragments%refpoint(sys%framecount, 37, 3))
        !ALLOCATE(natom_frag(37))
        ALLOCATE (sys%fragments%mass_tot_frag(sys%framecount, 37))
        ALLOCATE (com(sys%framecount, 37, 32, 3))
        !ALLOCATE(natom_frag(20))

        coord3(1) = 0.00_dp
        coord3(2) = 0.00_dp
        coord3(3) = 0.00_dp
        fragment = 0
        natom_frag = 0
        sys%fragments%refpoint = 0.00_dp
        com = 0.0_dp
        sys%fragments%mass_tot_frag = 0.0_dp

        sqrt3 = 1.73205080756887729352744634_dp

        a = 0.5_dp*(sys%cell%box_x + sys%cell%box_y)
        acosa = 0.5_dp*a
        asina = sqrt3*acosa
        hmat(1, 1) = a; hmat(1, 2) = acosa; hmat(1, 3) = 0.0_dp
        hmat(2, 1) = 0.0_dp; hmat(2, 2) = asina; hmat(2, 3) = 0.0_dp
        hmat(3, 1) = 0.0_dp; hmat(3, 2) = 0.0_dp; hmat(3, 3) = sys%cell%box_z

        PRINT *, sys%natom, sys%framecount

        !Assign fragments for solvent molecules!
        DO m = 1, sys%framecount
            j = 0
            !DO i=1,sys%natom
            DO i = 1, 98
                IF (sys%element(i).NE.'O') CYCLE
                IF (ANY(i==fragment(m, :, :))) CYCLE
                j = j + 1
                n = 0
                fragment(m, j, 1) = i
                DO k = 1, 98
                    IF (i==k) CYCLE
                    IF (ANY(k==fragment(m, :, :))) CYCLE
                    IF (sys%element(k).NE.'C') CYCLE
                    CALL pbc_hexagonal(md%coord_v(m, i, :), md%coord_v(m, k, :), sys)
                    IF (SQRT(DOT_PRODUCT(sys%cell%vec_pbc, sys%cell%vec_pbc))<1.781835_dp) THEN
                        IF (ANY(k==fragment(m, :, :))) CYCLE
                        n = n + 1
                        fragment(m, j, n + 1) = k
                        DO l = 1, 98
                            IF (k==l) CYCLE
                            IF (i==l) CYCLE
                            IF (sys%element(l)=='H') THEN
                                CALL pbc_hexagonal(md%coord_v(m, k, :), md%coord_v(m, l, :), sys)
                                IF (SQRT(DOT_PRODUCT(sys%cell%vec_pbc, sys%cell%vec_pbc))<1.451835_dp) THEN
                                    IF (ANY(l==fragment(m, :, :))) CYCLE
                                    n = n + 1
                                    fragment(m, j, n + 1) = l
                                END IF
                            END IF
                            IF (sys%element(l)=='C') THEN
                                CALL pbc_hexagonal(md%coord_v(m, k, :), md%coord_v(m, l, :), sys)
                                IF (SQRT(DOT_PRODUCT(sys%cell%vec_pbc, sys%cell%vec_pbc))<1.781835_dp) THEN
                                    IF (ANY(l==fragment(m, :, :))) CYCLE
                                    n = n + 1
                                    fragment(m, j, n + 1) = l

                                    DO o = 1, 98
                                        IF (l==o) CYCLE
                                        IF (i==o) CYCLE
                                        IF (k==o) CYCLE
                                        IF (sys%element(o)=='H') THEN
                                            CALL pbc_hexagonal(md%coord_v(m, l, :), md%coord_v(m, o, :), sys)
                                            IF (SQRT(DOT_PRODUCT(sys%cell%vec_pbc, sys%cell%vec_pbc))<1.451835_dp) THEN
                                                IF (ANY(o==fragment(m, :, :))) CYCLE
                                                n = n + 1
                                                fragment(m, j, n + 1) = o
                                            END IF
                                        END IF
                                        IF (sys%element(o)=='O') THEN
                                            CALL pbc_hexagonal(md%coord_v(m, l, :), md%coord_v(m, o, :), sys)
                                            IF (SQRT(DOT_PRODUCT(sys%cell%vec_pbc, sys%cell%vec_pbc))<1.781835_dp) THEN
                                                IF (ANY(o==fragment(m, :, :))) CYCLE
                                                n = n + 1
                                                fragment(m, j, n + 1) = o
                                            END IF
                                        END IF
                                    END DO
                                END IF
                            END IF
                        END DO
                    END IF
                END DO
            END DO
        END DO

        !Assign wannier centers to solvent fragments
        DO m = 1, sys%framecount
            !print*,m
            DO i = 1, 7
                p = 0
                DO j = 1, 14
                    IF (sys%element(fragment(m, i, j))=='H') CYCLE
                    DO k = 1, sys%natom
                        IF (sys%element(k).NE.'X') CYCLE
                        IF (sys%element(fragment(m, i, j))=='O' .OR. sys%element(fragment(m, i, j))=='C') THEN
                            CALL pbc_hexagonal(md%coord_v(m, fragment(m, i, j), :), md%coord_v(m, k, :), sys)
                            IF (SQRT(DOT_PRODUCT(sys%cell%vec_pbc, sys%cell%vec_pbc))<0.99_dp) THEN
                                IF (ANY(k==fragment(m, :, :))) CYCLE
                                p = p + 1
                                fragment(m, i, 14 + p) = k
                            END IF
                        END IF
                    END DO
                END DO
            END DO
        END DO

        p = 0

        !Assign fragments for the B-O ring!
        n = 0
        DO m = 1, sys%framecount
            j = 7
            !DO i=1,sys%natom
            DO i = 99, sys%natom
                IF (sys%element(i).NE.'B') CYCLE
                IF (ANY(i==fragment(m, :, :))) CYCLE
                j = j + 1
                n = 0
                fragment(m, j, 1) = i
                DO k = 99, sys%natom
                    IF (i==k) CYCLE
                    IF (ANY(k==fragment(m, :, :))) CYCLE
                    IF (sys%element(k).NE.'O') CYCLE
                    CALL pbc_hexagonal(md%coord_v(m, i, :), md%coord_v(m, k, :), sys)
                    IF (SQRT(DOT_PRODUCT(sys%cell%vec_pbc, sys%cell%vec_pbc))<1.551835_dp) THEN
                        IF (ANY(k==fragment(m, :, :))) CYCLE
                        n = n + 1
                        fragment(m, j, n + 1) = k
                        DO l = 99, sys%natom
                            IF (k==l) CYCLE
                            IF (i==l) CYCLE
                            IF (sys%element(l).NE.'B') CYCLE
                            CALL pbc_hexagonal(md%coord_v(m, k, :), md%coord_v(m, l, :), sys)
                            IF (SQRT(DOT_PRODUCT(sys%cell%vec_pbc, sys%cell%vec_pbc))<1.551835_dp) THEN
                                IF (ANY(l==fragment(m, :, :))) CYCLE
                                n = n + 1
                                fragment(m, j, n + 1) = l
                                DO o = 99, sys%natom
                                    IF (l==o) CYCLE
                                    IF (i==o) CYCLE
                                    IF (k==o) CYCLE
                                    IF (sys%element(o).NE.'O') CYCLE
                                    CALL pbc_hexagonal(md%coord_v(m, l, :), md%coord_v(m, o, :), sys)
                                    IF (SQRT(DOT_PRODUCT(sys%cell%vec_pbc, sys%cell%vec_pbc))<1.551835_dp) THEN
                                        IF (ANY(o==fragment(m, :, :))) CYCLE
                                        n = n + 1
                                        fragment(m, j, n + 1) = o
                                    END IF
                                END DO
                            END IF
                        END DO
                    END IF
                END DO
            END DO
        END DO

        !Assign Wannier centers to the B-O rings!
        DO m = 1, sys%framecount
            DO i = 8, 19
                n = 0
                DO j = 1, 6
                    IF (sys%element(fragment(m, i, j)).NE.'O') CYCLE
                    DO k = 1, sys%natom
                        IF (j==k) CYCLE
                        IF (sys%element(k).NE.'X') CYCLE
                        CALL pbc_hexagonal(md%coord_v(m, fragment(m, i, j), :), md%coord_v(m, k, :), sys)
                        IF (SQRT(DOT_PRODUCT(sys%cell%vec_pbc, sys%cell%vec_pbc))<0.7_dp) THEN
                            IF (ANY(k==fragment(m, :, :))) CYCLE
                            n = n + 1
                            fragment(m, i, n + 6) = k
                        END IF
                    END DO
                END DO
            END DO
        END DO

        PRINT *, fragment(5000, 19, :)

        !Assign fragments for the phenyl rings!
        n = 0
        DO m = 1, sys%framecount
            j = 19
            DO i = 99, sys%natom
                IF (sys%element(i).NE.'H') CYCLE
                IF (ANY(i==fragment(m, :, :))) CYCLE
                j = j + 1
                !        print*,j,'j',sys%element(i),i
                n = 0
                fragment(m, j, 1) = i
                DO k = 99, sys%natom
                    IF (i==k) CYCLE
                    ! IF (ANY(k==fragment(m,:,:))) CYCLE
                    IF (sys%element(k).NE.'C') CYCLE
                    CALL pbc_hexagonal(md%coord_v(m, i, :), md%coord_v(m, k, :), sys)
                    IF (SQRT(DOT_PRODUCT(sys%cell%vec_pbc, sys%cell%vec_pbc))<1.31835_dp) THEN
                        IF (ANY(k==fragment(m, :, :))) CYCLE
                        n = n + 1
                        fragment(m, j, n + 1) = k
                        !        print*,k,sys%element(k),i,sys%element(i), 'k ve i'
                        DO l = 99, sys%natom
                            !IF (ANY(l==fragment(m,:,:))) CYCLE
                            IF (sys%element(l).NE.'C') CYCLE
                            CALL pbc_hexagonal(md%coord_v(m, k, :), md%coord_v(m, l, :), sys)
                            IF (SQRT(DOT_PRODUCT(sys%cell%vec_pbc, sys%cell%vec_pbc))<1.55835_dp) THEN
                                IF (ANY(l==fragment(m, :, :))) CYCLE
                                n = n + 1
                                fragment(m, j, n + 1) = l
                                !        print*,l,sys%element(l),k,sys%element(k),'l ve k'
                                DO o = 99, sys%natom
                                    ! IF (ANY(o==fragment(m,:,:))) CYCLE
                                    IF (sys%element(o)=='X') CYCLE
                                    CALL pbc_hexagonal(md%coord_v(m, l, :), md%coord_v(m, o, :), sys)
                                    IF (SQRT(DOT_PRODUCT(sys%cell%vec_pbc, sys%cell%vec_pbc))<1.55835_dp) THEN
                                        IF (ANY(o==fragment(m, :, :))) CYCLE
                                        n = n + 1
                                        fragment(m, j, n + 1) = o
                                        ! print*,o,sys%element(o),l,sys%element(l),'o ve l'
                                        DO p = 99, sys%natom
                                            ! IF (ANY(p==fragment(m,:,:))) CYCLE
                                            IF (sys%element(p)=='X') CYCLE
                                            CALL pbc_hexagonal(md%coord_v(m, o, :), md%coord_v(m, p, :), sys)
                                            IF (SQRT(DOT_PRODUCT(sys%cell%vec_pbc, sys%cell%vec_pbc))<1.55835_dp) THEN
                                                IF (ANY(p==fragment(m, :, :))) CYCLE
                                                n = n + 1
                                                fragment(m, j, n + 1) = p
                                                !        print*,p,sys%element(p),o,sys%element(o),'p ve o'
                                                DO r = 99, sys%natom
                                                    !IF (ANY(r==fragment(m,:,:))) CYCLE
                                                    IF (sys%element(r)=='X') CYCLE
                                                    CALL pbc_hexagonal(md%coord_v(m, p, :), md%coord_v(m, r, :), sys)
                                                    IF (SQRT(DOT_PRODUCT(sys%cell%vec_pbc, sys%cell%vec_pbc))<1.55835_dp) THEN
                                                        IF (ANY(r==fragment(m, :, :))) CYCLE
                                                        n = n + 1
                                                        fragment(m, j, n + 1) = r
                                                        !        print*,r,sys%element(r),p,sys%element(p),'r ve p'
                                                        DO s = 99, sys%natom
                                                            !IF (ANY(s==fragment(m,:,:))) CYCLE
                                                            IF (sys%element(s)=='X') CYCLE
                                                            CALL pbc_hexagonal(md%coord_v(m, r, :), md%coord_v(m, s, :), sys)
                                                            IF (SQRT(DOT_PRODUCT(sys%cell%vec_pbc, sys%cell%vec_pbc))<1.55835_dp) THEN
                                                                IF (ANY(s==fragment(m, :, :))) CYCLE
                                                                n = n + 1
                                                                fragment(m, j, n + 1) = s
                                                                !        print*,s,sys%element(s),r,sys%element(r),'s ve r'
                                                                DO t = 99, sys%natom
                                                                    !IF (ANY(s==fragment(m,:,:))) CYCLE
                                                                    IF (sys%element(t)=='X') CYCLE
                                                                    CALL pbc_hexagonal(md%coord_v(m, s, :), md%coord_v(m, t, :), sys)
                                                                    IF (SQRT(DOT_PRODUCT(sys%cell%vec_pbc, sys%cell%vec_pbc))<1.55835_dp) THEN
                                                                        IF (ANY(t==fragment(m, :, :))) CYCLE
                                                                        n = n + 1
                                                                        fragment(m, j, n + 1) = t
                                                                        !         print*,t,sys%element(t),s,sys%element(s),'t ve s'
                                                                    END IF
                                                                END DO
                                                            END IF
                                                        END DO
                                                    END IF
                                                END DO
                                            END IF
                                        END DO
                                    END IF
                                END DO
                            END IF
                        END DO
                    END IF
                END DO
            END DO
        END DO

        PRINT *, fragment(1, 20, :)

        !Assign Wannier centers to the phenyl rings!

        DO m = 1, sys%framecount
            DO i = 20, 37
                n = 0
                DO j = 1, 10
                    IF (sys%element(fragment(m, i, j)).NE.'C') CYCLE
                    DO k = 1, sys%natom
                        IF (j==k) CYCLE
                        IF (sys%element(k).NE.'X') CYCLE
                        CALL pbc_hexagonal(md%coord_v(m, fragment(m, i, j), :), md%coord_v(m, k, :), sys)
                        IF (SQRT(DOT_PRODUCT(sys%cell%vec_pbc, sys%cell%vec_pbc))<1.0_dp) THEN
                            IF (ANY(k==fragment(m, i, :))) CYCLE
                            n = n + 1
                            fragment(m, i, n + 10) = k
                        END IF
                    END DO
                END DO
            END DO
        END DO

        !!get natom_frag
        DO i = 1, sys%framecount
            DO j = 1, 37
                natom_frag(j) = COUNT(fragment(i, j, :).NE.0)
            END DO
        END DO

        !!get total mass
        DO m = 1, sys%framecount
            DO i = 1, 37
                DO j = 1, natom_frag(i)
                    IF (sys%element(natom_frag(i))=='X') CYCLE
                    sys%fragments%mass_tot_frag(m, i) = sys%fragments%mass_tot_frag(m, i) + sys%mass_atom(fragment(m, i, j))
                END DO
            END DO
        END DO

        !!Find center of mass

        DO m = 1, sys%framecount
            DO i = 1, 37
                DO j = 2, natom_frag(i)
                    CALL pbc_hexagonal(md%coord_v(m, fragment(m, i, j), :), md%coord_v(m, fragment(m, i, 1), :), sys)
                    ! PRINT*,SQRT(DOT_PRODUCT(vec,vec)),SQRT(DOT_PRODUCT(sys%cell%vec_pbc,sys%cell%vec_pbc)),i,j,'0'
                    IF (sys%cell%vec(1)>3.0_dp .AND. sys%cell%vec(2)>5.2_dp) THEN
                        md%coord_v(m, fragment(m, i, j), 1) = md%coord_v(m, fragment(m, i, j), 1) - hmat(1, 2)
                        md%coord_v(m, fragment(m, i, j), 2) = md%coord_v(m, fragment(m, i, j), 2) - hmat(2, 2)
                        CALL pbc_hexagonal(md%coord_v(m, fragment(m, i, j), :), md%coord_v(m, fragment(m, i, 1), :), sys)
                        ! PRINT*,SQRT(DOT_PRODUCT(vec,vec)),SQRT(DOT_PRODUCT(sys%cell%vec_pbc,sys%cell%vec_pbc)),i,j,'3'
                    END IF
                    IF (sys%cell%vec(1)<-3.0_dp .AND. sys%cell%vec(2)>5.2_dp) THEN
                        md%coord_v(m, fragment(m, i, j), 1) = md%coord_v(m, fragment(m, i, j), 1) + hmat(1, 2)
                        md%coord_v(m, fragment(m, i, j), 2) = md%coord_v(m, fragment(m, i, j), 2) - hmat(2, 2)
                        CALL pbc_hexagonal(md%coord_v(m, fragment(m, i, j), :), md%coord_v(m, fragment(m, i, 1), :), sys)
                        ! PRINT*,SQRT(DOT_PRODUCT(vec,vec)),SQRT(DOT_PRODUCT(sys%cell%vec_pbc,sys%cell%vec_pbc)),i,j,'4'
                    END IF
                    IF (sys%cell%vec(1)>3.0_dp .AND. sys%cell%vec(2)<-5.2_dp) THEN
                        md%coord_v(m, fragment(m, i, j), 1) = md%coord_v(m, fragment(m, i, j), 1) - hmat(1, 2)
                        md%coord_v(m, fragment(m, i, j), 2) = md%coord_v(m, fragment(m, i, j), 2) + hmat(2, 2)
                        CALL pbc_hexagonal(md%coord_v(m, fragment(m, i, j), :), md%coord_v(m, fragment(m, i, 1), :), sys)
                        ! PRINT*,SQRT(DOT_PRODUCT(vec,vec)),SQRT(DOT_PRODUCT(sys%cell%vec_pbc,sys%cell%vec_pbc)),i,j,'5'
                    END IF
                    IF (sys%cell%vec(1)<-3.0_dp .AND. sys%cell%vec(2)<-5.2_dp) THEN
                        md%coord_v(m, fragment(m, i, j), 1) = md%coord_v(m, fragment(m, i, j), 1) + hmat(1, 2)
                        md%coord_v(m, fragment(m, i, j), 2) = md%coord_v(m, fragment(m, i, j), 2) + hmat(2, 2)
                        CALL pbc_hexagonal(md%coord_v(m, fragment(m, i, j), :), md%coord_v(m, fragment(m, i, 1), :), sys)
                        ! PRINT*,SQRT(DOT_PRODUCT(vec,vec)),SQRT(DOT_PRODUCT(sys%cell%vec_pbc,sys%cell%vec_pbc)),i,j,'6'
                    END IF
                    IF (sys%cell%vec(1)>4.9_dp .AND. sys%cell%vec(2)<5.2_dp .AND. sys%cell%vec(2)>-5.2_dp) THEN
                        md%coord_v(m, fragment(m, i, j), 1) = md%coord_v(m, fragment(m, i, j), 1) - hmat(1, 1)
                        CALL pbc_hexagonal(md%coord_v(m, fragment(m, i, j), :), md%coord_v(m, fragment(m, i, 1), :), sys)
                        ! PRINT*,SQRT(DOT_PRODUCT(vec,vec)),SQRT(DOT_PRODUCT(sys%cell%vec_pbc,sys%cell%vec_pbc)),i,j,'1'
                    END IF
                    IF (sys%cell%vec(1)<-4.9_dp .AND. sys%cell%vec(2)<5.2_dp .AND. sys%cell%vec(2)>-5.2_dp) THEN
                        md%coord_v(m, fragment(m, i, j), 1) = md%coord_v(m, fragment(m, i, j), 1) + hmat(1, 1)
                        CALL pbc_hexagonal(md%coord_v(m, fragment(m, i, j), :), md%coord_v(m, fragment(m, i, 1), :), sys)
                        ! PRINT*,SQRT(DOT_PRODUCT(vec,vec)),SQRT(DOT_PRODUCT(sys%cell%vec_pbc,sys%cell%vec_pbc)),i,j,'2'
                    END IF
                    IF (sys%cell%vec(3)<-4.9_dp) THEN
                        md%coord_v(m, fragment(m, i, j), 3) = md%coord_v(m, fragment(m, i, j), 3) + hmat(3, 3)
                        CALL pbc_hexagonal(md%coord_v(m, fragment(m, i, j), :), md%coord_v(m, fragment(m, i, 1), :), sys)
                        ! PRINT*,SQRT(DOT_PRODUCT(vec,vec)),SQRT(DOT_PRODUCT(sys%cell%vec_pbc,sys%cell%vec_pbc)),i,j,'7'
                    END IF
                    IF (sys%cell%vec(3)>5.0_dp) THEN
                        md%coord_v(m, fragment(m, i, j), 3) = md%coord_v(m, fragment(m, i, j), 3) - hmat(3, 3)
                        CALL pbc_hexagonal(md%coord_v(m, fragment(m, i, j), :), md%coord_v(m, fragment(m, i, 1), :), sys)
                        ! PRINT*,SQRT(DOT_PRODUCT(vec,vec)),SQRT(DOT_PRODUCT(sys%cell%vec_pbc,sys%cell%vec_pbc)),i,j,'8'
                    END IF
                END DO
            END DO
        END DO

        bec_pbc = 0.0_dp

        DO i = 1, 37
            DO j = 1, natom_frag(i)
                DO m = 1, sys%framecount
                    pox_x = 40.0_dp
                    pox_y = 40.0_dp
                    pox_z = 40.0_dp
                    bec(m, :) = md%coord_v(m, fragment(m, i, j), :)
                    bec_pbc(m, 1) = bec(m, 1) - pox_x*ANINT((1./pox_x)*bec(m, 1))
                    bec_pbc(m, 2) = bec(m, 2) - pox_y*ANINT((1./pox_y)*bec(m, 2))
                    bec_pbc(m, 3) = bec(m, 3) - pox_z*ANINT((1./pox_z)*bec(m, 3))
                    com(m, i, j, :) = bec_pbc(m, :)*sys%mass_atom(fragment(m, i, j))
                    sys%fragments%refpoint(m, i, :) = sys%fragments%refpoint(m, i, :) + com(m, i, j, :)
                END DO
            END DO
        END DO

        DO m = 1, sys%framecount
            DO i = 1, 37
                sys%fragments%refpoint(m, i, :) = sys%fragments%refpoint(m, i, :)/sys%fragments%mass_tot_frag(m, i)
            END DO
        END DO

        OPEN (UNIT=60, FILE=TRIM(filename)//'-dipole_result_7.xyz', STATUS='unknown', IOSTAT=stat)
        DO m = 1, sys%framecount
            WRITE (60, *) sys%natom + 37
            WRITE (60, *)
            DO i = 1, 37
                DO j = 1, natom_frag(i)
                    WRITE (60, *) sys%element(fragment(m, i, j)), md%coord_v(m, fragment(m, i, j), :)
                END DO
                WRITE (60, *) 'N', sys%fragments%refpoint(m, i, :)
            END DO
        END DO
        CLOSE (60)

        DEALLOCATE (com)

    END SUBROUTINE solv_frag_index

!!***************************************************************************************************
!!***************************************************************************************************
    SUBROUTINE wannier_frag(natom_frag, filename, dipole, fragment, gs, sys, md, dips)

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
        DO m = 1, sys%framecount
            DO i = 1, sys%mol_num
                DO j = 1, natom_frag(i)
                    IF (sys%system=='1') THEN
                        CALL pbc_hexagonal(md%coord_v(m, fragment(m, i, j), :), sys%fragments%refpoint(m, i, :), sys) !! If fragment approach!!!
                        dipole(m, i, :) = dipole(m, i, :) + sys%cell%vec_pbc*1.889725989_dp*sys%charge(fragment(m, i, j))
                    ELSEIF (dips%type_dipole=='1' .AND. sys%system=='2') THEN
                        CALL pbc_orthorombic_old(md%coord_v(m, fragment(m, i, j), :), sys%fragments%refpoint(m, 1, :), sys%cell%vec, sys%cell%vec_pbc, &
                                                 pox_all, pox_x, pox_y, pox_z) !! IF COM is for whole supercel!!!
                        dipole(m, 1, :) = dipole(m, 1, :) + sys%cell%vec_pbc*1.889725989_dp*sys%charge(fragment(m, i, j))
                    END IF
                END DO
            END DO
        END DO

        PRINT *, sys%fragments%mass_tot_cell, 'mass tot cell'

        DO m = 1, sys%framecount
            IF (sys%system=='2' .AND. dips%type_dipole=='1') THEN
                IF (dipole(m, 1, 1)>120) THEN
                    dipole(m, 1, 1) = dipole(m, 1, 1) - (hmat(1, 1)*3*1.889725989)
                END IF
!  dipole(m,1,1)=dipole(m,1,1)-(hmat(1,1)*3*1.889725989)
                dipole(m, 1, 2) = dipole(m, 1, 2) + (hmat(2, 2)*4*1.889725989)!+(hmat(1,2)*3*1.889725989)
                dipole(m, 1, :) = REAL((dipole(m, 1, :)*2.54174622741_dp)/sys%fragments%mass_tot_cell, kind=dp) !!Dividing by total mass, change it later
            ELSEIF (sys%system=='1') THEN
                DO i = 1, sys%mol_num !! define something for this
                    dipole(m, i, :) = REAL((dipole(m, i, :)*2.54174622741_dp)/sys%fragments%mass_tot_frag(m, i), kind=dp) !!Dividing by total mass, change it later
                END DO
            END IF
        END DO

!DO j=21,44 !! define something for this
! WRITE(filename_dip, '(i0)') j
! OPEN(UNIT=68,FILE='COF-1_dipoles-'//trim(filename_dip)//'.txt',STATUS='unknown',IOSTAT=stat)
! DO m=1,sys%framecount
!  WRITE(68,*) m,dipole(m,j,1:3)!,vec_pbc(3)
! ENDDO
!ENDDO
!CLOSE(68)

!OPEN(UNIT=68,FILE='dipole_result_final.xyz',STATUS='unknown',IOSTAT=stat)
!DO m=1,sys%framecount
!DO j=1,20
        !  WRITE(68,'(2X,A15,4X,I2,6X,F20.12)') "net dipole",j, SQRT(DOT_PRODUCT(dipole(m,j,:),dipole(m,j,:)))*2.54174622741_dp

! ENDDO
! ENDDO
!CLOSE(68)

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
    END SUBROUTINE wannier_frag
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

