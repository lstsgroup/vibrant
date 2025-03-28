MODULE dipole_calc

    USE kinds, ONLY: dp
    USE setup, ONLY: pbc_orthorombic, pbc_hexagonal

    IMPLICIT NONE

    PUBLIC :: wannier, center_mass, wannier_frag, solv_frag_index
    PRIVATE

CONTAINS

    SUBROUTINE center_mass(natom_frag, natom, refpoint, coord_v, filename, element, box_all, box_x, box_y, &
                           box_z, vec, vec_pbc, fragment, mass_atom, framecount, cell_type, mass_tot_frag, frag_type, mol_num, &
                           nfrag, type_dipole, system, mass_tot_cell)

        CHARACTER(LEN=40), INTENT(INOUT)                             :: filename, cell_type, type_dipole, system
        CHARACTER(LEN=2), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)     :: element
        CHARACTER(LEN=40), INTENT(IN)                                :: frag_type
        INTEGER, INTENT(INOUT)                                       :: natom, framecount, mol_num
        INTEGER, INTENT(OUT)                                         :: nfrag
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT)                :: natom_frag
        INTEGER, DIMENSION(:, :, :), ALLOCATABLE, INTENT(OUT)            :: fragment
        REAL(kind=dp), INTENT(INOUT)                                  :: box_all, box_x, box_y, box_z, vec(3), vec_pbc(3)
        REAL(kind=dp), INTENT(OUT)                                    :: mass_tot_cell
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)         :: mass_atom
        REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE, INTENT(OUT)         :: mass_tot_frag
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE, INTENT(OUT)       :: refpoint
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE, INTENT(INOUT)     :: coord_v

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

        mol_num = 44
        ALLOCATE (com(framecount, mol_num, 35, 3))
        ALLOCATE (coord_shifted(framecount, natom, 3))
        ALLOCATE (com2(framecount, natom, 3))
        ALLOCATE (mass_tot_frag(framecount, mol_num))
        ALLOCATE (refpoint(framecount, mol_num, 3))
        ALLOCATE (fragment(framecount, mol_num, 32))
        ALLOCATE (natom_frag(mol_num))
        ALLOCATE (element_shifted(framecount, natom))

        coord3(1) = 0.00_dp
        coord3(2) = 0.00_dp
        coord3(3) = 0.00_dp
        nfrag = 0
        fragment = 0
        natom_frag = 0
        refpoint = 0.00_dp
        com = 0.0_dp
        com2 = 0.0_dp
        mass_tot_frag = 0.0_dp
        n = 0
        l = 0
        q = 0
        sqrt3 = 1.73205080756887729352744634_dp

        a = 0.5_dp*(box_x + box_y)
        acosa = 0.5_dp*a
        asina = sqrt3*acosa
        hmat(1, 1) = a; hmat(1, 2) = acosa; hmat(1, 3) = 0.0_dp
        hmat(2, 1) = 0.0_dp; hmat(2, 2) = asina; hmat(2, 3) = 0.0_dp
        hmat(3, 1) = 0.0_dp; hmat(3, 2) = 0.0_dp; hmat(3, 3) = box_z

!Assign fragments for the B-O ring!
        n = 0
        p = 0
        DO m = 1, framecount
            j = 0
            DO i = 1, natom
                IF (element(i).NE.'B') CYCLE
                IF (ANY(i==fragment(m, :, :))) CYCLE
                j = j + 1
                n = 0
                fragment(m, j, 1) = i
                DO k = 1, natom
                    IF (i==k) CYCLE
                    IF (ANY(k==fragment(m, :, :))) CYCLE
                    IF (element(k).NE.'O') CYCLE
                    CALL pbc_hexagonal(coord_v(m, i, :), coord_v(m, k, :), &
                                       vec, vec_pbc, box_all, box_x, box_y, box_z)
                    IF (SQRT(DOT_PRODUCT(vec_pbc, vec_pbc))<1.551835_dp) THEN
                        IF (ANY(k==fragment(m, :, :))) CYCLE
                        n = n + 1
                        fragment(m, j, n + 1) = k
                        DO l = 1, natom
                            IF (k==l) CYCLE
                            IF (i==l) CYCLE
                            IF (element(l).NE.'B') CYCLE
                            CALL pbc_hexagonal(coord_v(m, k, :), coord_v(m, l, :), &
                                               vec, vec_pbc, box_all, box_x, box_y, box_z)
                            IF (SQRT(DOT_PRODUCT(vec_pbc, vec_pbc))<1.551835_dp) THEN
                                IF (ANY(l==fragment(m, :, :))) CYCLE
                                n = n + 1
                                fragment(m, j, n + 1) = l
                                DO o = 1, natom
                                    IF (l==o) CYCLE
                                    IF (i==o) CYCLE
                                    IF (k==o) CYCLE
                                    IF (element(o).NE.'O') CYCLE
                                    CALL pbc_hexagonal(coord_v(m, l, :), coord_v(m, o, :), &
                                                       vec, vec_pbc, box_all, box_x, box_y, box_z)
                                    IF (SQRT(DOT_PRODUCT(vec_pbc, vec_pbc))<1.551835_dp) THEN
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
        DO m = 1, framecount
            j = 8
            DO i = 1, natom
                IF (element(i).NE.'H') CYCLE
                IF (ANY(i==fragment(m, :, :))) CYCLE
                j = j + 1
                n = 0
                fragment(m, j, 1) = i
                DO k = 1, natom
                    IF (i==k) CYCLE
                    IF (element(k).NE.'C') CYCLE
                    CALL pbc_hexagonal(coord_v(m, i, :), coord_v(m, k, :), &
                                       vec, vec_pbc, box_all, box_x, box_y, box_z)
                    IF (SQRT(DOT_PRODUCT(vec_pbc, vec_pbc))<1.31835_dp) THEN
                        IF (ANY(k==fragment(m, :, :))) CYCLE
                        n = n + 1
                        fragment(m, j, n + 1) = k
                        DO l = 1, natom
                            IF (element(l).NE.'C') CYCLE
                            CALL pbc_hexagonal(coord_v(m, k, :), coord_v(m, l, :), &
                                               vec, vec_pbc, box_all, box_x, box_y, box_z)
                            IF (SQRT(DOT_PRODUCT(vec_pbc, vec_pbc))<1.55835_dp) THEN
                                IF (ANY(l==fragment(m, :, :))) CYCLE
                                n = n + 1
                                fragment(m, j, n + 1) = l
                                DO o = 1, natom
                                    IF (element(o)=='X') CYCLE
                                    CALL pbc_hexagonal(coord_v(m, l, :), coord_v(m, o, :), &
                                                       vec, vec_pbc, box_all, box_x, box_y, box_z)
                                    IF (SQRT(DOT_PRODUCT(vec_pbc, vec_pbc))<1.55835_dp) THEN
                                        IF (ANY(o==fragment(m, :, :))) CYCLE
                                        n = n + 1
                                        fragment(m, j, n + 1) = o
                                        DO p = 1, natom
                                            IF (element(p)=='X') CYCLE
                                            CALL pbc_hexagonal(coord_v(m, o, :), coord_v(m, p, :), &
                                                               vec, vec_pbc, box_all, box_x, box_y, box_z)
                                            IF (SQRT(DOT_PRODUCT(vec_pbc, vec_pbc))<1.55835_dp) THEN
                                                IF (ANY(p==fragment(m, :, :))) CYCLE
                                                n = n + 1
                                                fragment(m, j, n + 1) = p
                                                DO r = 1, natom
                                                    IF (element(r)=='X') CYCLE
                                                    CALL pbc_hexagonal(coord_v(m, p, :), coord_v(m, r, :), &
                                                                       vec, vec_pbc, box_all, box_x, box_y, box_z)
                                                    IF (SQRT(DOT_PRODUCT(vec_pbc, vec_pbc))<1.55835_dp) THEN
                                                        IF (ANY(r==fragment(m, :, :))) CYCLE
                                                        n = n + 1
                                                        fragment(m, j, n + 1) = r
                                                        DO s = 1, natom
                                                            IF (element(s)=='X') CYCLE
                                                            CALL pbc_hexagonal(coord_v(m, r, :), coord_v(m, s, :), &
                                                                               vec, vec_pbc, box_all, box_x, box_y, box_z)
                                                            IF (SQRT(DOT_PRODUCT(vec_pbc, vec_pbc))<1.55835_dp) THEN
                                                                IF (ANY(s==fragment(m, :, :))) CYCLE
                                                                n = n + 1
                                                                fragment(m, j, n + 1) = s
                                                                DO t = 1, natom
                                                                    IF (element(t)=='X') CYCLE
                                                                    CALL pbc_hexagonal(coord_v(m, s, :), coord_v(m, t, :), &
                                                                                       vec, vec_pbc, box_all, box_x, box_y, box_z)
                                                                    IF (SQRT(DOT_PRODUCT(vec_pbc, vec_pbc))<1.55835_dp) THEN
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
        DO m = 1, framecount
            j = 20
            DO i = 1, natom
                IF (element(i).NE.'B') CYCLE
                n = 0
                j = j + 1
                fragment(m, j, 1) = i
                DO k = 1, natom
                    IF (element(k).NE.'C') CYCLE
                    CALL pbc_hexagonal(coord_v(m, i, :), coord_v(m, k, :), vec, vec_pbc, box_all, box_x, box_y, box_z)
                    IF (SQRT(DOT_PRODUCT(vec_pbc, vec_pbc))<2.0_dp) THEN
                        n = n + 1
                        fragment(m, j, n + 1) = k
                    END IF
                END DO
            END DO
        END DO

!!! Assign wannier centers!!!

        DO m = 1, framecount
            DO i = 1, 8
                n = 0
                DO j = 1, 6
                    IF (element(fragment(m, i, j)).NE.'O') CYCLE
                    DO k = 1, natom
                        IF (j==k) CYCLE
                        IF (element(k).NE.'X') CYCLE
                        CALL pbc_hexagonal(coord_v(m, fragment(m, i, j), :), coord_v(m, k, :), vec, vec_pbc, &
                                           box_all, box_x, box_y, box_z)
                        IF (SQRT(DOT_PRODUCT(vec_pbc, vec_pbc))<0.7_dp) THEN
                            IF (ANY(k==fragment(m, :, :))) CYCLE
                            n = n + 1
                            fragment(m, i, n + 6) = k
                        END IF
                    END DO
                END DO
            END DO
        END DO

        l = 0

        DO m = 1, framecount
            DO i = 9, 20
                n = 0
                DO j = 1, 10
                    IF (element(fragment(m, i, j)).NE.'C') CYCLE
                    outer: DO k = 1, natom
                        IF (element(k).NE.'X') CYCLE
                        CALL pbc_hexagonal(coord_v(m, fragment(m, i, j), :), coord_v(m, k, :), vec, vec_pbc, &
                                           box_all, box_x, box_y, box_z)
                        IF (SQRT(DOT_PRODUCT(vec_pbc, vec_pbc))<1.0_dp) THEN
                            IF (system=='1') THEN
                                inner: DO l = 1, natom
                                    IF (element(l).NE.'B') CYCLE
                                    CALL pbc_hexagonal(coord_v(m, k, :), coord_v(m, l, :), vec, vec_pbc, &
                                                       box_all, box_x, box_y, box_z)
                                    IF (SQRT(DOT_PRODUCT(vec_pbc, vec_pbc))<1.1_dp) CYCLE outer
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

        DO m = 1, framecount
            DO i = 21, mol_num
                DO j = 1, 2
                    IF (element(fragment(m, i, j)).NE.'C') CYCLE
                    DO k = 1, natom
                        IF (element(k).NE.'X') CYCLE
                        IF (ANY(k==fragment(m, :, :))) CYCLE
                        CALL pbc_hexagonal(coord_v(m, fragment(m, i, j), :), coord_v(m, k, :), vec, vec_pbc, &
                                           box_all, box_x, box_y, box_z)
                        IF (SQRT(DOT_PRODUCT(vec_pbc, vec_pbc))<1.0_dp) THEN
                            fragment(m, i, 3) = k
                        END IF
                    END DO
                END DO
            END DO
        END DO

!!get natom_frag
        DO i = 1, framecount
            DO j = 1, mol_num
                natom_frag(j) = COUNT(fragment(i, j, :).NE.0)
            END DO
        END DO

        IF (frag_type=='1') THEN
            nfrag = 8
        ELSEIF (frag_type=='2') THEN
            nfrag = 12
        ELSEIF (frag_type=='3') THEN
            nfrag = 24
        ELSEIF (system=='2' .AND. type_dipole=='1') THEN
            nfrag = 1
        END IF

!!get total mass
        DO m = 1, framecount
            DO i = 1, mol_num
                DO j = 1, natom_frag(i)
                    mass_tot_frag(m, i) = mass_tot_frag(m, i) + mass_atom(fragment(m, i, j))
                END DO
            END DO
        END DO

!!Find center of mass

!!!For B-O and C-B fragments!!!
        IF (frag_type=='1' .OR. frag_type=='2' .OR. (type_dipole=='1' .AND. system=='2')) THEN
            DO m = 1, framecount
                DO i = 1, 20
                    DO j = 2, natom_frag(i)
                        CALL pbc_hexagonal(coord_v(m, fragment(m, i, j), :), coord_v(m, fragment(m, i, 1), :), vec, vec_pbc, &
                                           box_all, box_x, box_y, box_z)
                        IF (vec(1)>3.0_dp .AND. vec(2)>5.2_dp) THEN
                            coord_v(m, fragment(m, i, j), 1) = coord_v(m, fragment(m, i, j), 1) - hmat(1, 2)
                            coord_v(m, fragment(m, i, j), 2) = coord_v(m, fragment(m, i, j), 2) - hmat(2, 2)
                            CALL pbc_hexagonal(coord_v(m, fragment(m, i, j), :), coord_v(m, fragment(m, i, 1), :), vec, vec_pbc, &
                                               box_all, box_x, box_y, box_z)
                        END IF
                        IF (vec(1)<-3.0_dp .AND. vec(2)>5.2_dp) THEN
                            coord_v(m, fragment(m, i, j), 1) = coord_v(m, fragment(m, i, j), 1) + hmat(1, 2)
                            coord_v(m, fragment(m, i, j), 2) = coord_v(m, fragment(m, i, j), 2) - hmat(2, 2)
                            CALL pbc_hexagonal(coord_v(m, fragment(m, i, j), :), coord_v(m, fragment(m, i, 1), :), vec, vec_pbc, &
                                               box_all, box_x, box_y, box_z)
                        END IF
                        IF (vec(1)>3.0_dp .AND. vec(2)<-5.2_dp) THEN
                            coord_v(m, fragment(m, i, j), 1) = coord_v(m, fragment(m, i, j), 1) - hmat(1, 2)
                            coord_v(m, fragment(m, i, j), 2) = coord_v(m, fragment(m, i, j), 2) + hmat(2, 2)
                            CALL pbc_hexagonal(coord_v(m, fragment(m, i, j), :), coord_v(m, fragment(m, i, 1), :), vec, vec_pbc, &
                                               box_all, box_x, box_y, box_z)
                        END IF
                        IF (vec(1)<-3.0_dp .AND. vec(2)<-5.2_dp) THEN
                            coord_v(m, fragment(m, i, j), 1) = coord_v(m, fragment(m, i, j), 1) + hmat(1, 2)
                            coord_v(m, fragment(m, i, j), 2) = coord_v(m, fragment(m, i, j), 2) + hmat(2, 2)
                            CALL pbc_hexagonal(coord_v(m, fragment(m, i, j), :), coord_v(m, fragment(m, i, 1), :), vec, vec_pbc, &
                                               box_all, box_x, box_y, box_z)
                        END IF
                        IF (vec(1)>4.9_dp .AND. vec(2)<5.2_dp .AND. vec(2)>-5.2_dp) THEN
                            coord_v(m, fragment(m, i, j), 1) = coord_v(m, fragment(m, i, j), 1) - hmat(1, 1)
                            CALL pbc_hexagonal(coord_v(m, fragment(m, i, j), :), coord_v(m, fragment(m, i, 1), :), vec, vec_pbc, &
                                               box_all, box_x, box_y, box_z)
                        END IF
                        IF (vec(1)<-4.9_dp .AND. vec(2)<5.2_dp .AND. vec(2)>-5.2_dp) THEN
                            coord_v(m, fragment(m, i, j), 1) = coord_v(m, fragment(m, i, j), 1) + hmat(1, 1)
                            CALL pbc_hexagonal(coord_v(m, fragment(m, i, j), :), coord_v(m, fragment(m, i, 1), :), vec, vec_pbc, &
                                               box_all, box_x, box_y, box_z)
                        END IF
                        IF (vec(3)<-4.9_dp) THEN
                            coord_v(m, fragment(m, i, j), 3) = coord_v(m, fragment(m, i, j), 3) + hmat(3, 3)
                            CALL pbc_hexagonal(coord_v(m, fragment(m, i, j), :), coord_v(m, fragment(m, i, 1), :), vec, vec_pbc, &
                                               box_all, box_x, box_y, box_z)
                        END IF
                        IF (vec(3)>5.0_dp) THEN
                            coord_v(m, fragment(m, i, j), 3) = coord_v(m, fragment(m, i, j), 3) - hmat(3, 3)
                            CALL pbc_hexagonal(coord_v(m, fragment(m, i, j), :), coord_v(m, fragment(m, i, 1), :), vec, vec_pbc, &
                                               box_all, box_x, box_y, box_z)
                        END IF
                    END DO
                END DO
            END DO

!!!WRAPPING UP THE FRAGMENTS TOGETHER
            DO m = 1, framecount
                DO i = 1, 20
                    DO j = 1, natom_frag(i)
                        IF (coord_v(m, fragment(m, i, j), 3)>6.0_dp) THEN
                            coord_v(m, fragment(m, i, j), 3) = coord_v(m, fragment(m, i, j), 3) - hmat(3, 3)
                        END IF
                        IF (coord_v(m, fragment(m, i, j), 3)>0.0_dp .AND. coord_v(m, fragment(m, i, j), 3)<2.0_dp) THEN
                            IF (coord_v(m, fragment(m, i, j), 2)<-3.0_dp .AND. coord_v(m, fragment(m, i, j), 1)<-6.2_dp) THEN
                                coord_v(m, fragment(m, i, 1:natom_frag(i)), 1) = coord_v(m, fragment(m, i, 1:natom_frag(i)), 1) &
                                                                                 + hmat(1, 2)
                                coord_v(m, fragment(m, i, 1:natom_frag(i)), 2) = coord_v(m, fragment(m, i, 1:natom_frag(i)), 2) &
                                                                                 + hmat(2, 2)
                            END IF
                        END IF
                        IF (coord_v(m, fragment(m, i, j), 1)<-8.0_dp) THEN
                            coord_v(m, fragment(m, i, 1:natom_frag(i)), 1) = coord_v(m, fragment(m, i, 1:natom_frag(i)), 1) &
                                                                             + hmat(1, 1)
                        END IF
                    END DO
                END DO
            END DO

        ELSEIF (frag_type=='3') THEN
            DO m = 1, framecount
                DO i = 21, mol_num
                    DO j = 2, natom_frag(i)
                        CALL pbc_hexagonal(coord_v(m, fragment(m, i, j), :), coord_v(m, fragment(m, i, 1), :), vec, vec_pbc, &
                                           box_all, box_x, box_y, box_z)
                        IF (vec(1)>3.0_dp .AND. vec(2)>5.2_dp) THEN
                            coord_v(m, fragment(m, i, j), 1) = coord_v(m, fragment(m, i, j), 1) - hmat(1, 2)
                            coord_v(m, fragment(m, i, j), 2) = coord_v(m, fragment(m, i, j), 2) - hmat(2, 2)
                            CALL pbc_hexagonal(coord_v(m, fragment(m, i, j), :), coord_v(m, fragment(m, i, 1), :), vec, vec_pbc, &
                                               box_all, box_x, box_y, box_z)
                        END IF
                        IF (vec(1)<-3.0_dp .AND. vec(2)>5.2_dp) THEN
                            coord_v(m, fragment(m, i, j), 1) = coord_v(m, fragment(m, i, j), 1) + hmat(1, 2)
                            coord_v(m, fragment(m, i, j), 2) = coord_v(m, fragment(m, i, j), 2) - hmat(2, 2)
                            CALL pbc_hexagonal(coord_v(m, fragment(m, i, j), :), coord_v(m, fragment(m, i, 1), :), vec, vec_pbc, &
                                               box_all, box_x, box_y, box_z)
                        END IF
                        IF (vec(1)>3.0_dp .AND. vec(2)<-5.2_dp) THEN
                            coord_v(m, fragment(m, i, j), 1) = coord_v(m, fragment(m, i, j), 1) - hmat(1, 2)
                            coord_v(m, fragment(m, i, j), 2) = coord_v(m, fragment(m, i, j), 2) + hmat(2, 2)
                            CALL pbc_hexagonal(coord_v(m, fragment(m, i, j), :), coord_v(m, fragment(m, i, 1), :), vec, vec_pbc, &
                                               box_all, box_x, box_y, box_z)
                        END IF
                        IF (vec(1)<-3.0_dp .AND. vec(2)<-5.2_dp) THEN
                            coord_v(m, fragment(m, i, j), 1) = coord_v(m, fragment(m, i, j), 1) + hmat(1, 2)
                            coord_v(m, fragment(m, i, j), 2) = coord_v(m, fragment(m, i, j), 2) + hmat(2, 2)
                            CALL pbc_hexagonal(coord_v(m, fragment(m, i, j), :), coord_v(m, fragment(m, i, 1), :), vec, vec_pbc, &
                                               box_all, box_x, box_y, box_z)
                        END IF
                        IF (vec(1)>4.9_dp .AND. vec(2)<5.2_dp .AND. vec(2)>-5.2_dp) THEN
                            coord_v(m, fragment(m, i, j), 1) = coord_v(m, fragment(m, i, j), 1) - hmat(1, 1)
                            CALL pbc_hexagonal(coord_v(m, fragment(m, i, j), :), coord_v(m, fragment(m, i, 1), :), vec, vec_pbc, &
                                               box_all, box_x, box_y, box_z)
                        END IF
                        IF (vec(1)<-4.9_dp .AND. vec(2)<5.2_dp .AND. vec(2)>-5.2_dp) THEN
                            coord_v(m, fragment(m, i, j), 1) = coord_v(m, fragment(m, i, j), 1) + hmat(1, 1)
                            CALL pbc_hexagonal(coord_v(m, fragment(m, i, j), :), coord_v(m, fragment(m, i, 1), :), vec, vec_pbc, &
                                               box_all, box_x, box_y, box_z)
                        END IF
                        IF (vec(3)<-4.9_dp) THEN
                            coord_v(m, fragment(m, i, j), 3) = coord_v(m, fragment(m, i, j), 3) + hmat(3, 3)
                            CALL pbc_hexagonal(coord_v(m, fragment(m, i, j), :), coord_v(m, fragment(m, i, 1), :), vec, vec_pbc, &
                                               box_all, box_x, box_y, box_z)
                        END IF
                        IF (vec(3)>5.0_dp) THEN
                            coord_v(m, fragment(m, i, j), 3) = coord_v(m, fragment(m, i, j), 3) - hmat(3, 3)
                            CALL pbc_hexagonal(coord_v(m, fragment(m, i, j), :), coord_v(m, fragment(m, i, 1), :), vec, vec_pbc, &
                                               box_all, box_x, box_y, box_z)
                        END IF
                    END DO
                END DO
            END DO

!!!WRAPPING UP THE FRAGMENTS TOGETHER
            DO m = 1, framecount
                DO i = 21, mol_num
                    DO j = 1, natom_frag(i)
                        IF (coord_v(m, fragment(m, i, j), 3)>6.0_dp) THEN
                            coord_v(m, fragment(m, i, j), 3) = coord_v(m, fragment(m, i, j), 3) - hmat(3, 3)
                        END IF
                        IF (coord_v(m, fragment(m, i, j), 3)>0.0_dp .AND. coord_v(m, fragment(m, i, j), 3)<2.0_dp) THEN
                            IF (coord_v(m, fragment(m, i, j), 2)<-3.0_dp .AND. coord_v(m, fragment(m, i, j), 1)<-6.2_dp) THEN
                                coord_v(m, fragment(m, i, 1:natom_frag(i)), 1) = coord_v(m, fragment(m, i, 1:natom_frag(i)), 1) &
                                                                                 + hmat(1, 2)
                                coord_v(m, fragment(m, i, 1:natom_frag(i)), 2) = coord_v(m, fragment(m, i, 1:natom_frag(i)), 2) &
                                                                                 + hmat(2, 2)
                            END IF
                        END IF
                        IF (coord_v(m, fragment(m, i, j), 1)<-8.0_dp) THEN
                            coord_v(m, fragment(m, i, 1:natom_frag(i)), 1) = coord_v(m, fragment(m, i, 1:natom_frag(i)), 1) &
                                                                             + hmat(1, 1)
                        END IF
                    END DO
                END DO
            END DO

        END IF

        bec_pbc = 0.0_dp
        coord_shifted = 0.
        l = 0

        IF (system=='2' .AND. type_dipole=='1') THEN
            mol_num = 20
        END IF

        DO m = 1, framecount
            DO i = 1, mol_num
                DO j = 1, natom_frag(i)
                    pox_x = 50.0_dp
                    pox_y = 50.0_dp
                    pox_z = 50.0_dp
                    bec(m, :) = coord_v(m, fragment(m, i, j), :)
                    bec_pbc(m, 1) = bec(m, 1) - pox_x*ANINT((1./pox_x)*bec(m, 1))
                    bec_pbc(m, 2) = bec(m, 2) - pox_y*ANINT((1./pox_y)*bec(m, 2))
                    bec_pbc(m, 3) = bec(m, 3) - pox_z*ANINT((1./pox_z)*bec(m, 3))
                    com(m, i, j, :) = bec_pbc(m, :)*mass_atom(fragment(m, i, j))
                    IF (system=='1') THEN
                        refpoint(m, i, :) = refpoint(m, i, :) + com(m, i, j, :)
                    ELSEIF (type_dipole=='1' .AND. system=='2') THEN
                        refpoint(m, 1, :) = refpoint(m, 1, :) + com(m, i, j, :)  !!!For COM of whole system
                    END IF
                END DO
            END DO
        END DO

        mass_tot_cell = 0.0_dp

!!get total mass of the cell
        DO j = 1, natom
            mass_tot_cell = mass_tot_cell + mass_atom(j)
        END DO

        DO m = 1, framecount
            IF (type_dipole=='1' .AND. system=='2') THEN
                refpoint(m, 1, :) = refpoint(m, 1, :)/mass_tot_cell
            ELSEIF (system=='1') THEN
                DO i = 1, mol_num
                    refpoint(m, i, :) = refpoint(m, i, :)/mass_tot_frag(m, i)
                END DO
            END IF
        END DO

        OPEN (UNIT=60, FILE=TRIM(filename)//'-dipole_result_wholesystem.xyz', STATUS='unknown', IOSTAT=stat)
        DO m = 1, framecount
            WRITE (60, *) natom + 1
            WRITE (60, *)
            DO i = 1, 20
                DO j = 1, natom_frag(i)
                    WRITE (60, *) element(fragment(m, i, j)), coord_v(m, fragment(m, i, j), :)
                END DO
!WRITE(60,*) 'N', refpoint(m,i,:)
            END DO
            WRITE (60, *) 'N', refpoint(m, 1, :)
        END DO
        CLOSE (60)

        OPEN (UNIT=12, FILE=TRIM(filename)//'-dipole_result_Ph.xyz', STATUS='unknown', IOSTAT=stat)
        DO m = 1, framecount
            WRITE (12, *) 288
            WRITE (12, *)
            DO i = 9, 20
                DO j = 1, natom_frag(i)
                    WRITE (12, *) element(fragment(m, i, j)), coord_v(m, fragment(m, i, j), :)
                END DO
                WRITE (12, *) 'N', refpoint(m, i, :)
            END DO
        END DO
        CLOSE (12)

        DEALLOCATE (com)
        PRINT *, mass_atom(fragment(1, 9, 1)), mass_atom(fragment(1, 20, 2)), element(fragment(1, 9, 1)), &
            element(fragment(1, 20, 2))
        PRINT *, mass_atom(10), element(10), 'mass atoms'
    END SUBROUTINE center_mass

!***************************************************************************************************
!***************************************************************************************************

    SUBROUTINE solv_frag_index(natom, coord_v, filename, element, box_all, vec, vec_pbc, &
                               box_x, box_y, box_z, mass_atom, framecount, cell_type, refpoint, natom_frag, fragment, mass_tot_frag)

        INTEGER, INTENT(INOUT)                                       :: natom, framecount
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE, INTENT(INOUT)     :: coord_v
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)         :: mass_atom
        CHARACTER(LEN=40), INTENT(INOUT)                             :: filename, cell_type
        CHARACTER(LEN=2), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)     :: element
        REAL(kind=dp), INTENT(INOUT)                                  :: box_all, box_x, box_y, box_z
        REAL(kind=dp), INTENT(INOUT)                                  :: vec(3), vec_pbc(3)
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE, INTENT(OUT)       :: refpoint
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT)              :: natom_frag
        INTEGER, DIMENSION(:, :, :), ALLOCATABLE, INTENT(INOUT)          :: fragment
        REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE, INTENT(OUT)         :: mass_tot_frag

        CHARACTER(LEN=40)                                           :: length
        INTEGER                                                     :: stat   ! error status of OPEN statements
        INTEGER                                                     :: r, m, p, q, i, j, k, n, l, o, s, t
        REAL(kind=dp)                                                :: bec(5000, 3), bec_pbc(5000, 3)
        REAL(kind=dp)                                                :: hmat(3, 3), sqrt3, acosa, asina, a
        REAL(kind=dp)                                                :: pox_x, pox_y, pox_z
        REAL(kind=dp), DIMENSION(3)                                   :: coord3
        REAL(kind=dp), DIMENSION(:, :, :, :), ALLOCATABLE                 :: com

!ALLOCATE(fragment(framecount,37,32))
        ALLOCATE (refpoint(framecount, 37, 3))
!ALLOCATE(natom_frag(37))
        ALLOCATE (mass_tot_frag(framecount, 37))
        ALLOCATE (com(framecount, 37, 32, 3))
!ALLOCATE(natom_frag(20))

        coord3(1) = 0.00_dp
        coord3(2) = 0.00_dp
        coord3(3) = 0.00_dp
        fragment = 0
        natom_frag = 0
        refpoint = 0.00_dp
        com = 0.0_dp
        mass_tot_frag = 0.0_dp

        sqrt3 = 1.73205080756887729352744634_dp

        a = 0.5_dp*(box_x + box_y)
        acosa = 0.5_dp*a
        asina = sqrt3*acosa
        hmat(1, 1) = a; hmat(1, 2) = acosa; hmat(1, 3) = 0.0_dp
        hmat(2, 1) = 0.0_dp; hmat(2, 2) = asina; hmat(2, 3) = 0.0_dp
        hmat(3, 1) = 0.0_dp; hmat(3, 2) = 0.0_dp; hmat(3, 3) = box_z

        PRINT *, natom, framecount

!Assign fragments for solvent molecules!
        DO m = 1, framecount
            j = 0
!DO i=1,natom
            DO i = 1, 98
                IF (element(i).NE.'O') CYCLE
                IF (ANY(i==fragment(m, :, :))) CYCLE
                j = j + 1
                n = 0
                fragment(m, j, 1) = i
                DO k = 1, 98
                    IF (i==k) CYCLE
                    IF (ANY(k==fragment(m, :, :))) CYCLE
                    IF (element(k).NE.'C') CYCLE
                    CALL pbc_hexagonal(coord_v(m, i, :), coord_v(m, k, :), vec, vec_pbc, box_all, box_x, box_y, box_z)
                    IF (SQRT(DOT_PRODUCT(vec_pbc, vec_pbc))<1.781835_dp) THEN
                        IF (ANY(k==fragment(m, :, :))) CYCLE
                        n = n + 1
                        fragment(m, j, n + 1) = k
                        DO l = 1, 98
                            IF (k==l) CYCLE
                            IF (i==l) CYCLE
                            IF (element(l)=='H') THEN
                                CALL pbc_hexagonal(coord_v(m, k, :), coord_v(m, l, :), vec, vec_pbc, box_all, box_x, box_y, box_z)
                                IF (SQRT(DOT_PRODUCT(vec_pbc, vec_pbc))<1.451835_dp) THEN
                                    IF (ANY(l==fragment(m, :, :))) CYCLE
                                    n = n + 1
                                    fragment(m, j, n + 1) = l
                                END IF
                            END IF
                            IF (element(l)=='C') THEN
                                CALL pbc_hexagonal(coord_v(m, k, :), coord_v(m, l, :), vec, vec_pbc, box_all, box_x, box_y, box_z)
                                IF (SQRT(DOT_PRODUCT(vec_pbc, vec_pbc))<1.781835_dp) THEN
                                    IF (ANY(l==fragment(m, :, :))) CYCLE
                                    n = n + 1
                                    fragment(m, j, n + 1) = l

                                    DO o = 1, 98
                                        IF (l==o) CYCLE
                                        IF (i==o) CYCLE
                                        IF (k==o) CYCLE
                                        IF (element(o)=='H') THEN
                                            CALL pbc_hexagonal(coord_v(m, l, :), coord_v(m, o, :), vec, vec_pbc, &
                                                               box_all, box_x, box_y, box_z)
                                            IF (SQRT(DOT_PRODUCT(vec_pbc, vec_pbc))<1.451835_dp) THEN
                                                IF (ANY(o==fragment(m, :, :))) CYCLE
                                                n = n + 1
                                                fragment(m, j, n + 1) = o
                                            END IF
                                        END IF
                                        IF (element(o)=='O') THEN
                                            CALL pbc_hexagonal(coord_v(m, l, :), coord_v(m, o, :), vec, vec_pbc, &
                                                               box_all, box_x, box_y, box_z)
                                            IF (SQRT(DOT_PRODUCT(vec_pbc, vec_pbc))<1.781835_dp) THEN
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
        DO m = 1, framecount
!print*,m
            DO i = 1, 7
                p = 0
                DO j = 1, 14
                    IF (element(fragment(m, i, j))=='H') CYCLE
                    DO k = 1, natom
                        IF (element(k).NE.'X') CYCLE
                        IF (element(fragment(m, i, j))=='O' .OR. element(fragment(m, i, j))=='C') THEN
                            CALL pbc_hexagonal(coord_v(m, fragment(m, i, j), :), coord_v(m, k, :), vec, vec_pbc, &
                                               box_all, box_x, box_y, box_z)
                            IF (SQRT(DOT_PRODUCT(vec_pbc, vec_pbc))<0.99_dp) THEN
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
        DO m = 1, framecount
            j = 7
!DO i=1,natom
            DO i = 99, natom
                IF (element(i).NE.'B') CYCLE
                IF (ANY(i==fragment(m, :, :))) CYCLE
                j = j + 1
                n = 0
                fragment(m, j, 1) = i
                DO k = 99, natom
                    IF (i==k) CYCLE
                    IF (ANY(k==fragment(m, :, :))) CYCLE
                    IF (element(k).NE.'O') CYCLE
                    CALL pbc_hexagonal(coord_v(m, i, :), coord_v(m, k, :), vec, vec_pbc, box_all, &
                                       box_x, box_y, box_z)
                    IF (SQRT(DOT_PRODUCT(vec_pbc, vec_pbc))<1.551835_dp) THEN
                        IF (ANY(k==fragment(m, :, :))) CYCLE
                        n = n + 1
                        fragment(m, j, n + 1) = k
                        DO l = 99, natom
                            IF (k==l) CYCLE
                            IF (i==l) CYCLE
                            IF (element(l).NE.'B') CYCLE
                            CALL pbc_hexagonal(coord_v(m, k, :), coord_v(m, l, :), vec, vec_pbc, &
                                               box_all, box_x, box_y, box_z)
                            IF (SQRT(DOT_PRODUCT(vec_pbc, vec_pbc))<1.551835_dp) THEN
                                IF (ANY(l==fragment(m, :, :))) CYCLE
                                n = n + 1
                                fragment(m, j, n + 1) = l
                                DO o = 99, natom
                                    IF (l==o) CYCLE
                                    IF (i==o) CYCLE
                                    IF (k==o) CYCLE
                                    IF (element(o).NE.'O') CYCLE
                                    CALL pbc_hexagonal(coord_v(m, l, :), coord_v(m, o, :), vec, vec_pbc, &
                                                       box_all, box_x, box_y, box_z)
                                    IF (SQRT(DOT_PRODUCT(vec_pbc, vec_pbc))<1.551835_dp) THEN
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
        DO m = 1, framecount
            DO i = 8, 19
                n = 0
                DO j = 1, 6
                    IF (element(fragment(m, i, j)).NE.'O') CYCLE
                    DO k = 1, natom
                        IF (j==k) CYCLE
                        IF (element(k).NE.'X') CYCLE
                        CALL pbc_hexagonal(coord_v(m, fragment(m, i, j), :), coord_v(m, k, :), vec, vec_pbc, &
                                           box_all, box_x, box_y, box_z)
                        IF (SQRT(DOT_PRODUCT(vec_pbc, vec_pbc))<0.7_dp) THEN
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
        DO m = 1, framecount
            j = 19
            DO i = 99, natom
                IF (element(i).NE.'H') CYCLE
                IF (ANY(i==fragment(m, :, :))) CYCLE
                j = j + 1
!        print*,j,'j',element(i),i
                n = 0
                fragment(m, j, 1) = i
                DO k = 99, natom
                    IF (i==k) CYCLE
                    ! IF (ANY(k==fragment(m,:,:))) CYCLE
                    IF (element(k).NE.'C') CYCLE
                    CALL pbc_hexagonal(coord_v(m, i, :), coord_v(m, k, :), &
                                       vec, vec_pbc, box_all, box_x, box_y, box_z)
                    IF (SQRT(DOT_PRODUCT(vec_pbc, vec_pbc))<1.31835_dp) THEN
                        IF (ANY(k==fragment(m, :, :))) CYCLE
                        n = n + 1
                        fragment(m, j, n + 1) = k
!        print*,k,element(k),i,element(i), 'k ve i'
                        DO l = 99, natom
                            !IF (ANY(l==fragment(m,:,:))) CYCLE
                            IF (element(l).NE.'C') CYCLE
                            CALL pbc_hexagonal(coord_v(m, k, :), coord_v(m, l, :), &
                                               vec, vec_pbc, box_all, box_x, box_y, box_z)
                            IF (SQRT(DOT_PRODUCT(vec_pbc, vec_pbc))<1.55835_dp) THEN
                                IF (ANY(l==fragment(m, :, :))) CYCLE
                                n = n + 1
                                fragment(m, j, n + 1) = l
!        print*,l,element(l),k,element(k),'l ve k'
                                DO o = 99, natom
                                    ! IF (ANY(o==fragment(m,:,:))) CYCLE
                                    IF (element(o)=='X') CYCLE
                                    CALL pbc_hexagonal(coord_v(m, l, :), coord_v(m, o, :), &
                                                       vec, vec_pbc, box_all, box_x, box_y, box_z)
                                    IF (SQRT(DOT_PRODUCT(vec_pbc, vec_pbc))<1.55835_dp) THEN
                                        IF (ANY(o==fragment(m, :, :))) CYCLE
                                        n = n + 1
                                        fragment(m, j, n + 1) = o
                                        ! print*,o,element(o),l,element(l),'o ve l'
                                        DO p = 99, natom
                                            ! IF (ANY(p==fragment(m,:,:))) CYCLE
                                            IF (element(p)=='X') CYCLE
                                            CALL pbc_hexagonal(coord_v(m, o, :), coord_v(m, p, :), &
                                                               vec, vec_pbc, box_all, box_x, box_y, box_z)
                                            IF (SQRT(DOT_PRODUCT(vec_pbc, vec_pbc))<1.55835_dp) THEN
                                                IF (ANY(p==fragment(m, :, :))) CYCLE
                                                n = n + 1
                                                fragment(m, j, n + 1) = p
                                                !        print*,p,element(p),o,element(o),'p ve o'
                                                DO r = 99, natom
                                                    !IF (ANY(r==fragment(m,:,:))) CYCLE
                                                    IF (element(r)=='X') CYCLE
                                                    CALL pbc_hexagonal(coord_v(m, p, :), coord_v(m, r, :), &
                                                                       vec, vec_pbc, box_all, box_x, box_y, box_z)
                                                    IF (SQRT(DOT_PRODUCT(vec_pbc, vec_pbc))<1.55835_dp) THEN
                                                        IF (ANY(r==fragment(m, :, :))) CYCLE
                                                        n = n + 1
                                                        fragment(m, j, n + 1) = r
                                                        !        print*,r,element(r),p,element(p),'r ve p'
                                                        DO s = 99, natom
                                                            !IF (ANY(s==fragment(m,:,:))) CYCLE
                                                            IF (element(s)=='X') CYCLE
                                                            CALL pbc_hexagonal(coord_v(m, r, :), coord_v(m, s, :), &
                                                                               vec, vec_pbc, box_all, box_x, box_y, box_z)
                                                            IF (SQRT(DOT_PRODUCT(vec_pbc, vec_pbc))<1.55835_dp) THEN
                                                                IF (ANY(s==fragment(m, :, :))) CYCLE
                                                                n = n + 1
                                                                fragment(m, j, n + 1) = s
                                                                !        print*,s,element(s),r,element(r),'s ve r'
                                                                DO t = 99, natom
                                                                    !IF (ANY(s==fragment(m,:,:))) CYCLE
                                                                    IF (element(t)=='X') CYCLE
                                                                    CALL pbc_hexagonal(coord_v(m, s, :), coord_v(m, t, :), &
                                                                                       vec, vec_pbc, box_all, box_x, box_y, box_z)
                                                                    IF (SQRT(DOT_PRODUCT(vec_pbc, vec_pbc))<1.55835_dp) THEN
                                                                        IF (ANY(t==fragment(m, :, :))) CYCLE
                                                                        n = n + 1
                                                                        fragment(m, j, n + 1) = t
                                                                        !         print*,t,element(t),s,element(s),'t ve s'
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

        DO m = 1, framecount
            DO i = 20, 37
                n = 0
                DO j = 1, 10
                    IF (element(fragment(m, i, j)).NE.'C') CYCLE
                    DO k = 1, natom
                        IF (j==k) CYCLE
                        IF (element(k).NE.'X') CYCLE
                        CALL pbc_hexagonal(coord_v(m, fragment(m, i, j), :), coord_v(m, k, :), vec, vec_pbc, &
                                           box_all, box_x, box_y, box_z)
                        IF (SQRT(DOT_PRODUCT(vec_pbc, vec_pbc))<1.0_dp) THEN
                            IF (ANY(k==fragment(m, i, :))) CYCLE
                            n = n + 1
                            fragment(m, i, n + 10) = k
                        END IF
                    END DO
                END DO
            END DO
        END DO

!!get natom_frag
        DO i = 1, framecount
            DO j = 1, 37
                natom_frag(j) = COUNT(fragment(i, j, :).NE.0)
            END DO
        END DO

!!get total mass
        DO m = 1, framecount
            DO i = 1, 37
                DO j = 1, natom_frag(i)
                    IF (element(natom_frag(i))=='X') CYCLE
                    mass_tot_frag(m, i) = mass_tot_frag(m, i) + mass_atom(fragment(m, i, j))
                END DO
            END DO
        END DO

!!Find center of mass

        DO m = 1, framecount
            DO i = 1, 37
                DO j = 2, natom_frag(i)
                    CALL pbc_hexagonal(coord_v(m, fragment(m, i, j), :), coord_v(m, fragment(m, i, 1), :), &
                                       vec, vec_pbc, box_all, box_x, box_y, box_z)
                    ! PRINT*,SQRT(DOT_PRODUCT(vec,vec)),SQRT(DOT_PRODUCT(vec_pbc,vec_pbc)),i,j,'0'
                    IF (vec(1)>3.0_dp .AND. vec(2)>5.2_dp) THEN
                        coord_v(m, fragment(m, i, j), 1) = coord_v(m, fragment(m, i, j), 1) - hmat(1, 2)
                        coord_v(m, fragment(m, i, j), 2) = coord_v(m, fragment(m, i, j), 2) - hmat(2, 2)
                        CALL pbc_hexagonal(coord_v(m, fragment(m, i, j), :), coord_v(m, fragment(m, i, 1), :), &
                                           vec, vec_pbc, box_all, box_x, box_y, box_z)
                        ! PRINT*,SQRT(DOT_PRODUCT(vec,vec)),SQRT(DOT_PRODUCT(vec_pbc,vec_pbc)),i,j,'3'
                    END IF
                    IF (vec(1)<-3.0_dp .AND. vec(2)>5.2_dp) THEN
                        coord_v(m, fragment(m, i, j), 1) = coord_v(m, fragment(m, i, j), 1) + hmat(1, 2)
                        coord_v(m, fragment(m, i, j), 2) = coord_v(m, fragment(m, i, j), 2) - hmat(2, 2)
                        CALL pbc_hexagonal(coord_v(m, fragment(m, i, j), :), coord_v(m, fragment(m, i, 1), :), &
                                           vec, vec_pbc, box_all, box_x, box_y, box_z)
                        ! PRINT*,SQRT(DOT_PRODUCT(vec,vec)),SQRT(DOT_PRODUCT(vec_pbc,vec_pbc)),i,j,'4'
                    END IF
                    IF (vec(1)>3.0_dp .AND. vec(2)<-5.2_dp) THEN
                        coord_v(m, fragment(m, i, j), 1) = coord_v(m, fragment(m, i, j), 1) - hmat(1, 2)
                        coord_v(m, fragment(m, i, j), 2) = coord_v(m, fragment(m, i, j), 2) + hmat(2, 2)
                        CALL pbc_hexagonal(coord_v(m, fragment(m, i, j), :), coord_v(m, fragment(m, i, 1), :), &
                                           vec, vec_pbc, box_all, box_x, box_y, box_z)
                        ! PRINT*,SQRT(DOT_PRODUCT(vec,vec)),SQRT(DOT_PRODUCT(vec_pbc,vec_pbc)),i,j,'5'
                    END IF
                    IF (vec(1)<-3.0_dp .AND. vec(2)<-5.2_dp) THEN
                        coord_v(m, fragment(m, i, j), 1) = coord_v(m, fragment(m, i, j), 1) + hmat(1, 2)
                        coord_v(m, fragment(m, i, j), 2) = coord_v(m, fragment(m, i, j), 2) + hmat(2, 2)
                        CALL pbc_hexagonal(coord_v(m, fragment(m, i, j), :), coord_v(m, fragment(m, i, 1), :), &
                                           vec, vec_pbc, box_all, box_x, box_y, box_z)
                        ! PRINT*,SQRT(DOT_PRODUCT(vec,vec)),SQRT(DOT_PRODUCT(vec_pbc,vec_pbc)),i,j,'6'
                    END IF
                    IF (vec(1)>4.9_dp .AND. vec(2)<5.2_dp .AND. vec(2)>-5.2_dp) THEN
                        coord_v(m, fragment(m, i, j), 1) = coord_v(m, fragment(m, i, j), 1) - hmat(1, 1)
                        CALL pbc_hexagonal(coord_v(m, fragment(m, i, j), :), coord_v(m, fragment(m, i, 1), :), &
                                           vec, vec_pbc, box_all, box_x, box_y, box_z)
                        ! PRINT*,SQRT(DOT_PRODUCT(vec,vec)),SQRT(DOT_PRODUCT(vec_pbc,vec_pbc)),i,j,'1'
                    END IF
                    IF (vec(1)<-4.9_dp .AND. vec(2)<5.2_dp .AND. vec(2)>-5.2_dp) THEN
                        coord_v(m, fragment(m, i, j), 1) = coord_v(m, fragment(m, i, j), 1) + hmat(1, 1)
                        CALL pbc_hexagonal(coord_v(m, fragment(m, i, j), :), coord_v(m, fragment(m, i, 1), :), &
                                           vec, vec_pbc, box_all, box_x, box_y, box_z)
                        ! PRINT*,SQRT(DOT_PRODUCT(vec,vec)),SQRT(DOT_PRODUCT(vec_pbc,vec_pbc)),i,j,'2'
                    END IF
                    IF (vec(3)<-4.9_dp) THEN
                        coord_v(m, fragment(m, i, j), 3) = coord_v(m, fragment(m, i, j), 3) + hmat(3, 3)
                        CALL pbc_hexagonal(coord_v(m, fragment(m, i, j), :), coord_v(m, fragment(m, i, 1), :), &
                                           vec, vec_pbc, box_all, box_x, box_y, box_z)
                        ! PRINT*,SQRT(DOT_PRODUCT(vec,vec)),SQRT(DOT_PRODUCT(vec_pbc,vec_pbc)),i,j,'7'
                    END IF
                    IF (vec(3)>5.0_dp) THEN
                        coord_v(m, fragment(m, i, j), 3) = coord_v(m, fragment(m, i, j), 3) - hmat(3, 3)
                        CALL pbc_hexagonal(coord_v(m, fragment(m, i, j), :), coord_v(m, fragment(m, i, 1), :), &
                                           vec, vec_pbc, box_all, box_x, box_y, box_z)
                        ! PRINT*,SQRT(DOT_PRODUCT(vec,vec)),SQRT(DOT_PRODUCT(vec_pbc,vec_pbc)),i,j,'8'
                    END IF
                END DO
            END DO
        END DO

        bec_pbc = 0.0_dp

        DO i = 1, 37
            DO j = 1, natom_frag(i)
                DO m = 1, framecount
                    pox_x = 40.0_dp
                    pox_y = 40.0_dp
                    pox_z = 40.0_dp
                    bec(m, :) = coord_v(m, fragment(m, i, j), :)
                    bec_pbc(m, 1) = bec(m, 1) - pox_x*ANINT((1./pox_x)*bec(m, 1))
                    bec_pbc(m, 2) = bec(m, 2) - pox_y*ANINT((1./pox_y)*bec(m, 2))
                    bec_pbc(m, 3) = bec(m, 3) - pox_z*ANINT((1./pox_z)*bec(m, 3))
                    com(m, i, j, :) = bec_pbc(m, :)*mass_atom(fragment(m, i, j))
                    refpoint(m, i, :) = refpoint(m, i, :) + com(m, i, j, :)
                END DO
            END DO
        END DO

        DO m = 1, framecount
            DO i = 1, 37
                refpoint(m, i, :) = refpoint(m, i, :)/mass_tot_frag(m, i)
            END DO
        END DO

        OPEN (UNIT=60, FILE=TRIM(filename)//'-dipole_result_7.xyz', STATUS='unknown', IOSTAT=stat)
        DO m = 1, framecount
            WRITE (60, *) natom + 37
            WRITE (60, *)
            DO i = 1, 37
                DO j = 1, natom_frag(i)
                    WRITE (60, *) element(fragment(m, i, j)), coord_v(m, fragment(m, i, j), :)
                END DO
                WRITE (60, *) 'N', refpoint(m, i, :)
            END DO
        END DO
        CLOSE (60)

        DEALLOCATE (com)

    END SUBROUTINE solv_frag_index
!***************************************************************************************************
!***************************************************************************************************

    SUBROUTINE wannier_frag(natom_frag, filename, natom, element, coord_v, box_all, box_x, box_y, box_z, vec, vec_pbc, &
                            dipole, refpoint, fragment, framecount, mass_tot_frag, mol_num, system, type_dipole, charge, &
                            mass_tot_cell)

        CHARACTER(LEN=40), INTENT(IN)                                :: filename, system, type_dipole
        CHARACTER(LEN=2), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)     :: element
        INTEGER, INTENT(INOUT)                                       :: natom, framecount, mol_num
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT)              :: natom_frag
        INTEGER, DIMENSION(:, :, :), ALLOCATABLE, INTENT(INOUT)          :: fragment
        REAL(kind=dp), INTENT(INOUT)                                  :: box_all, box_x, box_y, box_z, mass_tot_cell
        REAL(kind=dp), DIMENSION(3), INTENT(INOUT)                     :: vec, vec_pbc
        REAL(kind=dp), DIMENSION(5000, 3)                              :: bec, bec_pbc
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)         :: charge
        REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE, INTENT(INOUT)       :: mass_tot_frag
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE, INTENT(INOUT)     :: coord_v, refpoint
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE, INTENT(OUT)       :: dipole

        CHARACTER(LEN=40)                                           :: length, filename_dip
        INTEGER                                                     :: stat   ! error status of OPEN statements
        INTEGER                                                     :: i, j, l, n, m
        REAL(kind=dp)                                                :: dist!,mass_tot(20)
        REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE                     :: mass
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE                       :: coord2
        REAL(kind=dp)                                                :: pox_x, pox_y, pox_z, pox_all
        REAL(kind=dp)                                                :: hmat(3, 3), sqrt3, acosa, asina, a

        ALLOCATE (dipole(framecount, mol_num, 3), coord2(1))

        dist = 1.2_dp
        dipole = 0.0_dp

        pox_x = 40.0_dp
        pox_y = 40.0_dp
        pox_z = 40.0_dp
        pox_all = 40.0_dp

        sqrt3 = 1.73205080756887729352744634_dp
        a = 0.5_dp*(box_x + box_y)
        acosa = 0.5_dp*a
        asina = sqrt3*acosa
        hmat(1, 1) = a; hmat(1, 2) = acosa; hmat(1, 3) = 0.0_dp
        hmat(2, 1) = 0.0_dp; hmat(2, 2) = asina; hmat(2, 3) = 0.0_dp
        hmat(3, 1) = 0.0_dp; hmat(3, 2) = 0.0_dp; hmat(3, 3) = box_z
        PRINT *, mol_num, 'check mol num'
        DO m = 1, framecount
            DO i = 1, mol_num
                DO j = 1, natom_frag(i)
                    IF (system=='1') THEN
                        CALL pbc_hexagonal(coord_v(m, fragment(m, i, j), :), refpoint(m, i, :), vec, vec_pbc, &
                                           box_all, box_x, box_y, box_z) !! If fragment approach!!!
                        dipole(m, i, :) = dipole(m, i, :) + vec_pbc*1.889725989_dp*charge(fragment(m, i, j))
                    ELSEIF (type_dipole=='1' .AND. system=='2') THEN
                        CALL pbc_orthorombic(coord_v(m, fragment(m, i, j), :), refpoint(m, 1, :), vec, vec_pbc, &
                                             pox_all, pox_x, pox_y, pox_z) !! IF COM is for whole supercel!!!
                        dipole(m, 1, :) = dipole(m, 1, :) + vec_pbc*1.889725989_dp*charge(fragment(m, i, j))
                    END IF
                END DO
            END DO
        END DO

        PRINT *, mass_tot_cell, 'mass tot cell'

        DO m = 1, framecount
            IF (system=='2' .AND. type_dipole=='1') THEN
                IF (dipole(m, 1, 1)>120) THEN
                    dipole(m, 1, 1) = dipole(m, 1, 1) - (hmat(1, 1)*3*1.889725989)
                END IF
!  dipole(m,1,1)=dipole(m,1,1)-(hmat(1,1)*3*1.889725989)
                dipole(m, 1, 2) = dipole(m, 1, 2) + (hmat(2, 2)*4*1.889725989)!+(hmat(1,2)*3*1.889725989)
                dipole(m, 1, :) = REAL((dipole(m, 1, :)*2.54174622741_dp)/mass_tot_cell, kind=dp) !!Dividing by total mass, change it later
            ELSEIF (system=='1') THEN
                DO i = 1, mol_num !! define something for this
                    dipole(m, i, :) = REAL((dipole(m, i, :)*2.54174622741_dp)/mass_tot_frag(m, i), kind=dp) !!Dividing by total mass, change it later
                END DO
            END IF
        END DO

!DO j=21,44 !! define something for this
! WRITE(filename_dip, '(i0)') j
! OPEN(UNIT=68,FILE='COF-1_dipoles-'//trim(filename_dip)//'.txt',STATUS='unknown',IOSTAT=stat)
! DO m=1,framecount
!  WRITE(68,*) m,dipole(m,j,1:3)!,vec_pbc(3)
! ENDDO
!ENDDO
!CLOSE(68)

!OPEN(UNIT=68,FILE='dipole_result_final.xyz',STATUS='unknown',IOSTAT=stat)
!DO m=1,framecount
!DO j=1,20
        !  WRITE(68,'(2X,A15,4X,I2,6X,F20.12)') "net dipole",j, SQRT(DOT_PRODUCT(dipole(m,j,:),dipole(m,j,:)))*2.54174622741_dp

! ENDDO
! ENDDO
!CLOSE(68)

!OPEN(UNIT=61,FILE='dipole_result_vmd.xyz',STATUS='unknown',IOSTAT=stat)
!DO m=1,framecount
! DO j=1,20
        ! WRITE(61,'(2X,A15,4X,I2,6X,3F20.12)') "center of mass", j, refpoint(m,j,:)

        !WRITE(61,'(2X,A15,4X,I2,6X,3F20.12)') "com+net dipole", j,  refpoint(m,j,:)+(dipole(m,j,:)*2.541746227414447_dp)

! ENDDO
! ENDDO
!CLOSE(61)

!OPEN(UNIT=51,FILE='dipole_result_vmd_final.xyz',STATUS='unknown',IOSTAT=stat)
!DO m=1,framecount
!WRITE(51,*) framecount
!WRITE(51,*)
! DO j=1,20
        ! WRITE(51,*) "draw arrow     ", "{",refpoint(m,j,:)&
        !        ,"}      ","{",refpoint(m,j,:)+(dipole(m,j,:)*2.541746227414447_dp),"}"
! ENDDO
!ENDDO
!CLOSE(51)
        DEALLOCATE (mass_tot_frag)
        DEALLOCATE (refpoint)
    END SUBROUTINE wannier_frag

!***************************************************************************************************
!***************************************************************************************************

    SUBROUTINE wannier(element, filename, natom, box_all, box_x, box_y, box_z, vec, vec_pbc, debye, mol_num, &
                       periodic, mass_tot, framecount, mass_atom, coord_v, dip)

        CHARACTER(LEN=40), INTENT(INOUT)                                  :: filename, periodic
        INTEGER, INTENT(INOUT)                                            :: natom, framecount, mol_num
        CHARACTER(LEN=2), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)          :: element
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE, INTENT(INOUT)          :: coord_v
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE, INTENT(OUT)            :: dip
        REAL(kind=dp), INTENT(INOUT)                                       :: box_all, box_x, box_y, box_z
        REAL(kind=dp), INTENT(INOUT)                                       :: debye, mass_tot
        REAL(kind=dp), DIMENSION(3), INTENT(INOUT)                          :: vec, vec_pbc
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)              :: mass_atom

        INTEGER                                                          :: stat                            ! error status of OPEN statements
        INTEGER                                                          :: k, i, j, l, n
        REAL(kind=dp)                                                     :: dist, coord1(3)
        CHARACTER(LEN=40)                                                :: length, chara
        REAL, DIMENSION(:), ALLOCATABLE                                   :: dip_tot !--> in your case of 1 molecule the dimension is one, i.e., it is a simple real
        REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE                          :: refpoint2
        REAL, DIMENSION(:, :, :), ALLOCATABLE                               :: com

        dist = 1.2_dp
        coord1(:) = 0.00_dp

        ALLOCATE (dip(framecount, mol_num, 3))
        ALLOCATE (dip_tot(mol_num))
        ALLOCATE (refpoint2(framecount, 3))
        ALLOCATE (com(framecount, natom, 3))

        dip = 0.0_dp
        l = 0

        IF (periodic=='y' .OR. periodic=='yes') THEN

            DO i = 1, natom
                IF (element(i).NE."O") CYCLE
                l = l + 1
                DO j = 1, natom
                    IF (element(j)=="O") CYCLE
                    DO k = 1, framecount
                        CALL pbc_orthorombic(coord_v(k, j, :), coord_v(k, i, :), vec, vec_pbc, box_all, box_x, box_y, box_z)

                        IF (SQRT(DOT_PRODUCT(vec_pbc, vec_pbc))>dist) CYCLE
                        IF (element(j)=="X") THEN
                            dip(k, l, :) = REAL(dip(k, l, :) + vec_pbc*1.889725989_dp*(-2.0_dp)/debye, kind=dp)
                        ELSEIF (element(j)=="H") THEN
                            dip(k, l, :) = REAL(dip(k, l, :) + vec_pbc*1.889725989_dp*(1.0_dp)/debye, kind=dp)
                        END IF
                    END DO
                    !dip_tot(l)=SQRT(dipole2(l,1)**2+dipole2(l,2)**2+dipole2(l,3)**2)
                END DO
            END DO

        ELSE IF (periodic=='n' .OR. periodic=='no') THEN

            refpoint2 = 0.0_dp

            DO j = 1, natom
                DO k = 1, framecount

                    CALL pbc_orthorombic(coord_v(k, j, :), coord1, vec, vec_pbc, box_all, box_x, box_y, box_z)

                    com(k, j, :) = vec_pbc*mass_atom(j)
                    refpoint2(k, :) = refpoint2(k, :) + com(k, j, :)
                END DO
            END DO

            refpoint2 = REAL(refpoint2/mass_tot, kind=dp)

            DO j = 1, natom
                l = 1
                DO k = 1, framecount
                    CALL pbc_orthorombic(coord_v(k, j, :), refpoint2(k, :), vec, vec_pbc, box_all, box_x, box_y, box_z)

                    IF (element(j)=="X") THEN
                        dip(k, l, :) = REAL(dip(k, l, :) + vec_pbc*1.889725989_dp*(-2.0_dp)/debye, kind=dp)
                    ELSEIF (element(j)=="H") THEN
                        dip(k, l, :) = REAL(dip(k, l, :) + vec_pbc*1.889725989_dp*(1.0_dp)/debye, kind=dp)
                    ELSEIF (element(j)=="O") THEN
                        dip(k, l, :) = REAL(dip(k, l, :) + vec_pbc*1.889725989_dp*(6.0_dp)/debye, kind=dp)
                    ELSEIF (element(j)=="B") THEN
                        dip(k, l, :) = REAL(dip(k, l, :) + vec_pbc*1.889725989_dp*(3.0_dp)/debye, kind=dp)
                    END IF
                    !dip_tot(l)=SQRT(dipole2(l,1)**2+dipole2(l,2)**2+dipole2(l,3)**2)
                END DO
            END DO
        END IF

        OPEN (UNIT=60, FILE='dipole_result.xyz', STATUS='unknown', IOSTAT=stat)

        DO i = 1, framecount
            WRITE (60, *) framecount
            WRITE (60, *) mol_num

            DO j = 1, mol_num
                WRITE (60, *) j, dip(i, j, :)
            END DO
        END DO
        CLOSE (60)

        DEALLOCATE (dip_tot, com, refpoint2)

    END SUBROUTINE wannier

END MODULE dipole_calc

