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

!> @brief Module containing all procedures involving calculation of dipole moments
!>        and correcting dipole moment jumps
MODULE dipole_calc

    USE kinds, ONLY: dp, str_len
    USE constants, ONLY: debye, bohr2ang, pi
    USE ISO_FORTRAN_ENV, ONLY: output_unit, error_unit
    USE vib_types, ONLY: global_settings, systems, molecular_dynamics, static, dipoles, fragment_group
    USE cell_types, ONLY: build_hmat, pbc, pbc_vec, invert3x3, determinant3x3

    USE OMP_LIB

    IMPLICIT NONE

    PUBLIC :: compute_dipole, check_jumps, assign_wannier, compute_dipole_frag, check_lattice_parameters
    PRIVATE

CONTAINS

    !> @brief Removes discontinuities ("jumps") in dipole moment trajectories due to
    !>        polarization quantum crossings in periodic systems.
    !>
    !> This subroutine corrects the spurious discontinuities in the dipole moment
    !> trajectory that arise from periodic boundary conditions in systems where
    !> polarization is defined modulo a polarization quantum. It does the correction
    !> by subtracting integer multiples of polarization quantum from each dipole array.
    !>
    !> @param[inout] sys     -- System data structure (provides cell and frame information)
    !> @param[in,out] natom  -- Number of atoms for cartesian coordinates, or just `1`
    !>                          if the quantity (e.g. dipole or polarizability) refers
    !>                          to the entire system.
    !> @param[inout] dipole  -- 3D array of dipole moments for each frame (nframes × 1 × 3)
    !>
    SUBROUTINE check_jumps(dipole, natom, sys)
        TYPE(systems), INTENT(INOUT) :: sys
        INTEGER, INTENT(IN) :: natom  ! (nframes,1,3)
        REAL(dp), DIMENSION(:, :, :), INTENT(INOUT) :: dipole  ! (nframes,1,3)

        INTEGER :: m, i, j, stat
        REAL(dp) :: hmat(3, 3)
        REAL(dp) :: pol_quantum(3)

        ! --- lattice
        CALL build_hmat(sys, hmat)

   !!!! Calculate the polarization quantum
        DO i = 1, 3
            pol_quantum(i) = SUM(hmat(i, :))/(bohr2ang*debye)  ! e·Bohr
        END DO

   !!Subtract multiples of polarization quantum
        DO j = 1, natom
            DO i = 1, 3
                DO m = 2, sys%framecount
                    dipole(m, j, i) = dipole(m, j, i) - NINT((dipole(m, j, i) - dipole(m - 1, j, i))/pol_quantum(i))*pol_quantum(i)
                END DO
            END DO
        END DO

    END SUBROUTINE check_jumps

!******************************************************************************************************************!
!******************************************************************************************************************!

    !> @brief Dynamically assigns Wannier centers to individual fragments in each MD frame based on
    !>        a distance criteria.
    !>
    !> This subroutine extends each fragment definition by identifying Wannier centers
    !> that are spatially close to atoms of the fragment. If Wannier center is also close to another
    !> fragment, it is not considered (e.g. Wannier centers on bonds between two fragments).
    !>
    !> @param[inout] sys  -- System data structure (provides fragment definitions, atomic
    !>                       elements, and cell parameters. The updated fragment lists with
    !>                       Wannier centers are stored in `sys%fragments%type_frag(:)%fragment_frame`.)
    !> @param[in]    md   -- MD data structure containing cartesian coordinates for all frames.
    !>
    SUBROUTINE assign_wannier(sys, md)

        USE kinds, ONLY: dp
        TYPE(systems), INTENT(INOUT)         :: sys
        TYPE(molecular_dynamics), INTENT(IN) :: md

        INTEGER :: i, j, m, i_group, l, i_atom, frag_count, natom_frag
        REAL(dp) :: hmat(3, 3), h_inv(3, 3), dr(3), dr2(3), dist
        INTEGER, ALLOCATABLE :: frag_atoms_frame(:)
        REAL(kind=8), allocatable :: coords(:, :), distances(:)
        INTEGER, allocatable :: n_new_wanniers(:, :)
        INTEGER, ALLOCATABLE :: new_wannier_centers(:, :, :)
        INTEGER :: n_wannier_tot, n_frag_atoms_tot, n_other_tot
        INTEGER, ALLOCATABLE :: idx_wannier_centers(:), idx_fragment_atoms(:), idx_other_atoms(:)
        INTEGER, ALLOCATABLE :: map_idx_frag_to_group(:), map_idx_frag_to_frag(:)
        INTEGER :: i_wannier, max_nfrag, my_group, my_fragment
        LOGICAL :: belongs_to_exactly_one

        !Cut-off distance between the atom and the Wannier center  
        dist = 1.1_dp 

        ! --- lattice
        CALL build_hmat(sys, hmat)
        CALL invert3x3(hmat, h_inv)

        ! get the indices of wannier centers 
        ALLOCATE(idx_wannier_centers(sys%natom))
        n_wannier_tot = 0
        DO i_atom = 1, sys%natom
            IF (sys%element(i_atom) == 'X') THEN
                ! count number of wannier centers
                n_wannier_tot = n_wannier_tot + 1
                idx_wannier_centers(n_wannier_tot) = i_atom
            END IF
        END DO

        ! get the indices of fragment atoms
        ALLOCATE (idx_fragment_atoms(sys%natom))
        ALLOCATE (map_idx_frag_to_group(sys%natom))
        ALLOCATE (map_idx_frag_to_frag(sys%natom))
        n_frag_atoms_tot = 0
        DO i_group = 1, sys%fragments%ngroup
            frag_count = sys%fragments%type_frag(i_group)%nfrag
            DO j = 1, frag_count
                natom_frag = SIZE(sys%fragments%type_frag(i_group)%fragment(j)%frag_atoms)
                DO i = 1, natom_frag
                    n_frag_atoms_tot = n_frag_atoms_tot + 1
                    idx_fragment_atoms(n_frag_atoms_tot) = sys%fragments%type_frag(i_group)%fragment(j)%frag_atoms(i)
                    map_idx_frag_to_group(n_frag_atoms_tot) = i_group
                    map_idx_frag_to_frag(n_frag_atoms_tot) = j
                END DO
            END DO
        END DO

        ! get the indices of other atoms
        ALLOCATE(idx_other_atoms(sys%natom))
        n_other_tot = 0
        DO i_atom = 1, sys%natom
            IF (sys%element(i_atom) /= 'X') THEN
                IF (.NOT. ANY(i_atom==idx_fragment_atoms(1:n_frag_atoms_tot))) THEN
                    n_other_tot = n_other_tot + 1
                    idx_other_atoms(n_other_tot) = i_atom
                END IF
            END IF
        END DO

       max_nfrag = 0
       DO i_group = 1, sys%fragments%ngroup
            frag_count = sys%fragments%type_frag(i_group)%nfrag
            max_nfrag = MAX(max_nfrag, frag_count)
            ! --- allocate storage for augmented fragment lists
            IF (.NOT. ALLOCATED(sys%fragments%type_frag(i_group)%fragment_frame)) THEN
                ALLOCATE (sys%fragments%type_frag(i_group)%fragment_frame(sys%framecount, frag_count))
            END IF
        END DO

        ALLOCATE (coords(3, n_frag_atoms_tot + n_other_tot))
        ALLOCATE (new_wannier_centers(sys%natom, max_nfrag, sys%fragments%ngroup))
        ALLOCATE (n_new_wanniers(max_nfrag, sys%fragments%ngroup))
        ALLOCATE (distances(max(n_frag_atoms_tot, n_other_tot)))

        ! --- loop over frames
        DO i = 1, sys%framecount
            ! fist tabulate fragment atoms 
            DO l = 1, n_frag_atoms_tot
                coords(:, l) = md%coord_v(i, idx_fragment_atoms(l), :)
            END DO
            ! next the other atoms
            DO l = 1, n_other_tot
                coords(:, n_frag_atoms_tot + l) = md%coord_v(i, idx_other_atoms(l), :)
            END DO

            ! now go over every wannier center and look to which fragment it belongs
            n_new_wanniers(:, :) = 0
            DO i_wannier = 1, n_wannier_tot
                i_atom = idx_wannier_centers(i_wannier)
                
                ! sort out wannier centers that are too close to other atoms
                CALL pbc_vec(n_other_tot, coords(:, n_frag_atoms_tot+1:), &
                     md%coord_v(i, i_atom, :), sys, distances)
                IF (ANY(distances(1:n_other_tot)<dist)) CYCLE

                ! now check to which fragment it belongs
                CALL pbc_vec(n_frag_atoms_tot, coords(:, 1:n_frag_atoms_tot), &
                     md%coord_v(i, i_atom, :), sys, distances)
                
                belongs_to_exactly_one = .false.
                my_fragment = -1
                my_group = -1
                DO l = 1, n_frag_atoms_tot
                    IF (distances(l)<dist) THEN
                        IF (my_fragment == -1) THEN
                            my_fragment = map_idx_frag_to_frag(l)
                            my_group = map_idx_frag_to_group(l)
                            belongs_to_exactly_one = .true.
                        ELSEIF (my_fragment /= map_idx_frag_to_frag(l) .OR. my_group /= map_idx_frag_to_group(l)) THEN
                            ! belongs to more than one fragment, skip
                            belongs_to_exactly_one = .false.
                            EXIT
                        END IF
                    END IF
                END DO

                ! now store the new wannier center if it belongs to exactly one fragment
                IF (belongs_to_exactly_one) THEN
                    n_new_wanniers(my_fragment, my_group) = n_new_wanniers(my_fragment, my_group) + 1
                    new_wannier_centers(n_new_wanniers(my_fragment, my_group), my_fragment, my_group) = i_atom
                END IF
            END DO

            ! now store it into the type
             ! --- loop over number of fragment sections
            DO i_group = 1, sys%fragments%ngroup
                frag_count = sys%fragments%type_frag(i_group)%nfrag
                ! --- loop over the number of atom lists in each fragment section
                DO j = 1, frag_count
                    natom_frag = SIZE(sys%fragments%type_frag(i_group)%fragment(j)%frag_atoms)

                    ! store augmented list into sys
                    IF (ALLOCATED(sys%fragments%type_frag(i_group)%fragment_frame(i, j)%frag_atoms)) &
                        DEALLOCATE (sys%fragments%type_frag(i_group)%fragment_frame(i, j)%frag_atoms)
                    ALLOCATE (sys%fragments%type_frag(i_group)%fragment_frame(i, j)%frag_atoms(natom_frag + n_new_wanniers(j, i_group)))
                    sys%fragments%type_frag(i_group)%fragment_frame(i, j)%frag_atoms(1:natom_frag) = &
                        sys%fragments%type_frag(i_group)%fragment(j)%frag_atoms(:)
                    sys%fragments%type_frag(i_group)%fragment_frame(i, j)%frag_atoms(natom_frag+1:) = &
                        new_wannier_centers(1:n_new_wanniers(j, i_group), j, i_group)
                END DO
            END DO
        END DO

    END SUBROUTINE assign_wannier

!******************************************************************************************************************!
!******************************************************************************************************************!

    !> @brief Computes dipole moments for all individual fragments across MD frames.
    !>
    !> This subroutine calculates the dipole moments of individual molecular fragments
    !> for each frame in an MD trajectory. Fragments are organized into
    !> groups (`sys%fragments%type_frag`), each of which may contain multiple fragments
    !> composed of specific atoms.
    !>
    !> @param[inout] sys     -- System data structure (provides fragment definitions,
    !>                          atomic masses, charges, and cell information)
    !> @param[in]    md      -- MD data structure with coordinates for each frame
    !> @param[out]   dipole  -- 4D array (ngroup × nframes × nfrag × 3) containing dipole
    !>                          vectors for each fragment of each group at every frame.
    !>
    SUBROUTINE compute_dipole_frag(dipole, sys, md)

        TYPE(systems), INTENT(INOUT)        :: sys
        TYPE(molecular_dynamics), INTENT(IN) :: md
        TYPE(fragment_group), ALLOCATABLE :: type_frag(:)
        REAL(dp), DIMENSION(:, :, :, :), ALLOCATABLE, INTENT(OUT) :: dipole
        REAL(dp), DIMENSION(:, :, :, :), ALLOCATABLE :: refpoint

        CHARACTER(len=40)                                         :: output_fname
        CHARACTER(LEN=str_len)                                          :: msg
        INTEGER :: m, i, i_group, j, idx, stat, natom_frag, runit
        REAL(dp) :: hmat(3, 3), h_inv(3, 3)
        REAL(dp) :: COM(3), dr(3), r_ref(3), r_unwrapped(3)
        REAL(dp) :: mass_tot

        !!Initialize
        sys%fragments%max_frag = 0
       
        !! Find the maximum number of fragments a group can have
        DO i_group = 1, sys%fragments%ngroup
            sys%fragments%max_frag = MAX(sys%fragments%max_frag, sys%fragments%type_frag(i_group)%nfrag)
        END DO

        !!Allocate
        ALLOCATE (refpoint(sys%fragments%ngroup, sys%framecount, sys%fragments%max_frag, 3))
        ALLOCATE (dipole(sys%fragments%ngroup, sys%framecount, sys%fragments%max_frag, 3))

        !!Initialize
        refpoint = 0.0_dp
        dipole = 0.0_dp

        ! --- lattice (if needed later)
        CALL build_hmat(sys, hmat)
        CALL invert3x3(hmat, h_inv)

        ! --- loop over frames
        DO m = 1, sys%framecount
            ! --- loop over number of fragment sections
            DO i_group = 1, sys%fragments%ngroup
                ! --- loop over the number of atom lists in each fragment section
                DO j = 1, sys%fragments%type_frag(i_group)%nfrag
                    natom_frag = SIZE(sys%fragments%type_frag(i_group)%fragment_frame(m, j)%frag_atoms)

                    ! pick reference atom (first atom in fragment)
                    idx = sys%fragments%type_frag(i_group)%fragment_frame(m, j)%frag_atoms(1)
                    r_ref = md%coord_v(m, idx, :)

                    COM = 0.0_dp
                    mass_tot = 0.0_dp

                    ! loop over all atoms in fragment
                    DO i = 1, natom_frag
                        idx = sys%fragments%type_frag(i_group)%fragment_frame(m, j)%frag_atoms(i)

                        ! skip Wannier centers (X)
                        IF (TRIM(sys%element(idx))=='X') CYCLE

                        ! unwrap atom relative to reference
                        CALL pbc(md%coord_v(m, idx, :), r_ref, sys, dr)
                        r_unwrapped = r_ref + dr

                        ! accumulate mass-weighted sum
                        COM = COM + r_unwrapped*sys%mass_atom(idx)
                        mass_tot = mass_tot + sys%mass_atom(idx)
                    END DO

                    ! finalize COM (for each fragment)
                    IF (mass_tot>0.0_dp) COM = COM/mass_tot

                    ! save reference point for all fragment groups
                    refpoint(i_group, m, j, :) = COM

                    ! Calculate dipoles
                    DO i = 1, natom_frag
                        idx = sys%fragments%type_frag(i_group)%fragment_frame(m, j)%frag_atoms(i)

                        ! unwrap atom relative to reference
                        CALL pbc(md%coord_v(m, idx, :), refpoint(i_group, m, j, :), sys, dr)
                        dipole(i_group, m, j, :) = dipole(i_group, m, j, :) + sys%charge(idx)*dr(:)/bohr2ang

                    END DO

                END DO
            END DO
        END DO
        
        !!Check for jumps in dipole moments and correct if any
        DO i_group = 1, sys%fragments%ngroup
            CALL check_jumps(dipole(i_group, :, 1:sys%fragments%type_frag(i_group)%nfrag, :), sys%fragments%type_frag(i_group)%nfrag, sys)
        END DO

        DEALLOCATE (refpoint)
    END SUBROUTINE compute_dipole_frag

!******************************************************************************************************************!
!******************************************************************************************************************!

    !> @brief Computes the total dipole moment for each frame of a molecular dynamics trajectory.
    !>
    !> This subroutine calculates the time-dependent dipole moment of the system by summing
    !> over all atomic charges multiplied by their displacement vectors (relative to the
    !> center of mass) under periodic boundary conditions.
    !>
    !> @param[inout] sys    --  System data structure (provides atomic and cell information)
    !> @param[in]    md     --  Molecular dynamics data (provides coordinates per frame)
    !> @param[out]   dipole --  3D array (nframes × 1 × 3) containing the dipole vector for each frame
    !>
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

        !Initialize
        dipole = 0.0_dp

        ! --- lattice
        CALL build_hmat(sys, hmat)
        CALL invert3x3(hmat, h_inv)

        ! --- Determine the center of mass of the whole system
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

            ! Compute dipole
            DO i = 1, sys%natom
                CALL pbc(md%coord_v(m, i, :), sys%fragments%refpoint(m, 1, :), sys, dr)
                dipole(m, 1, :) = dipole(m, 1, :) + sys%charge(i)*dr/bohr2ang
            END DO

            ! Convert to Debye
            dipole(m, 1, :) = dipole(m, 1, :)/debye
        END DO

        !!Check for jumps in dipole moments and correct if any
        CALL check_jumps(dipole, sys%mol_num, sys)

        DEALLOCATE (sys%fragments%refpoint)
    END SUBROUTINE compute_dipole

    !> @brief Computes cell parameters from lattice vectors.
    !>
    !> This subroutine calculates the cell parameters of the system from the provided lattice vectors.
    !>
    !> @param[inout] sys    --  System data structure (provides atomic and cell information)
    !>
    SUBROUTINE check_lattice_parameters(sys)
        
        TYPE(systems), INTENT(INOUT) :: sys
        REAL(kind=dp) :: dot_ab, dot_ac, dot_bc
        IF (sys%cell%box_x<0 .OR. sys%cell%box_y<0 .OR.sys%cell%box_z<0 .OR.sys%cell%angle_alpha<0 .OR. sys%cell%angle_beta<0 .OR. sys%cell%angle_gamma<0) THEN
            IF (ALLOCATED(sys%cell%lattice_x) .AND. ALLOCATED(sys%cell%lattice_y) .AND. ALLOCATED(sys%cell%lattice_z)) THEN
                ! Lengths
                sys%cell%box_x = sqrt(sum(sys%cell%lattice_x**2))
                sys%cell%box_y = sqrt(sum(sys%cell%lattice_y**2))
                sys%cell%box_z = sqrt(sum(sys%cell%lattice_z**2))
            
                ! Dot products
                dot_ab = sum(sys%cell%lattice_x * sys%cell%lattice_y)
                dot_ac = sum(sys%cell%lattice_x * sys%cell%lattice_z)
                dot_bc = sum(sys%cell%lattice_y * sys%cell%lattice_z)
            
                ! Angles (convert from radians to degrees)
                sys%cell%angle_alpha = acos(dot_bc / (sys%cell%box_y * sys%cell%box_z)) * 180.0d0 / pi
                sys%cell%angle_beta  = acos(dot_ac / (sys%cell%box_x * sys%cell%box_z)) * 180.0d0 / pi
                sys%cell%angle_gamma = acos(dot_ab / (sys%cell%box_x * sys%cell%box_y)) * 180.0d0 / pi
                
                WRITE (*, '(4X,A, T60, A)') "Found lattice vector and calculating cell parameters:"
                WRITE (*, '(6X,A, T60, F0.6, 1X, F0.6, 1X,F0.6)') "ell parameter: ", sys%cell%box_x, sys%cell%box_y, sys%cell%box_z
                WRITE (*, '(6X,A, T60, F0.6,1X,F0.6,1X,F0.6)') "cell angles: ", sys%cell%angle_alpha, sys%cell%angle_beta,  sys%cell%angle_gamma
            ELSE
                WRITE (error_unit, '(4X,"[ERROR] ",A)') 'Cell parameter not correctly defined' 
                STOP
            END IF
        END IF
        
    END SUBROUTINE check_lattice_parameters

END MODULE dipole_calc
