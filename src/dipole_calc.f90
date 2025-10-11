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
    USE constants, ONLY: debye, bohr2ang
    USE ISO_FORTRAN_ENV, ONLY: output_unit, error_unit
    USE vib_types, ONLY: global_settings, systems, molecular_dynamics, static, dipoles, fragment_group
    USE cell_types, ONLY: build_hmat, pbc, invert3x3, determinant3x3
    USE output_io, ONLY: check_file_open
    USE OMP_LIB

    IMPLICIT NONE

    PUBLIC :: compute_dipole, check_jumps, assign_wannier, compute_dipole_frag
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
        
        !Cut-off distance between the atom and the Wannier center  
        dist = 1.1_dp 
        ! --- lattice
        CALL build_hmat(sys, hmat)
        CALL invert3x3(hmat, h_inv)

        !!Allocate
        DO i_group = 1, sys%fragments%ngroup
            frag_count = sys%fragments%type_frag(i_group)%nfrag
            ! --- allocate storage for augmented fragment lists
            IF (.NOT. ALLOCATED(sys%fragments%type_frag(i_group)%fragment_frame)) THEN
                ALLOCATE (sys%fragments%type_frag(i_group)%fragment_frame(sys%framecount, frag_count))
            END IF
        END DO

        ! --- loop over frames
        DO i = 1, sys%framecount
            ! --- loop over number of fragment sections
            DO i_group = 1, sys%fragments%ngroup
                frag_count = sys%fragments%type_frag(i_group)%nfrag
                ! --- loop over the number of atom lists in each fragment section
                DO j = 1, frag_count
                    natom_frag = SIZE(sys%fragments%type_frag(i_group)%fragment(j)%frag_atoms)

                    ! start from static fragment definition
                    frag_atoms_frame = sys%fragments%type_frag(i_group)%fragment(j)%frag_atoms
                    ! --- loop over the number of atoms in each atom_list
                    DO m = 1, natom_frag
                        outer: DO i_atom = 1, sys%natom
                            IF (sys%element(i_atom)/='X') CYCLE   ! only Wannier centers

                            CALL pbc(md%coord_v(i, i_atom, :), &
                                     md%coord_v(i, sys%fragments%type_frag(i_group)%fragment(j)%frag_atoms(m), :), &
                                     sys, dr)

                            IF (SQRT(DOT_PRODUCT(dr, dr))<dist) THEN
                                ! skip if this Wannier is already in the fragment
                                IF (ANY(i_atom==frag_atoms_frame)) CYCLE outer

                                ! check if Wannier is too close to any atom outside the fragment
                                inner: DO l = 1, sys%natom
                                    IF (sys%element(l)=='X') CYCLE inner         ! skip other Wanniers
                                    IF (ANY(l==sys%fragments%type_frag(i_group)%fragment(j)%frag_atoms)) CYCLE inner  ! skip current fragment atoms

                                    CALL pbc(md%coord_v(i, l, :), md%coord_v(i, i_atom, :), sys, dr2)
                                    IF (SQRT(DOT_PRODUCT(dr2, dr2))<dist) CYCLE outer
                                END DO inner

                                ! append if it survived the filters
                                frag_atoms_frame = [frag_atoms_frame, i_atom]
                            END IF
                        END DO outer
                    END DO

                    ! store augmented list into sys
                    IF (ALLOCATED(sys%fragments%type_frag(i_group)%fragment_frame(i, j)%frag_atoms)) &
                        DEALLOCATE (sys%fragments%type_frag(i_group)%fragment_frame(i, j)%frag_atoms)
                    ALLOCATE (sys%fragments%type_frag(i_group)%fragment_frame(i, j)%frag_atoms(SIZE(frag_atoms_frame)))
                    sys%fragments%type_frag(i_group)%fragment_frame(i, j)%frag_atoms = frag_atoms_frame

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

END MODULE dipole_calc
