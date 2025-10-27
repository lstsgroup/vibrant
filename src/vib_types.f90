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

!> @brief Module containing all type definitions and allocation/deallocation procedures.
MODULE vib_types

    USE kinds, ONLY: dp

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: global_settings, systems, molecular_dynamics, static, dipoles, raman, static_property, init_global_settings, &
              init_systems, init_molecular_dynamics, init_static, init_dipoles, init_raman, deallocate_types, fragment_type, fragment_group

    !***************************************************************************
    TYPE spectral_type
        CHARACTER(LEN=40)                               :: read_function ! type of the spectral calculation !
    END TYPE spectral_type

    !***************************************************************************
    TYPE fragment_frame_type
        INTEGER, ALLOCATABLE :: frag_atoms(:)   ! per-frame augmented list (atoms + Wannier)
    END TYPE fragment_frame_type

    !***************************************************************************
    TYPE fragment_type
        INTEGER, ALLOCATABLE :: frag_atoms(:)  ! augmented list (atoms + wannier)
    END TYPE fragment_type
    !***************************************************************************
    TYPE fragments
        LOGICAL                                          ::  frag !yes/no for fragments
        INTEGER                                          :: ngroup !number of fragment sections
        INTEGER                                             :: max_frag !Maximum number of fragments in a fragment group
        INTEGER                                             :: nfrag !Number of fragments in a fragment group 
        INTEGER                                             :: natom_frag !Number of atoms in a fragment
        INTEGER, DIMENSION(:), ALLOCATABLE               :: frag_atoms !Atom list of the fragment
        TYPE(fragment_type), ALLOCATABLE                  :: fragment(:) !(frag)
        TYPE(fragment_group), ALLOCATABLE                 :: type_frag(:) !type of the fragment group (ngroup)
        TYPE(fragment_frame_type), ALLOCATABLE           :: fragment_frame(:, :) ! (frame, frag)
        REAL(kind=dp)                                       :: mass_tot_cell !total mass of the cell
        REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE         :: mass_tot_frag !total mass of the fragment
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE      :: refpoint !ref point for the center of mass
    END TYPE fragments

    !***************************************************************************
    TYPE cell
        CHARACTER(LEN=40)                              :: cell_type
        REAL(kind=dp)                                  :: box_all, box_x, box_y, box_z, vec(3), vec_pbc(3)
        REAL(kind=dp)                                  :: angle_alpha, angle_beta, angle_gamma
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE       :: lattice_x, lattice_y, lattice_z
    END TYPE cell
    !***************************************************************************
    TYPE frame_type
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE        :: frame !frame definition for fragments
    END TYPE frame_type
    !***************************************************************************
    TYPE electric_field
        CHARACTER(LEN=40)                                   :: wannier_xyz !wannier file
        CHARACTER(LEN=40)                                   :: static_dip_file_xyz !dipole file
        INTEGER, DIMENSION(:), ALLOCATABLE                         :: natom_frag_xyz !number of atoms in fragments
        INTEGER, DIMENSION(:, :, :), ALLOCATABLE                     :: fragment_xyz 
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE                :: alpha_xyz !polarizability
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE                :: dip_xyz !dipole
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE                :: alpha_diff_xyz !time derivatives of alphas
    END TYPE electric_field
    !***************************************************************************
    TYPE static_property
        TYPE(atom_information), DIMENSION(:), ALLOCATABLE :: atom !atom information for static calculations
    CONTAINS
        PROCEDURE :: init_staticp
        PROCEDURE :: dealloc_staticp
    END TYPE static_property
    !***************************************************************************
    TYPE resonant_raman
        CHARACTER(LEN=40)                                   :: check_pade !!pade interpolation flag
        INTEGER                                             :: framecount_rtp !rtp framecount
        INTEGER                                             :: framecount_rtp_pade !framecount after Pade interpolation
        REAL(kind=dp)                                       :: dt_rtp !RTP time step in fs
        REAL(kind=dp)                                       :: freq_range_rtp ! RTP freq range
        REAL(kind=dp)                                       :: damping_constant !damping constant 
        TYPE(static_property), DIMENSION(3)  ::  static_dip_rtp ! static field-free dipoles
        TYPE(static_property), DIMENSION(3)  ::  static_dip_x_rtp !static dipoles under electric field in x-direction
        TYPE(static_property), DIMENSION(3)  ::  static_dip_y_rtp !static dipoles under electric field in y-direction
        TYPE(static_property), DIMENSION(3)  ::  static_dip_z_rtp !static dipoles under electric field in z-direction
        TYPE(static_property), DIMENSION(3, 3)  ::  pol_rtp !polarizabilities in time domain
        COMPLEX(kind=dp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: zhat_pol_rtp !Polarizabilities in frequency domain

        !MD RR
        COMPLEX(kind=dp), DIMENSION(:, :), ALLOCATABLE      :: z_iso_resraman !time domain isotropic MD-based polarizabilities
        COMPLEX(kind=dp), DIMENSION(:, :), ALLOCATABLE      :: z_aniso_resraman !time domain anisotropic MD-based polarizabilities
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE      :: dip_x_rtp !dipoles under electric field in x-direction
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE      :: dip_y_rtp !dipoles under electric field in y-direction
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE      :: dip_z_rtp !dipoles under electric field in z-direction
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE      :: alpha_resraman_x_re !real component of the freq-domain x-field pol.
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE      :: alpha_resraman_x_im !imaginary component of the freq-domain x-field pol.
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE      :: alpha_resraman_y_re !real component of the freq-domain y-field pol.
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE      :: alpha_resraman_y_im !imaginary component of the freq-domain y-field pol.
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE      :: alpha_resraman_z_re !real component of the freq-domain z-field pol.
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE      :: alpha_resraman_z_im !imaginary component of the freq-domain z-field pol.
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE      :: alpha_resraman_x_diff_re !time-derivative of alpha_resraman_x_re
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE      :: alpha_resraman_x_diff_im !time-derivative of alpha_resraman_x_im
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE      :: alpha_resraman_y_diff_re ! time-derivative of alpha_resraman_y_re
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE      :: alpha_resraman_y_diff_im ! time-derivative of alpha_resraman_y_im
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE      :: alpha_resraman_z_diff_re !time-derivative of alpha_resraman_z_re
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE      :: alpha_resraman_z_diff_im !time-derivative of alpha_resraman_z_im
    CONTAINS
        PROCEDURE :: init_rr
        PROCEDURE :: init_rr_static_dip   ! initializes static_dip*_rtp
        PROCEDURE :: init_rr_pol          ! initializes pol_rtp (3x3)
        PROCEDURE :: dealloc_rr_all
    END TYPE resonant_raman

    !***************************************************************************
    TYPE global_settings
        LOGICAL                                          ::  md !yes/no
        TYPE(spectral_type)                              ::  spectral_type ! global setting of spectral type 'P' , 'IR' , 'R' etc.
        REAL(kind=dp)                                       ::  temp !temperature
        REAL(kind=dp)                                       ::  fwhm !full width half-maximum
        CHARACTER(LEN=40)                                   :: spectra_verbosity !printing of orthogonal/parallel etc. Raman spectra
    END TYPE global_settings

    !***************************************************************************
    TYPE systems
        INTEGER                                             :: natom               ! number of atoms
        INTEGER                                             :: framecount          ! number of frames
        INTEGER                                             :: mol_num             ! number of molecules
        CHARACTER(LEN=40)                                   :: filename            ! read in file
        CHARACTER(LEN=40)                                   :: type_traj           !  trajectory type for power spec
        CHARACTER(LEN=40)                                   :: input_mass          ! mass weighting (y/n)
        CHARACTER(LEN=2), DIMENSION(:), ALLOCATABLE         :: element              ! ALLOCATE sys%element(sys%natom)
        REAL(kind=dp)                                       :: mass_tot             ! total mass of the system
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE            :: atom_mass_inv_sqrt   ! ALLOCATE sys%atom_mass_inv_sqrt(sys%natom)
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE            :: charge               ! ALLOCATE sys%charge(sys%natom)
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE            :: mass_atom            ! ALLOCATE sys%mass_atom(sys%natom)
        REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE         :: coord                ! ALLOCATE sys%coord(sys%natom, 3)
        REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE         :: mass_mat             ! ALLOCATE sys%mass_mat(sys%natom, sys%natom)
        TYPE(CELL)                                          :: cell        
        TYPE(fragments)                                     :: fragments
    END TYPE systems

    !***************************************************************************
    TYPE molecular_dynamics
        INTEGER                                             :: t_cor               ! correlation depth
        CHARACTER(LEN=40)                                   :: trajectory_file !maybe not needed should be in system type
        CHARACTER(LEN=40)                                   :: velocity_file   !maybe not needed should be in system type
        REAL(kind=dp)                                       :: snapshot_time_step ! snapshots_time_step equal to md%dt ?
        REAL(kind=dp)                                       :: dt   ! MD time step in fs
        REAL(kind=dp)                                       :: freq_range ! maximum frequency range based on dt
        REAL(kind=dp)                                       :: freq_res ! frequency resolution 
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE            :: z       ! autocorrelation function ALLOCATE z(2*md%t_cor-1)
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE      :: v ! coordinates vector ALLOCATE coord_v(sys%framecount, natom, 3)
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE      :: coord_v ! coordinates vector ALLOCATE coord_v(sys%framecount, natom, 3)
        COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE         :: zhat !md%z in frequency domain
        TYPE(fragment_group), ALLOCATABLE                   :: type_frag(:)
    CONTAINS
    END TYPE molecular_dynamics
    !***************************************************************************
    TYPE static
        LOGICAL                                             :: write_mol !yes/no
        INTEGER                                             :: nmodes               ! number of normal modes
        CHARACTER(LEN=40)                                   :: diag_hessian         ! for IR/Raman, yes or no
        CHARACTER(LEN=40)                                   :: normal_freq_file     ! file
        CHARACTER(LEN=40)                                   :: normal_displ_file    ! file
        CHARACTER(LEN=40)                                   :: force_file           ! name of force file
        REAL(kind=dp)                                       :: dx                   ! atom displacement
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE            :: freq                 ! ALLOCATE (stats%freq(stats%nmodes))
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE      :: disp                 ! ALLOCATE (stats%disp(stats%nmodes, sys%natom, 3))
        TYPE(static_property), DIMENSION(:, :), ALLOCATABLE  ::  force
    CONTAINS
        PROCEDURE :: init_force
        PROCEDURE :: dealloc_force
    END TYPE static
    !***************************************************************************
    TYPE dipoles
        CHARACTER(LEN=40)                                   :: dip_file   !file
        CHARACTER(LEN=40)                                   :: dip_x_file !file 
        CHARACTER(LEN=40)                                   :: dip_y_file !file
        CHARACTER(LEN=40)                                   :: dip_z_file !file ! Dipolemoments mybe move to dipole class differenes static_dip_free_file and dip_file
        CHARACTER(LEN=40)                                   :: type_dipole !berry/wannier/dfpt  ! From IR calc what are those`? !!we can add these to static and dipoles
        REAL(kind=dp)                                       :: e_field !electric field
        REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE         :: dip_dq  ! derivatives of the dipoles wrt normal coordinates
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE      :: dipole   ! ALLOCATE static_dip(sys%natom, 3, 2, 3) is this neeeded?
        TYPE(static_property), DIMENSION(3)  ::  static_dip
        TYPE(fragment_group), ALLOCATABLE :: type_frag(:)
    CONTAINS
        PROCEDURE :: init_dip
        PROCEDURE :: dealloc_dip
    END TYPE dipoles
    !***************************************************************************
TYPE fragment_group
    INTEGER :: nfrag = 0
    TYPE(fragment_type), ALLOCATABLE :: fragment(:)
    TYPE(fragment_frame_type), ALLOCATABLE :: fragment_frame(:, :)
END TYPE fragment_group
    !***************************************************************************

    TYPE displacement_type
        TYPE(frame_type), DIMENSION(3)  :: XYZ
    END TYPE displacement_type
    !***************************************************************************

    TYPE atom_information
        TYPE(displacement_type), DIMENSION(2) :: displacement
    END TYPE atom_information
    !***************************************************************************
    TYPE raman
        CHARACTER(LEN=40)                                   :: static_dip_free_file !file
        CHARACTER(LEN=40)                                   :: static_pol_file !file
        CHARACTER(LEN=40)                                   :: wannier_free !file 
        CHARACTER(LEN=40)                                   :: wannier_x !file
        CHARACTER(LEN=40)                                   :: wannier_y !file 
        CHARACTER(LEN=40)                                   :: wannier_z !file
        CHARACTER(LEN=40)                                   :: averaging !isotropic/orientational
        CHARACTER(LEN=40)                                   :: direction !direction for isotropic averaging
        REAL(kind=dp), DIMENSION(:, :, :, :), ALLOCATABLE   :: static_dip_free ! Dipolemoments mybe move to dipole class
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE   :: laser_in
        REAL(kind=dp), DIMENSION(:, :, :, :), ALLOCATABLE   :: static_dip_x !Dipolemoments mybe move to dipole class
        REAL(kind=dp), DIMENSION(:, :, :, :), ALLOCATABLE   :: static_dip_y  ! Dipolemoments mybe move to dipole class
        REAL(kind=dp), DIMENSION(:, :, :, :), ALLOCATABLE   :: static_dip_z  ! Dipolemoments mybe move to dipole class
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE            :: z_iso !isotropic  
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE            :: z_aniso !anisotropic
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE            :: z_ortho !orthogonal
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE            :: z_para !parallel
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE            :: raman_int 
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE      :: pol_dq !pol derivative wrt normal modes  !ALLOCATE (rams%pol_dq(stats%nmodes, 3, 3))
        TYPE(static_property), DIMENSION(3, 3)              :: pol
        TYPE(resonant_raman)                                :: RR
        TYPE(electric_field), DIMENSION(3)                  :: e_field
    CONTAINS
        PROCEDURE :: init_pol
        PROCEDURE :: dealloc_pol
    END TYPE raman


CONTAINS

    !**************************************************************************************!
    !**************************************************************************************!

    !> @brief Initializes the static_property object for a given number of atoms and frames.
    !>
    !> Allocates the displacement arrays for each atom, displacement direction, and
    !> Cartesian component, and sets all frame values to zero.
    !>
    !> @param[inout] this    -- Static property object to initialize
    !> @param[in]    natoms  -- Number of atoms in the system
    !> @param[in]    nframes -- Number of frames (time steps) to allocate per displacement
    !>
    SUBROUTINE init_staticp(this, natoms, nframes)
        CLASS(static_property), INTENT(INOUT) :: this
        INTEGER, INTENT(IN) :: natoms, nframes
        INTEGER :: i, d, xyz

        IF (.NOT. ALLOCATED(this%atom)) ALLOCATE (this%atom(natoms))

        DO i = 1, natoms
            DO d = 1, 2
                DO xyz = 1, 3
                    ALLOCATE (this%atom(i)%displacement(d)%XYZ(xyz)%frame(nframes))
                    this%atom(i)%displacement(d)%XYZ(xyz)%frame = 0.0_dp
                END DO
            END DO
        END DO
    END SUBROUTINE

    !**************************************************************************************!
    !**************************************************************************************!

    !> @brief Deallocates all memory associated with a static_property object.
    !>
    !> Frees displacement arrays for each atom, direction, and Cartesian component,
    !> and deallocates the atom array itself if allocated.
    !>
    !> @param[inout] this -- Static property object to be deallocated
    !>
    SUBROUTINE dealloc_staticp(this)
        CLASS(static_property), INTENT(INOUT) :: this
        INTEGER :: i, d, xyz

        IF (ALLOCATED(this%atom)) THEN
            DO i = 1, SIZE(this%atom)
                DO d = 1, 2
                    DO xyz = 1, 3
                        IF (ALLOCATED(this%atom(i)%displacement(d)%XYZ(xyz)%frame)) THEN
                            DEALLOCATE (this%atom(i)%displacement(d)%XYZ(xyz)%frame)
                        END IF
                    END DO
                END DO
            END DO
            DEALLOCATE (this%atom)
        END IF
    END SUBROUTINE dealloc_staticp

    !**************************************************************************************!
    !**************************************************************************************!

    !> @brief Initializes the force container for all atoms and Cartesian components.
    !>
    !> Allocates the 2D force array (natoms × 3) and initializes each element’s
    !> static_property structure for the given number of atoms and frames.
    !>
    !> @param[inout] this    -- Force object to initialize
    !> @param[in]    natoms  -- Number of atoms in the system
    !> @param[in]    nframes -- Number of frames to allocate
    !> 
    SUBROUTINE init_force(this, natoms, nframes)
        CLASS(static), INTENT(INOUT) :: this
        INTEGER, INTENT(IN) :: natoms, nframes
        INTEGER :: j, k

        ALLOCATE (this%force(natoms, 3))
        DO j = 1, natoms
            DO k = 1, 3
                CALL this%force(j, k)%init_staticp(natoms, nframes)
            END DO
        END DO
    END SUBROUTINE

    !**************************************************************************************!
    !**************************************************************************************!
 
    !> @brief Deallocates all memory associated with the force container.
    !>
    !> Calls `dealloc_staticp` for each atom and Cartesian component,
    !> then deallocates the overall force array if allocated.
    !>
    !> @param[inout] this -- Force object to be deallocated
    !>
    SUBROUTINE dealloc_force(this)
        CLASS(static), INTENT(INOUT) :: this
        INTEGER :: j, k

        IF (ALLOCATED(this%force)) THEN
            DO j = 1, SIZE(this%force, 1)
                DO k = 1, SIZE(this%force, 2)
                    CALL this%force(j, k)%dealloc_staticp()
                END DO
            END DO
            DEALLOCATE (this%force)
        END IF
    END SUBROUTINE

    !**************************************************************************************!
    !**************************************************************************************!

    !> @brief Initializes and allocates the dipole moment containers for all Cartesian components.
    !>
    !> @param[inout] this    -- Dipole object to initialize
    !> @param[in]    natoms  -- Number of atoms in the system
    !> @param[in]    nframes -- Number of frames to allocate
    !>
    SUBROUTINE init_dip(this, natoms, nframes)
        CLASS(dipoles), INTENT(INOUT) :: this
        INTEGER, INTENT(IN) :: natoms, nframes
        INTEGER :: comp

        DO comp = 1, 3
            CALL this%static_dip(comp)%init_staticp(natoms, nframes)
        END DO
    END SUBROUTINE

    !**************************************************************************************!
    !**************************************************************************************!

    !> @brief Deallocates memory for all Cartesian dipole components.
    !>
    !> @param[inout] this -- Dipole object to be deallocated
    !>
    SUBROUTINE dealloc_dip(this)
        CLASS(dipoles), INTENT(INOUT) :: this
        INTEGER :: c
        DO c = 1, 3
            CALL this%static_dip(c)%dealloc_staticp()
        END DO
    END SUBROUTINE

    !**************************************************************************************!
    !**************************************************************************************!

    !> @brief Initializes the polarizability tensor components for all atoms and frames.
    !>
    !> @param[inout] this    -- Raman (polarizability) object to initialize
    !> @param[in]    natoms  -- Number of atoms in the system
    !> @param[in]    nframes -- Number of frames to allocate
    !>
    SUBROUTINE init_pol(this, natoms, nframes)
        CLASS(raman), INTENT(INOUT) :: this
        INTEGER, INTENT(IN) :: natoms, nframes
        INTEGER :: i, j

        DO i = 1, 3
            DO j = 1, 3
                CALL this%pol(i, j)%init_staticp(natoms, nframes)
            END DO
        END DO
    END SUBROUTINE

    !**************************************************************************************!
    !**************************************************************************************!

    !> @brief Initializes the resonance Raman (RR) control parameters to default values.
    !>
    !> @param[inout] this -- RR object to initialize
    !>
    SUBROUTINE init_rr(this)
        CLASS(resonant_raman), INTENT(INOUT) :: this

        this%check_pade = ''
        this%framecount_rtp = -1
        this%framecount_rtp_pade = -1
        this%dt_rtp = -1.0_dp
        this%freq_range_rtp = -1.0_dp
        this%damping_constant = -1.0_dp
    END SUBROUTINE

    !**************************************************************************************!
    !**************************************************************************************!

    !> @brief Deallocates all elements of the 3×3 polarizability tensor.
    !>
    !> @param[inout] this -- Raman (polarizability) object to be deallocated
    SUBROUTINE dealloc_pol(this)
        CLASS(raman), INTENT(INOUT) :: this
        INTEGER :: i, j
        DO i = 1, 3
            DO j = 1, 3
                CALL this%pol(i, j)%dealloc_staticp()
            END DO
        END DO
    END SUBROUTINE

    !**************************************************************************************!
    !**************************************************************************************!

    !> @brief Initializes an array of static dipole properties for resonance Raman (RR) calculations.
    !>
    !> @param[inout] this    -- RR object associated with the dipole array
    !> @param[inout] arr     -- Array of static_property objects to initialize
    !> @param[in]    natoms  -- Number of atoms in the system
    !> @param[in]    nframes -- Number of frames to allocate
    !>
    SUBROUTINE init_rr_static_dip(this, arr, natoms, nframes)
        CLASS(resonant_raman), INTENT(INOUT) :: this
        TYPE(static_property), DIMENSION(:), INTENT(INOUT) :: arr
        INTEGER, INTENT(IN)   :: natoms, nframes
        INTEGER :: c
        DO c = 1, SIZE(arr)
            CALL arr(c)%init_staticp(natoms, nframes)
        END DO
    END SUBROUTINE

    !**************************************************************************************!
    !**************************************************************************************!

    !> @brief Initializes and allocates the polarizability tensor components for resonance 
    !>        Raman (RR) calculations.
    !> @param[inout] this    -- RR object to initialize
    !> @param[in]    natoms  -- Number of atoms in the system
    !> @param[in]    nframes -- Number of frames to allocate
    !>
    SUBROUTINE init_rr_pol(this, natoms, nframes)
        CLASS(resonant_raman), INTENT(INOUT) :: this
        INTEGER, INTENT(IN)    :: natoms, nframes
        INTEGER :: i, j

        DO i = 1, 3
            DO j = 1, 3
                CALL this%pol_rtp(i, j)%init_staticp(natoms, nframes)
            END DO
        END DO
    END SUBROUTINE init_rr_pol

    !**************************************************************************************!
    !**************************************************************************************!

    !> @brief Deallocates all static dipole and polarizability components in the resonance
    !>        Raman (RR) object.
    !> @param[inout] this -- RR object to be deallocated
    !>
    SUBROUTINE dealloc_rr_all(this)
        CLASS(resonant_raman), INTENT(INOUT) :: this
        INTEGER :: c, i, j
        DO c = 1, 3
            CALL this%static_dip_rtp(c)%dealloc_staticp()
            CALL this%static_dip_x_rtp(c)%dealloc_staticp()
            CALL this%static_dip_y_rtp(c)%dealloc_staticp()
            CALL this%static_dip_z_rtp(c)%dealloc_staticp()
        END DO
        DO i = 1, 3
            DO j = 1, 3
                CALL this%pol_rtp(i, j)%dealloc_staticp()
            END DO
        END DO
    END SUBROUTINE dealloc_rr_all

    !**************************************************************************************!
    !**************************************************************************************!

    !> @brief Initializes global simulation settings with default values.
    !>
    !> @param[out] gs  -- Global settings object to initialize
    !>
    SUBROUTINE init_global_settings(gs)
        TYPE(global_settings), INTENT(out) :: gs
        gs%spectra_verbosity = 'normal'
        gs%temp = -1.0_dp
        gs%spectral_type%read_function = ''
    END SUBROUTINE init_global_settings

    !**************************************************************************************!
    !**************************************************************************************!

    !> @brief Initializes the system structure with default values.
    !>
    !> @param[out] sys -- System structure to initialize
    !>
    SUBROUTINE init_systems(sys)
        TYPE(systems), INTENT(out) :: sys

        sys%natom = -1
        sys%framecount = -1
        sys%mol_num = -1
        sys%filename = ''
        sys%type_traj = ''
        sys%input_mass = ''
        sys%mass_tot = -1.0_dp

        !cell
        sys%cell%box_x = -1.0_dp
        sys%cell%box_y = -1.0_dp
        sys%cell%box_z = -1.0_dp

        sys%cell%angle_alpha    = -1.0_dp
        sys%cell%angle_beta     = -1.0_dp
        sys%cell%angle_gamma    = -1.0_dp
    END SUBROUTINE init_systems

    !**************************************************************************************!
    !**************************************************************************************!
 
    !> @brief Initializes molecular dynamics control parameters with default values.
    !>
    !> @param[out] md  -- MD control structure to initialize
    !>
    SUBROUTINE init_molecular_dynamics(md)
        TYPE(molecular_dynamics), INTENT(out) :: md
        md%t_cor = -1
        md%trajectory_file = ''
        md%velocity_file = ''
        md%snapshot_time_step = -1.0_dp
        md%dt = -1.0_dp
        md%freq_range = -1.0_dp
        md%freq_res = -1.0_dp
    END SUBROUTINE init_molecular_dynamics

    !**************************************************************************************!
    !**************************************************************************************!

    !> @brief Initializes the static calculation settings to default values.
    !>
    !> @param[out] stats  -- Static calculation data structure to initialize
    !>
    SUBROUTINE init_static(stats)
        TYPE(static), INTENT(out) :: stats
        stats%nmodes = -1
        stats%diag_hessian = ''
        stats%normal_freq_file = ''
        stats%normal_displ_file = ''
        stats%force_file = ''
        stats%dx = -1.0_dp
    END SUBROUTINE init_static

    !**************************************************************************************!
    !**************************************************************************************!

    !> @brief Initializes the dipole configuration with default values.
    !>
    !> @param[out] dip  -- Dipole data structure to initialize
    !>
    SUBROUTINE init_dipoles(dip)
        TYPE(dipoles), INTENT(out) :: dip
        dip%type_dipole = ''
        dip%dip_file = ''
    END SUBROUTINE init_dipoles

    !**************************************************************************************!
    !**************************************************************************************!

    !> @brief Initializes the Raman calculation structure with default values.
    !>
    !> @param[out] ram  Raman data structure to initialize
    !> 
    SUBROUTINE init_raman(ram)
        TYPE(raman), INTENT(out) :: ram
        !numeric
        ram%static_dip_free_file = ''
        ram%static_pol_file = ''
        !analytic
        ram%wannier_free = ''
        ram%wannier_x = ''
        ram%wannier_y = ''
        ram%wannier_z = ''
        ram%averaging = ''
        ram%direction = ''

        CALL ram%RR%init_rr()
    END SUBROUTINE init_raman

    !**************************************************************************************!
    !**************************************************************************************!

    !> @brief Deallocates all major data structures used in the simulation framework.
    !>
    !> @param[inout, optional] gs    -- Global settings structure to deallocate  
    !> @param[inout, optional] sys   -- System data structure to deallocate  
    !> @param[inout, optional] md    -- Molecular dynamics structure to deallocate  
    !> @param[inout, optional] stats -- Static calculation structure to deallocate  
    !> @param[inout, optional] ram   -- Raman structure to deallocate  
    !> @param[inout, optional] dip   -- Dipole structure to deallocate
    !>
    SUBROUTINE deallocate_types(gs, sys, md, stats, ram, dip)

        TYPE(global_settings), INTENT(INOUT), OPTIONAL :: gs
        TYPE(systems), INTENT(INOUT), OPTIONAL :: sys
        TYPE(molecular_dynamics), INTENT(INOUT), OPTIONAL:: md
        TYPE(static), INTENT(INOUT), OPTIONAL:: stats
        TYPE(raman), INTENT(INOUT), OPTIONAL:: ram
        TYPE(dipoles), INTENT(INOUT), OPTIONAL:: dip

        ! systems
        IF (PRESENT(sys)) THEN
            CALL deallocate_system(sys)
        END IF
        ! molecular_dynamics
        IF (PRESENT(md)) THEN
            CALL deallocate_md(md)
        END IF
        ! static
        IF (PRESENT(stats)) THEN
            CALL deallocate_stats(stats)
        END IF
        !raman
        IF (PRESENT(ram)) THEN
            CALL deallocate_raman(ram)
        END IF
        !dipoles
        IF (PRESENT(dip)) THEN
            CALL deallocate_dipoles(dip)
        END IF
    END SUBROUTINE deallocate_types

    !**************************************************************************************!
    !**************************************************************************************!

    !> @brief Deallocates all memory associated with the system structure.
    !>
    !> @param[inout] sys -- System data structure to deallocate
    !>
     SUBROUTINE deallocate_system(sys)
        TYPE(systems), INTENT(INOUT) :: sys

        IF (ALLOCATED(sys%element)) DEALLOCATE (sys%element)
        IF (ALLOCATED(sys%atom_mass_inv_sqrt)) DEALLOCATE (sys%atom_mass_inv_sqrt)
        IF (ALLOCATED(sys%charge)) DEALLOCATE (sys%charge)
        IF (ALLOCATED(sys%mass_atom)) DEALLOCATE (sys%mass_atom)
        IF (ALLOCATED(sys%coord)) DEALLOCATE (sys%coord)
        IF (ALLOCATED(sys%mass_mat)) DEALLOCATE (sys%mass_mat)

        IF (ALLOCATED(sys%fragments%fragment)) DEALLOCATE (sys%fragments%fragment)
        IF (ALLOCATED(sys%fragments%mass_tot_frag)) DEALLOCATE (sys%fragments%mass_tot_frag)
        IF (ALLOCATED(sys%fragments%refpoint)) DEALLOCATE (sys%fragments%refpoint)

    END SUBROUTINE deallocate_system

    !**************************************************************************************!
    !**************************************************************************************!

    !> @brief Deallocates arrays associated with molecular dynamics data.
    !>
    !>
    !> @param[inout] md -- MD structure to deallocate
    !>
    SUBROUTINE deallocate_md(md)
        TYPE(molecular_dynamics), INTENT(INOUT) :: md

        IF (ALLOCATED(md%z)) DEALLOCATE (md%z)
        IF (ALLOCATED(md%v)) DEALLOCATE (md%v)
        IF (ALLOCATED(md%coord_v)) DEALLOCATE (md%coord_v)
        IF (ALLOCATED(md%zhat)) DEALLOCATE (md%zhat)

    END SUBROUTINE deallocate_md

    !**************************************************************************************!
    !**************************************************************************************!

    !> @brief Deallocates memory used in static calculations.
    !>
    !> @param[inout] stats -- Static calculation structure to deallocate
    SUBROUTINE deallocate_stats(stats)
        TYPE(static), INTENT(INOUT):: stats

        INTEGER                     :: xyz, i

        IF (ALLOCATED(stats%freq)) DEALLOCATE (stats%freq)
        IF (ALLOCATED(stats%disp)) DEALLOCATE (stats%disp)
        ! Deallocating Static Property
        CALL stats%dealloc_force()
    END SUBROUTINE deallocate_stats

    !**************************************************************************************!
    !**************************************************************************************!

    !> @brief Deallocates all memory associated with the dipole data structure.
    !>
    !> @param[inout] dips -- Dipole structure to deallocate
    !>
    SUBROUTINE deallocate_dipoles(dips)
        TYPE(dipoles), INTENT(INOUT) :: dips
        !IF (ALLOCATED(dips%dip_file)) DEALLOCATE(dips%dip_file)
        IF (ALLOCATED(dips%dip_dq)) DEALLOCATE (dips%dip_dq)
        IF (ALLOCATED(dips%dipole)) DEALLOCATE (dips%dipole)
        CALL dips%dealloc_dip()

    END SUBROUTINE deallocate_dipoles

    !**************************************************************************************!
    !**************************************************************************************!

    !> @brief Deallocates all memory associated with the Raman data structure.
    !>
    !> @param[inout] rams -- Raman structure to deallocate
    !>
    SUBROUTINE deallocate_raman(rams)
        TYPE(raman), INTENT(INOUT):: rams

        INTEGER                     :: xyz

        IF (ALLOCATED(rams%laser_in)) DEALLOCATE (rams%laser_in)
        IF (ALLOCATED(rams%static_dip_free)) DEALLOCATE (rams%static_dip_free)
        DO xyz = 1, 3
            !IF (ALLOCATED(rams%e_field(xyz)%static_dip_xyz)) DEALLOCATE(rams%e_field(xyz)%static_dip_xyz)
            IF (ALLOCATED(rams%e_field(xyz)%alpha_xyz)) DEALLOCATE (rams%e_field(xyz)%alpha_xyz)
            IF (ALLOCATED(rams%e_field(xyz)%dip_xyz)) DEALLOCATE (rams%e_field(xyz)%dip_xyz)
            IF (ALLOCATED(rams%e_field(xyz)%alpha_diff_xyz)) DEALLOCATE (rams%e_field(xyz)%alpha_diff_xyz)
            IF (ALLOCATED(rams%e_field(xyz)%fragment_xyz)) DEALLOCATE (rams%e_field(xyz)%fragment_xyz)
            IF (ALLOCATED(rams%e_field(xyz)%natom_frag_xyz)) DEALLOCATE (rams%e_field(xyz)%natom_frag_xyz)
        END DO

        CALL rams%dealloc_pol()

        !IF (ALLOCATED(rams%static_pol_file)) DEALLOCATE(rams%static_pol_file)
        IF (ALLOCATED(rams%z_iso)) DEALLOCATE (rams%z_iso)
        IF (ALLOCATED(rams%z_aniso)) DEALLOCATE (rams%z_aniso)
        IF (ALLOCATED(rams%z_ortho)) DEALLOCATE (rams%z_ortho)
        IF (ALLOCATED(rams%z_para)) DEALLOCATE (rams%z_para)
        IF (ALLOCATED(rams%raman_int)) DEALLOCATE (rams%raman_int)
        IF (ALLOCATED(rams%pol_dq)) DEALLOCATE (rams%pol_dq)


        CALL rams%rr%dealloc_rr_all()
        IF (ALLOCATED(rams%RR%zhat_pol_rtp)) DEALLOCATE (rams%RR%zhat_pol_rtp)
        IF (ALLOCATED(rams%RR%z_iso_resraman)) DEALLOCATE (rams%RR%z_iso_resraman)
        IF (ALLOCATED(rams%RR%z_aniso_resraman)) DEALLOCATE (rams%RR%z_aniso_resraman)
        IF (ALLOCATED(rams%RR%dip_x_rtp)) DEALLOCATE (rams%RR%dip_x_rtp)
        IF (ALLOCATED(rams%RR%dip_y_rtp)) DEALLOCATE (rams%RR%dip_y_rtp)
        IF (ALLOCATED(rams%RR%dip_z_rtp)) DEALLOCATE (rams%RR%dip_z_rtp)
        IF (ALLOCATED(rams%RR%alpha_resraman_x_re)) DEALLOCATE (rams%RR%alpha_resraman_x_re)
        IF (ALLOCATED(rams%RR%alpha_resraman_x_im)) DEALLOCATE (rams%RR%alpha_resraman_x_im)
        IF (ALLOCATED(rams%RR%alpha_resraman_y_re)) DEALLOCATE (rams%RR%alpha_resraman_y_re)
        IF (ALLOCATED(rams%RR%alpha_resraman_y_im)) DEALLOCATE (rams%RR%alpha_resraman_y_im)
        IF (ALLOCATED(rams%RR%alpha_resraman_z_re)) DEALLOCATE (rams%RR%alpha_resraman_z_re)
        IF (ALLOCATED(rams%RR%alpha_resraman_z_im)) DEALLOCATE (rams%RR%alpha_resraman_z_im)
        IF (ALLOCATED(rams%RR%alpha_resraman_x_diff_re)) DEALLOCATE (rams%RR%alpha_resraman_x_diff_re)
        IF (ALLOCATED(rams%RR%alpha_resraman_x_diff_im)) DEALLOCATE (rams%RR%alpha_resraman_x_diff_im)
        IF (ALLOCATED(rams%RR%alpha_resraman_y_diff_re)) DEALLOCATE (rams%RR%alpha_resraman_y_diff_re)
        IF (ALLOCATED(rams%RR%alpha_resraman_y_diff_im)) DEALLOCATE (rams%RR%alpha_resraman_y_diff_im)
        IF (ALLOCATED(rams%RR%alpha_resraman_z_diff_re)) DEALLOCATE (rams%RR%alpha_resraman_z_diff_re)
        IF (ALLOCATED(rams%RR%alpha_resraman_z_diff_im)) DEALLOCATE (rams%RR%alpha_resraman_z_diff_im)

    END SUBROUTINE deallocate_raman

END MODULE vib_types
