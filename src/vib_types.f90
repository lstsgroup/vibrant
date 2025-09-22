MODULE vib_types

    USE kinds, ONLY: dp

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: global_settings, systems, molecular_dynamics, static, dipoles, raman, static_property, init_global_settings, &
              init_systems, init_molecular_dynamics, init_static, init_dipoles, init_raman, deallocate_types

    !***************************************************************************
    TYPE spectral_type
        CHARACTER(LEN=40)                               :: read_function
    END TYPE spectral_type

    !***************************************************************************
    TYPE fragments
        INTEGER                                             :: nfrag
        INTEGER, DIMENSION(:), ALLOCATABLE                  :: natom_frag
        INTEGER, DIMENSION(:, :, :), ALLOCATABLE            :: fragment
        REAL(kind=dp)                                       :: mass_tot_cell
        REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE         :: mass_tot_frag
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE      :: refpoint
    END TYPE fragments

    !***************************************************************************
    TYPE cell
        CHARACTER(LEN=40)                              :: cell_type
        REAL(kind=dp)                                  :: box_all, box_x, box_y, box_z, vec(3), vec_pbc(3)
    END TYPE cell
    !***************************************************************************
    TYPE frame_type
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: frame
    END TYPE frame_type
    !***************************************************************************
    TYPE electric_field
        CHARACTER(LEN=40)                                   :: wannier_xyz
        CHARACTER(LEN=40)                                   :: static_dip_file_xyz

        INTEGER, DIMENSION(:), ALLOCATABLE                         :: natom_frag_xyz
        INTEGER, DIMENSION(:, :, :), ALLOCATABLE                     :: fragment_xyz
        !REAL(kind=dp), DIMENSION(:, :, :, :), ALLOCATABLE   :: static_dip_xyz
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE                :: alpha_xyz
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE                :: dip_xyz
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE                :: alpha_diff_xyz
    END TYPE electric_field
    !***************************************************************************
    TYPE static_property
        TYPE(atom_information), DIMENSION(:), ALLOCATABLE :: atom
    CONTAINS
        PROCEDURE :: init_staticp
        PROCEDURE :: dealloc_staticp
    END TYPE static_property
    !***************************************************************************
    TYPE resonant_raman
        CHARACTER(LEN=40)                                   :: check_pade
        INTEGER                                             :: framecount_rtp
        INTEGER                                             :: framecount_rtp_pade
        REAL(kind=dp)                                       :: dt_rtp
        REAL(kind=dp)                                       :: freq_range_rtp
        REAL(kind=dp)                                       :: damping_constant
        TYPE(static_property), DIMENSION(3)  ::  static_dip_rtp
        TYPE(static_property), DIMENSION(3)  ::  static_dip_x_rtp
        TYPE(static_property), DIMENSION(3)  ::  static_dip_y_rtp
        TYPE(static_property), DIMENSION(3)  ::  static_dip_z_rtp
        TYPE(static_property), DIMENSION(3, 3)  ::  pol_rtp
        COMPLEX(kind=dp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: zhat_pol_rtp
        !CHARACTER(LEN=40)                                   :: rtp_dipole_x, rtp_dipole_y, rtp_dipole_z
        !COMPLEX(kind=dp), DIMENSION(:, :), ALLOCATABLE      :: z_iso_resraman, z_aniso_resraman
    CONTAINS
        PROCEDURE :: init_rr_static_dip   ! initializes static_dip*_rtp
        PROCEDURE :: init_rr_pol          ! initializes pol_rtp (3x3)
        PROCEDURE :: dealloc_rr_all
    END TYPE resonant_raman

    !***************************************************************************
    TYPE global_settings
        LOGICAL                                          ::  md !yes/no
        TYPE(spectral_type)                              ::  spectral_type ! global setting of spectral type 'P' , 'IR' , 'R' etc.
        REAL(kind=dp)                                       ::  temp
        REAL(kind=dp)                                       ::  fwhm
    END TYPE global_settings

    !***************************************************************************
    TYPE systems
        INTEGER                                             :: natom               ! number of atoms
        INTEGER                                             :: framecount          ! number of frames
        INTEGER                                             :: mol_num             ! number of moleces ?
        CHARACTER(LEN=40)                                   :: filename            ! read in file
        CHARACTER(LEN=40)                                   :: frag_type           !  ???
        CHARACTER(LEN=40)                                   :: type_traj           !  trajectory type for power spec
        CHARACTER(LEN=40)                                   :: input_mass          ! mass weighting (y/n)
        CHARACTER(LEN=40)                                   :: system              ! fragment approach (1) or molecular approach? (2)
        CHARACTER(LEN=40)                                   :: periodic !          ! contain more than one molecule? (y/n)
        CHARACTER(LEN=2), DIMENSION(:), ALLOCATABLE         :: element              ! ALLOCATE sys%element(sys%natom)
        REAL(kind=dp)                                       :: mass_tot
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE            :: atom_mass_inv_sqrt   ! ALLOCATE sys%atom_mass_inv_sqrt(sys%natom)
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE            :: charge               ! ALLOCATE sys%charge(sys%natom)
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE            :: mass_atom            ! ALLOCATE sys%mass_atom(sys%natom)
        REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE         :: coord                ! ALLOCATE sys%coord(sys%natom, 3)
        REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE         :: mass_mat             ! ALLOCATE sys%mass_mat(sys%natom, sys%natom)
        TYPE(CELL)                                          :: cell
        TYPE(fragments)                                     :: fragments! <--- NEEDED?
    END TYPE systems

    !***************************************************************************
    TYPE molecular_dynamics
        INTEGER                                             :: t_cor               ! correlation depth
        CHARACTER(LEN=40)                                   :: trajectory_file !maybe not needed should be in system type
        CHARACTER(LEN=40)                                   :: velocity_file   !maybe not needed should be in system type
        REAL(kind=dp)                                       :: snapshot_time_step ! snapshots_time_step equal to dt ?
        REAL(kind=dp)                                       :: dt   ! not quite sure ?
        REAL(kind=dp)                                       :: freq_range ! not quite sure ?
        REAL(kind=dp)                                       :: freq_res ! not sure if right here ?
        REAL(kind=dp)                                       :: sinc_const !<---- CONSTANT ? Is this needed here?
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE            :: z       ! correlations vector ?
        REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE         :: velos_v ! velocity vector
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE      :: v ! coordinates vector
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE      :: coord_v ! coordinates vector ALLOCATE coord_v(sys%framecount, natom, 3)
        COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE         :: zhat ! correlations vector hat ?
    END TYPE molecular_dynamics
    !***************************************************************************
    TYPE static
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
        CHARACTER(LEN=40)                                   :: dip_file      !
        CHARACTER(LEN=40)                                   :: dip_x_file
        CHARACTER(LEN=40)                                   :: dip_y_file
        CHARACTER(LEN=40)                                   :: dip_z_file        ! Dipolemoments mybe move to dipole class differenes static_dip_free_file and dip_file
        CHARACTER(LEN=40)                                   :: type_dipole ! From IR calc what are those`? !!we can add these to static and dipoles
        REAL(kind=dp)                                       :: e_field !electric field
        REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE         :: dip_dq               !
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE      :: dipole               ! ALLOCATE static_dip(sys%natom, 3, 2, 3) is this neeeded?
        TYPE(static_property), DIMENSION(3)  ::  static_dip
    CONTAINS
        PROCEDURE :: init_dip
        PROCEDURE :: dealloc_dip
        !  REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE      :: dip
        !LOGICAL                                            ::  fragment !<yes/no>
    END TYPE dipoles
    !***************************************************************************

    TYPE displacement_type
        TYPE(frame_type), DIMENSION(3)  :: XYZ
    END TYPE displacement_type
    !***************************************************************************

    TYPE atom_information
        TYPE(displacement_type), DIMENSION(2) :: displacement
        !REAL(kind=dp), DIMENSION(3, 3), ALLOCATABLE :: pol
    END TYPE atom_information
    !***************************************************************************
    TYPE raman
        !LOGICAL                                             :: polarizability_type ! <numeric/analytic>
        !numeric
        CHARACTER(LEN=40)                                   :: static_dip_free_file
        REAL(kind=dp), DIMENSION(:, :, :, :), ALLOCATABLE   :: static_dip_free          ! Dipolemoments mybe move to dipole class
        REAL(kind=dp)                                       :: laser_in
        REAL(kind=dp), DIMENSION(:, :, :, :), ALLOCATABLE   :: static_dip_x             ! Dipolemoments mybe move to dipole class
        REAL(kind=dp), DIMENSION(:, :, :, :), ALLOCATABLE   :: static_dip_y             ! Dipolemoments mybe move to dipole class
        REAL(kind=dp), DIMENSION(:, :, :, :), ALLOCATABLE   :: static_dip_z             ! Dipolemoments mybe move to dipole class
        !end numeric
        ! analytic
        CHARACTER(LEN=40)                                   :: static_pol_file
        !end analytic
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE            :: z_iso, z_aniso, z_ortho, z_para
        CHARACTER(LEN=40)                                   :: wannier_free, wannier_x, wannier_y, wannier_z ! same as static_dip_free_file
        CHARACTER(LEN=40)                                   :: averaging, direction
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE            :: raman_int
        !REAL(kind=dp), DIMENSION(:, :, :, :, :), ALLOCATABLE:: pol ! ALLOCATE rams%pol(sys%natom, 3, 2, 3, 3)
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE      :: pol_dq !ALLOCATE (rams%pol_dq(stats%nmodes, 3, 3))
        TYPE(static_property), DIMENSION(3, 3)  :: pol
        TYPE(resonant_raman)                                :: RR
        TYPE(electric_field), DIMENSION(3)                  :: e_field
    CONTAINS
        PROCEDURE :: init_pol
        PROCEDURE :: dealloc_pol
    END TYPE raman

    !***************************************************************************

CONTAINS
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

    SUBROUTINE init_dip(this, natoms, nframes)
        CLASS(dipoles), INTENT(INOUT) :: this
        INTEGER, INTENT(IN) :: natoms, nframes
        INTEGER :: comp

        DO comp = 1, 3
            CALL this%static_dip(comp)%init_staticp(natoms, nframes)
        END DO
    END SUBROUTINE

    SUBROUTINE dealloc_dip(this)
        CLASS(dipoles), INTENT(INOUT) :: this
        INTEGER :: c
        DO c = 1, 3
            CALL this%static_dip(c)%dealloc_staticp()
        END DO
    END SUBROUTINE

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
    SUBROUTINE dealloc_pol(this)
        CLASS(raman), INTENT(INOUT) :: this
        INTEGER :: i, j
        DO i = 1, 3
            DO j = 1, 3
                CALL this%pol(i, j)%dealloc_staticp()
            END DO
        END DO
    END SUBROUTINE

    SUBROUTINE init_rr_static_dip(this, arr, natoms, nframes)
        CLASS(resonant_raman), INTENT(INOUT) :: this
        TYPE(static_property), DIMENSION(:), INTENT(INOUT) :: arr
        INTEGER, INTENT(IN)   :: natoms, nframes
        INTEGER :: c
        DO c = 1, SIZE(arr)
            CALL arr(c)%init_staticp(natoms, nframes)
        END DO
    END SUBROUTINE

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

    SUBROUTINE init_global_settings(gs)
        TYPE(global_settings), INTENT(out) :: gs

        gs%temp = -1.0_dp
        gs%spectral_type%read_function = ''
    END SUBROUTINE init_global_settings

    SUBROUTINE init_systems(sys)
        TYPE(systems), INTENT(out) :: sys

        sys%natom = -1
        sys%framecount = -1
        sys%mol_num = -1
        sys%filename = ''
        sys%frag_type = ''
        sys%type_traj = ''
        sys%input_mass = ''
        sys%system = ''
        sys%periodic = ''
        sys%mass_tot = -1.0_dp
    END SUBROUTINE init_systems

    SUBROUTINE init_molecular_dynamics(md)
        TYPE(molecular_dynamics), INTENT(out) :: md
        md%t_cor = -1
        md%trajectory_file = ''
        md%velocity_file = ''
        md%snapshot_time_step = -1.0_dp
        md%dt = -1.0_dp
        md%freq_range = -1.0_dp
        md%freq_res = -1.0_dp
        md%sinc_const = -1.0_dp
    END SUBROUTINE init_molecular_dynamics

    SUBROUTINE init_static(stats)
        TYPE(static), INTENT(out) :: stats
        stats%nmodes = -1
        stats%diag_hessian = ''
        stats%normal_freq_file = ''
        stats%normal_displ_file = ''
        stats%force_file = ''
        stats%dx = -1.0_dp
    END SUBROUTINE init_static

    SUBROUTINE init_dipoles(dip)
        TYPE(dipoles), INTENT(out) :: dip
        dip%type_dipole = ''
        dip%dip_file = ''
    END SUBROUTINE init_dipoles

    SUBROUTINE init_raman(ram)
        TYPE(raman), INTENT(out) :: ram
        ram%static_pol_file = ''
        ram%laser_in = -1.0_dp
    END SUBROUTINE init_raman

    SUBROUTINE deallocate_types(gs, sys, md, stats, ram, dip)

        TYPE(global_settings), INTENT(INOUT), OPTIONAL :: gs
        TYPE(systems), INTENT(INOUT), OPTIONAL :: sys
        TYPE(molecular_dynamics), INTENT(INOUT), OPTIONAL:: md
        TYPE(static), INTENT(INOUT), OPTIONAL:: stats
        TYPE(raman), INTENT(INOUT), OPTIONAL:: ram
        TYPE(dipoles), INTENT(INOUT), OPTIONAL:: dip

        ! global settings
        IF (PRESENT(gs)) THEN
            CALL deallocate_global_settings(gs)
        END IF
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

    SUBROUTINE deallocate_global_settings(gs)
        TYPE(global_settings), INTENT(INOUT) :: gs

        !IF (ALLOCATED(gs%md)) DEALLOCATE(gs%md)
        !IF (ALLOCATED(gs%temp)) DEALLOCATE(gs%temp)
        !IF (ALLOCATED(gs%laser_in)) DEALLOCATE(gs%laser_in)

        !IF (ALLOCATED(gs%spectral_type%read_function)) DEALLOCATE(gs%spectral_type%read_function)
        !IF (ALLOCATED(gs%spectral_type%type_input)) DEALLOCATE(gs%spectral_type%type_input)
        !IF (ALLOCATED(gs%spectral_type%type_static)) DEALLOCATE(gs%spectral_type%type_static)
        !IF (ALLOCATED(gs%spectral_type%type_dipole)) DEALLOCATE(gs%spectral_type%type_dipole)

    END SUBROUTINE deallocate_global_settings

    SUBROUTINE deallocate_system(sys)
        TYPE(systems), INTENT(INOUT) :: sys

        !IF (ALLOCATED(sys%natom)) DEALLOCATE(sys%natom)
        !IF (ALLOCATED(sys%framecount)) DEALLOCATE(sys%framecount)
        !IF (ALLOCATED(sys%mol_num)) DEALLOCATE(sys%mol_num)
        !IF (ALLOCATED(sys%filename)) DEALLOCATE(sys%filename)
        !IF (ALLOCATED(sys%frag_type)) DEALLOCATE(sys%frag_type)
        !IF (ALLOCATED(sys%input_mass)) DEALLOCATE(sys%input_mass)
        !IF (ALLOCATED(sys%system)) DEALLOCATE(sys%system)
        !IF (ALLOCATED(sys%periodic)) DEALLOCATE(sys%periodic)
        IF (ALLOCATED(sys%element)) DEALLOCATE (sys%element)
        !IF (ALLOCATED(sys%mass_tot)) DEALLOCATE(sys%mass_tot)
        IF (ALLOCATED(sys%atom_mass_inv_sqrt)) DEALLOCATE (sys%atom_mass_inv_sqrt)
        IF (ALLOCATED(sys%charge)) DEALLOCATE (sys%charge)
        IF (ALLOCATED(sys%mass_atom)) DEALLOCATE (sys%mass_atom)
        IF (ALLOCATED(sys%coord)) DEALLOCATE (sys%coord)
        IF (ALLOCATED(sys%mass_mat)) DEALLOCATE (sys%mass_mat)

        !IF (ALLOCATED(sys%cell%cell_type)) DEALLOCATE(sys%cell%cell_type)
        !IF (ALLOCATED(sys%cell%box_all)) DEALLOCATE(sys%cell%box_all)
        !IF (ALLOCATED(sys%cell%box_x)) DEALLOCATE(sys%cell%box_x)
        !IF (ALLOCATED(sys%cell%box_y)) DEALLOCATE(sys%cell%box_y)
        !IF (ALLOCATED(sys%cell%box_z)) DEALLOCATE(sys%cell%box_z)
        !IF (ALLOCATED(sys%cell%vec(3))) DEALLOCATE(sys%cell%vec(3))
        !IF (ALLOCATED(sys%cell%vec_pbc(3))) DEALLOCATE(sys%cell%vec_pbc(3))

        !IF (ALLOCATED(sys%fragments%nfrag)) DEALLOCATE(sys%fragments%nfrag)
        IF (ALLOCATED(sys%fragments%natom_frag)) DEALLOCATE (sys%fragments%natom_frag)
        IF (ALLOCATED(sys%fragments%fragment)) DEALLOCATE (sys%fragments%fragment)
        !IF (ALLOCATED(sys%fragments%mass_tot_cell)) DEALLOCATE(sys%fragments%mass_tot_cell)
        IF (ALLOCATED(sys%fragments%mass_tot_frag)) DEALLOCATE (sys%fragments%mass_tot_frag)
        IF (ALLOCATED(sys%fragments%refpoint)) DEALLOCATE (sys%fragments%refpoint)

    END SUBROUTINE deallocate_system

    SUBROUTINE deallocate_md(md)
        TYPE(molecular_dynamics), INTENT(INOUT) :: md

        !IF (ALLOCATED(md%trajectory_file)) DEALLOCATE(md%trajectory_file)
        !IF (ALLOCATED(md%velocity_file)) DEALLOCATE(md%velocity_file)
        !IF (ALLOCATED(md%snapshot_time_step)) DEALLOCATE(md%snapshot_time_step)
        !IF (ALLOCATED(md%correlation_depth)) DEALLOCATE(md%correlation_depth)
        !IF (ALLOCATED(md%dt)) DEALLOCATE(md%dt)
        !IF (ALLOCATED(md%freq_range)) DEALLOCATE(md%freq_range)
        !IF (ALLOCATED(md%freq_res)) DEALLOCATE(md%freq_res)
        !IF (ALLOCATED(md%sinc_const)) DEALLOCATE(md%sinc_const)
        IF (ALLOCATED(md%z)) DEALLOCATE (md%z)
        IF (ALLOCATED(md%velos_v)) DEALLOCATE (md%velos_v)
        IF (ALLOCATED(md%v)) DEALLOCATE (md%v)
        IF (ALLOCATED(md%coord_v)) DEALLOCATE (md%coord_v)
        IF (ALLOCATED(md%zhat)) DEALLOCATE (md%zhat)

    END SUBROUTINE deallocate_md

    SUBROUTINE deallocate_stats(stats)
        TYPE(static), INTENT(INOUT):: stats

        INTEGER                     :: xyz, i

        !IF (ALLOCATED(stats%nmodes)) DEALLOCATE(stats%nmodes)
        !IF (ALLOCATED(stats%normal_freq_file)) DEALLOCATE(stats%normal_freq_file)
        !IF (ALLOCATED(stats%normal_displ_file)) DEALLOCATE(stats%normal_displ_file)
        !IF (ALLOCATED(stats%force_file)) DEALLOCATE(stats%force_file)
        !IF (ALLOCATED(stats%dx)) DEALLOCATE(stats%dx)
        IF (ALLOCATED(stats%freq)) DEALLOCATE (stats%freq)
        IF (ALLOCATED(stats%disp)) DEALLOCATE (stats%disp)
        ! Deallocating Static Property
        CALL stats%dealloc_force()
    END SUBROUTINE deallocate_stats

    SUBROUTINE deallocate_dipoles(dips)
        TYPE(dipoles), INTENT(INOUT) :: dips
        !IF (ALLOCATED(dips%dip_file)) DEALLOCATE(dips%dip_file)
        IF (ALLOCATED(dips%dip_dq)) DEALLOCATE (dips%dip_dq)
        IF (ALLOCATED(dips%dipole)) DEALLOCATE (dips%dipole)
        CALL dips%dealloc_dip()

    END SUBROUTINE deallocate_dipoles

    SUBROUTINE deallocate_raman(rams)
        TYPE(raman), INTENT(INOUT):: rams

        INTEGER                     :: xyz

        !IF (ALLOCATED(rams%static_dip_free_file)) DEALLOCATE(rams%static_dip_free_file)
        !IF (ALLOCATED(rams%dip_x_file)) DEALLOCATE(rams%dip_x_file)
        !IF (ALLOCATED(rams%dip_y_file)) DEALLOCATE(rams%dip_y_file)
        !IF (ALLOCATED(rams%dip_z_file)) DEALLOCATE(rams%dip_z_file)
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
        !IF (ALLOCATED(rams%wannier_free)) DEALLOCATE(rams%wannier_free)
        !IF (ALLOCATED(rams%wannier_x)) DEALLOCATE(rams%wannier_x)
        !IF (ALLOCATED(rams%wannier_y)) DEALLOCATE(rams%wannier_y)
        !IF (ALLOCATED(rams%wannier_z))  DEALLOCATE(rams%wannier_z)
        !IF (ALLOCATED(rams%averaging)) DEALLOCATE(rams%averaging)
        !IF (ALLOCATED(rams%direction)) DEALLOCATE(rams%direction)
        IF (ALLOCATED(rams%raman_int)) DEALLOCATE (rams%raman_int)
        !IF (ALLOCATED(rams%pol)) DEALLOCATE (rams%pol)
        IF (ALLOCATED(rams%pol_dq)) DEALLOCATE (rams%pol_dq)

        !IF (ALLOCATED(rams%RR%check_pade)) DEALLOCATE(rams%RR%check_pade)
        !IF (ALLOCATED(rams%RR%framecount_rtp)) DEALLOCATE(rams%RR%framecount_rtp)
        !IF (ALLOCATED(rams%RR%framecount_rtp_pade)) DEALLOCATE(rams%RR%framecount_rtp_pade)
        !IF (ALLOCATED(rams%RR%dt_rtp)) DEALLOCATE(rams%RR%dt_rtp)
        !IF (ALLOCATED(rams%RR%freq_range_rtp)) DEALLOCATE(rams%RR%freq_range_rtp)

        CALL rams%rr%dealloc_rr_all()
        IF (ALLOCATED(rams%RR%zhat_pol_rtp)) DEALLOCATE (rams%RR%zhat_pol_rtp)
    END SUBROUTINE deallocate_raman

END MODULE vib_types
