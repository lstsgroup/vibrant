MODULE vib_types

    USE kinds, ONLY: dp

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: global_settings, systems, molecular_dynamics, static, dipoles, raman, deallocate_all_structures

    !***************************************************************************
    TYPE spectral_type
        CHARACTER(LEN=40)                               :: read_function
        CHARACTER(LEN=40)                               :: type_input, type_static, type_dipole ! From IR calc what are those`?
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

    TYPE resonant_raman
    CHARACTER(LEN=40)                                   :: check_pade
    INTEGER                                             :: framecount_rtp
    INTEGER                                             :: framecount_rtp_pade
    REAL(kind=dp)                                       :: dt_rtp
    REAL(kind=dp)                                       :: dom_rtp
    REAL(kind=dp), DIMENSION(:, :, :, :, :), ALLOCATABLE:: static_dip_rtp
    REAL(kind=dp), DIMENSION(:, :, :, :, :), ALLOCATABLE:: static_dip_x_rtp
    REAL(kind=dp), DIMENSION(:, :, :, :, :), ALLOCATABLE:: static_dip_y_rtp
    REAL(kind=dp), DIMENSION(:, :, :, :, :), ALLOCATABLE:: static_dip_z_rtp
    REAL(kind=dp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE:: pol_rtp
    COMPLEX(kind=dp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: zhat_pol_rtp
    !CHARACTER(LEN=40)                                   :: rtp_dipole_x, rtp_dipole_y, rtp_dipole_z
    !COMPLEX(kind=dp), DIMENSION(:, :), ALLOCATABLE      :: z_iso_resraman, z_aniso_resraman
END TYPE resonant_raman


    !***************************************************************************
    TYPE global_settings
        LOGICAL                                          ::  md !yes/no
        REAL(kind=dp)                                    ::  temp
        REAL(kind=dp)                                    ::  laser_in
        TYPE(spectral_type)                              ::  spectral_type ! global setting of spectral type 'P' , 'IR' , 'R' etc.
    END TYPE global_settings

    !***************************************************************************
    TYPE systems 
        INTEGER                                             :: natom               ! number of atoms
        INTEGER                                             :: framecount          ! number of frames
        INTEGER                                             :: mol_num             ! number of moleces ? 
        CHARACTER(LEN=40)                                   :: filename            ! read in file
        CHARACTER(LEN=40)                                   :: frag_type           !  ???
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
        CHARACTER(LEN=40)                                   :: trajectory_file !maybe not needed should be in system type
        CHARACTER(LEN=40)                                   :: velocity_file   !maybe not needed should be in system type
        REAL(kind=dp)                                       :: snapshot_time_step ! snapshots_time_step equal to dt ?
        REAL(kind=dp)                                       :: correlation_depth ! t_cor needed!
        REAL(kind=dp)                                       :: dt   ! not quite sure ?
        REAL(kind=dp)                                       :: dom ! not quite sure ?
        REAL(kind=dp)                                       :: freq_range ! not sure if right here ?
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
        CHARACTER(LEN=40)                                   :: normal_freq_file     ! file
        CHARACTER(LEN=40)                                   :: normal_displ_file    ! file
        CHARACTER(LEN=40)                                   :: force_file           ! name of force file
        REAL(kind=dp)                                       :: dx                   ! atom displacement
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE            :: freq                 ! ALLOCATE (stats%freq(stats%nmodes))
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE      :: disp                 ! ALLOCATE (stats%disp(stats%nmodes, sys%natom, 3))                 
        REAL(kind=dp), DIMENSION(:, :, :, :, :), ALLOCATABLE:: force                ! ALLOCATE (stats%force(2, sys%natom, 3, sys%natom, 3)) with +-, N atoms, 3 shifts, N atoms, xyz
    END TYPE static
    !***************************************************************************
    TYPE dipoles
        CHARACTER(LEN=40)                                   :: static_dip_file      !
        REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE         :: dip_dq               !
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE      :: dipole               ! ALLOCATE static_dip(sys%natom, 3, 2, 3) is this neeeded?
        REAL(kind=dp), DIMENSION(:, :, :, :), ALLOCATABLE   :: static_dip           ! field free dipole moment
        !REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE      :: dip
        !LOGICAL                                            ::  fragment !<yes/no>
    END TYPE dipoles
    !***************************************************************************
    TYPE raman
        !LOGICAL                                             :: polarizability_type ! <numeric/analytic>
        !numeric
        CHARACTER(LEN=40)                                   :: static_dip_free_file
        CHARACTER(LEN=40)                                   :: static_dip_x_file 
        CHARACTER(LEN=40)                                   :: static_dip_y_file 
        CHARACTER(LEN=40)                                   :: static_dip_z_file        ! Dipolemoments mybe move to dipole class differenes static_dip_free_file and static_dip_file
        REAL(kind=dp), DIMENSION(:, :, :, :), ALLOCATABLE   :: static_dip_free          ! Dipolemoments mybe move to dipole class
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
        REAL(kind=dp), DIMENSION(:, :, :, :, :), ALLOCATABLE:: pol ! ALLOCATE rams%pol(sys%natom, 3, 2, 3, 3)
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE      :: pol_dq !ALLOCATE (rams%pol_dq(stats%nmodes, 3, 3))
        TYPE(resonant_raman)                                :: RR
    END TYPE raman

    !***************************************************************************



CONTAINS
    SUBROUTINE deallocate_all_structures(sys, md, statik, ram, dip)
        IMPLICIT NONE
        TYPE(systems), INTENT(INOUT) :: sys
        TYPE(molecular_dynamics), INTENT(INOUT) :: md
        TYPE(static), INTENT(INOUT) :: statik
        TYPE(raman), INTENT(INOUT) :: ram
        TYPE(dipoles), INTENT(INOUT) :: dip

        ! systems
        IF (ALLOCATED(sys%coord)) DEALLOCATE (sys%coord)
        IF (ALLOCATED(sys%element)) DEALLOCATE (sys%element)
        IF (ALLOCATED(sys%mass_atom)) DEALLOCATE (sys%mass_atom)
        IF (ALLOCATED(sys%mass_mat)) DEALLOCATE (sys%mass_mat)
        IF (ALLOCATED(sys%atom_mass_inv_sqrt)) DEALLOCATE (sys%atom_mass_inv_sqrt)
        IF (ALLOCATED(sys%charge)) DEALLOCATE (sys%charge)

        ! molecular_dynamics
        IF (ALLOCATED(md%z)) DEALLOCATE (md%z)
        IF (ALLOCATED(md%zhat)) DEALLOCATE (md%zhat)
        IF (ALLOCATED(md%coord_v)) DEALLOCATE (md%coord_v)
        IF (ALLOCATED(md%v)) DEALLOCATE (md%v)
        IF (ALLOCATED(md%velos_v)) DEALLOCATE (md%velos_v)

        ! static
        IF (ALLOCATED(statik%freq)) DEALLOCATE (statik%freq)
        IF (ALLOCATED(statik%disp)) DEALLOCATE (statik%disp)
        IF (ALLOCATED(statik%force)) DEALLOCATE (statik%force)

        ! raman
        !IF (ALLOCATED(ram%static_dip_free)) DEALLOCATE (ram%static_dip_free)
        !IF (ALLOCATED(ram%static_dip_x)) DEALLOCATE (ram%static_dip_x)
        !IF (ALLOCATED(ram%static_dip_y)) DEALLOCATE (ram%static_dip_y)
        !IF (ALLOCATED(ram%static_dip_z)) DEALLOCATE (ram%static_dip_z)
        !IF (ALLOCATED(ram%z_iso)) DEALLOCATE (ram%z_iso)
        !IF (ALLOCATED(ram%z_aniso)) DEALLOCATE (ram%z_aniso)
        !IF (ALLOCATED(ram%z_ortho)) DEALLOCATE (ram%z_ortho)
        !IF (ALLOCATED(ram%z_para)) DEALLOCATE (ram%z_para)
        !IF (ALLOCATED(ram%raman_int)) DEALLOCATE (ram%raman_int)
        !IF (ALLOCATED(ram%pol)) DEALLOCATE (ram%pol)
        !IF (ALLOCATED(ram%pol_dq)) DEALLOCATE (ram%pol_dq)
        !IF (ALLOCATED(ram%static_dip_rtp)) DEALLOCATE (ram%static_dip_rtp)
        !IF (ALLOCATED(ram%pol_rtp)) DEALLOCATE (ram%pol_rtp)
        !IF (ALLOCATED(ram%RR%static_dipole_x_rtp)) DEALLOCATE (ram%static_dipole_x_rtp)
        !IF (ALLOCATED(ram%RR%static_dipole_y_rtp)) DEALLOCATE (ram%static_dipole_y_rtp)
        !IF (ALLOCATED(ram%RR%static_dipole_z_rtp)) DEALLOCATE (ram%static_dipole_z_rtp)
        !IF (ALLOCATED(ram%zhat_pol_rtp)) DEALLOCATE (ram%zhat_pol_rtp)

        ! dipoles
        IF (ALLOCATED(dip%static_dip)) DEALLOCATE (dip%static_dip)
        IF (ALLOCATED(dip%dip_dq)) DEALLOCATE (dip%dip_dq)
        IF (ALLOCATED(dip%dipole)) DEALLOCATE (dip%dipole)

    END SUBROUTINE deallocate_all_structures

END MODULE vib_types
