MODULE vib_types

    USE kinds, ONLY: dp

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: global_settings, systems, md, static, dipoles, raman

    !***************************************************************************
    TYPE spectral_type
        CHARACTER(LEN=40)                               :: read_function
        CHARACTER(LEN=40)                               :: type_input, type_static, type_dipole ! From IR calc what are those`?

    END TYPE spectral_type

    !***************************************************************************
    TYPE fragments

    END TYPE fragments

    !***************************************************************************
    TYPE global_settings
        TYPE(spectral_type)                              ::  spectral_type ! global setting of spectral type 'P' , 'IR' , 'R' etc.
        LOGICAL                                          ::  md !yes/no
        REAL(kind=dp)                                    ::  temp
        REAL(kind=dp)                                    ::  laser_in
    END TYPE global_settings

    !***************************************************************************
    TYPE systems
        CHARACTER(LEN=40)                                :: filename
        CHARACTER(LEN=40)                               :: system, periodic
        INTEGER                                          :: natom, framecount, mol_num, framecount_rtp
        !LOGICAL                                          ::  periodic !yes/no
        CHARACTER(LEN=40)                                :: cell_type
        REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE      :: coord
        CHARACTER(LEN=2), DIMENSION(:), ALLOCATABLE      :: element
        REAL(kind=dp):: mass_tot
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: mass_atom
        REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE   :: mass_mat
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE     :: atom_mass_inv_sqrt, charge
        REAL(kind=dp)                                  :: box_all, box_x, box_y, box_z, vec(3), vec_pbc(3)
        TYPE(fragments)                    ::  fragments
    END TYPE systems

    !***************************************************************************
    TYPE md
        CHARACTER(LEN=40)             :: trajectory_file !maybe not needed should be in system type
        CHARACTER(LEN=40)             :: velocity_file   !maybe not needed should be in system type
        REAL(kind=dp)             :: snapshot_time_step ! snapshots_time_step equal to dt ?
        REAL(kind=dp)             :: correlation_depth ! t_cor needed!
        REAL(kind=dp)                          :: dt, dom ! not quite sure ?
        REAL(kind=dp)::freq_range                           ! not sure if right here ?
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: z       ! correlations vector ?
        COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE    :: zhat ! correlations vector hat ?
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE     :: coord_v ! coordinates vector
        REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE        :: velos_v ! velocity vector
    END TYPE md
    !***************************************************************************
    TYPE static
        CHARACTER(LEN=40)                                   :: normal_freq_file, normal_displ_file, force_file
        REAL(kind=dp)                                       :: dx
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE            :: freq
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE              :: disp
        REAL(kind=dp), DIMENSION(:, :, :, :, :), ALLOCATABLE           :: force
        INTEGER                                             :: nmodes

        INTEGER :: framecount_rtp
    END TYPE static
    !***************************************************************************
    TYPE dipoles
        CHARACTER(LEN=40)                               :: static_dip_file
        REAL(kind=dp), DIMENSION(:, :, :, :), ALLOCATABLE    :: static_dip
        LOGICAL                                          ::  fragment !<yes/no>
        REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE     :: dip_dq
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE       :: dipole

    END TYPE dipoles
    !***************************************************************************
    TYPE raman
        LOGICAL :: polarizability_type ! <numeric/analytic>
        !numeric
        CHARACTER(LEN=40)  :: static_dip_free_file, static_dip_x_file, static_dip_y_file, static_dip_z_file
        REAL(kind=dp), DIMENSION(:, :, :, :), ALLOCATABLE     :: static_dip_free, static_dip_x, static_dip_y, static_dip_z
        !end numeric
        ! analytic
        CHARACTER(LEN=40)                               :: static_pol_file
        !end analytic

        REAL(kind=dp), DIMENSION(:), ALLOCATABLE      :: z_iso, z_aniso, z_ortho, z_para
        CHARACTER(LEN=40)                       :: wannier_free, wannier_x, wannier_y, wannier_z
        CHARACTER(LEN=40)                       :: averaging, direction
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: raman_int
        REAL(kind=dp), DIMENSION(:, :, :, :, :), ALLOCATABLE  :: pol
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE :: pol_dq

        !resonant_raman
        REAL(kind=dp), DIMENSION(:, :, :, :, :), ALLOCATABLE :: static_dip_rtp
        REAL(kind=dp)                                  ::  dt_rtp
        REAL(kind=dp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE  :: pol_rtp
        REAL(kind=dp), DIMENSION(:, :, :, :, :), ALLOCATABLE :: static_dipole_x_rtp, static_dipole_y_rtp, static_dipole_z_rtp
        CHARACTER(LEN=40)                               :: check_pade
        REAL(kind=dp)                                       :: dom_rtp
        COMPLEX(kind=dp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE          :: zhat_pol_rtp
        !end resonant_raman
    END TYPE raman

    !***************************************************************************

END MODULE vib_types
