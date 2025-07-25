PROGRAM vib2d

    USE, INTRINSIC           :: ISO_C_BINDING
    USE kinds, ONLY: dp, str_len
    USE constants, ONLY: speed_light, const_planck, const_permit, pi, const_charge, const_boltz, damping_constant, joule_unit, &
                         debye, ev_unit, action_unit, bohr2ang, hartreebohr2evang, at_u, ang, fs2s, reccm2ev, t_cor, temp, &
                         hessian_factor
    USE read_input, ONLY:  parse_command_line,  parse_input, check_input
    USE vib_types, ONLY: global_settings, systems, molecular_dynamics, static, dipoles, raman,init_global_settings,init_systems, init_molecular_dynamics,init_static, deallocate_types
    USE setup, ONLY: read_input, masses_charges, conversion, pbc_orthorombic, pbc_hexagonal
    USE read_traj, ONLY: read_coord, read_coord_frame, read_normal_modes, read_static, read_static_resraman
    USE dipole_calc, ONLY: center_mass, wannier, wannier_frag, solv_frag_index
    USE vel_cor, ONLY: cvv, cvv_iso, cvv_aniso, cvv_only_x, cvv_resraman
    USE fin_diff, ONLY: central_diff, forward_diff, finite_diff_static, finite_diff_static_resraman
    USE calc_spectra, ONLY: spec_power, normal_mode_analysis, spec_static_ir, spec_static_raman, spec_ir, spec_raman, spec_abs, spec_static_resraman!, spec_resraman
    USE omp_lib, ONLY: omp_get_num_threads
    USE config_info, ONLY: output_config_info

    IMPLICIT NONE

    INCLUDE 'fftw3.f03'

    INTEGER                                         :: b, i, j, k, natom, framecount, framecount_rtp_pade, t0, t1
    INTEGER                                         :: frm, nu, tau, stat, mol_num, nmodes, framecount_rtp, nfrag
    INTEGER                                         :: count_0, count_1, count_rate, count_max, num_threads
    INTEGER(kind=dp)                                 :: plan
    INTEGER, DIMENSION(:), ALLOCATABLE                :: natom_frag, natom_frag_x, natom_frag_free, nfrag_BO, nfrag_BC, nfrag_Ph
    INTEGER, DIMENSION(:, :, :), ALLOCATABLE            :: fragment
    CHARACTER(LEN=40)                               :: system, filename, static_pol_file, read_function, length, type_input, output_dip
    CHARACTER(LEN=40)                               :: filename_dip, type_dipole, direction, averaging, cell_type, coord_file
    CHARACTER(LEN=40)                               :: normal_freq_file, normal_displ_file, frag_type, type_static, force_file
    CHARACTER(LEN=40)                               :: static_dip_free_file, static_dip_x_file, static_dip_y_file, static_dip_z_file
    CHARACTER(LEN=40)                               :: rtp_dipole_x, rtp_dipole_y, rtp_dipole_z, check_pade, charac
    CHARACTER(LEN=40)                               :: output_dip_free, output_dip_x, output_dip_y, output_dip_z
    CHARACTER(LEN=40)                               :: dipole_free, dipole_x, dipole_y, dipole_z, output_findif_dip
    CHARACTER(LEN=40)                               :: output_findif_x, output_findif_y, output_findif_z, input_mass
    CHARACTER(LEN=40)                               :: wannier_free, wannier_x, wannier_y, wannier_z, periodic
    CHARACTER(LEN=2), DIMENSION(:), ALLOCATABLE       :: element
    REAL(kind=dp)                                    :: dist, box_all, box_x, box_y, box_z, mass_tot
    REAL(kind=dp)                                    :: freq_range, dom, ce, co, h_kbT, raman_eq, a, dt_rtp, dom_rtp
    REAL(kind=dp)                                    ::  laser_in
    REAL(kind=dp)                                    :: f, tmax, omega, theta, sinth, costh, sinsq
    REAL(kind=dp)                                    :: cossq, thsq, thcub, alpha, beta, gamma0, dt, multiplier, dx
    REAL(kind=dp)                                    :: time_init, time_final, elapsed_time
    REAL(kind=dp)                                    ::  sinc_const, mass_tot_cell
    REAL(kind=dp), DIMENSION(3)                       :: vec, vec_pbc, coord2, coord1
    REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE       :: refpoint, refpoint_free, refpoint_x, refpoint_y, refpoint_z
    REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE       :: alpha_resraman_x, alpha_resraman_y, alpha_resraman_z
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE           :: kissfft, z, norm, mass_atom, z_aniso, z_iso, z_ortho, z_para, zhat_depol
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE           :: atom_mass_inv_sqrt, charge
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE           :: zhat_para_all, zhat_depol_x, zhat_unpol_x, freq, raman_int, test_x
    REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE         :: test_x2, mass_mat
    COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE        :: zhat, zhat_aniso, zhat_iso, zhat_para, zhat_ortho, zhat_unpol
    COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE        :: integral, zhat_test, y_out
    COMPLEX(kind=dp), DIMENSION(:, :), ALLOCATABLE      :: zhat_test2, z_iso_resraman, z_aniso_resraman
    COMPLEX(kind=dp), DIMENSION(:, :, :), ALLOCATABLE    :: zhat_resraman_x
    COMPLEX(kind=dp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE        :: zhat_pol_rtp
    REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE         :: coord, mass, dipole2, refpoint2, mass_tot_frag, dip_dq
    REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE       :: dip, dip_free, dip_x, dip_y, dip_z
    REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE       :: disp, pol_dq
    REAL(kind=dp), DIMENSION(:, :, :, :), ALLOCATABLE     :: com, pol_dq_rtp
    REAL(kind=dp), DIMENSION(:, :, :, :), ALLOCATABLE     :: static_dip_x, static_dip_y, static_dip_z, static_dip_free
    REAL(kind=dp), DIMENSION(:, :, :, :, :), ALLOCATABLE   :: pol, force
    REAL(kind=dp), DIMENSION(:, :, :, :, :), ALLOCATABLE       :: static_dip_x_rtp, static_dip_y_rtp, static_dip_z_rtp
    REAL(kind=dp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: pol_rtp
    REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE       :: coord_v, coord_dip, coord_f, coord_x, coord_y, coord_z, dipole
    REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE       :: coord_v_x, coord_v_free, alpha_x, alpha_y, alpha_z, v
    REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE       :: alpha_diff_x, alpha_diff_y, alpha_diff_z
    
    ! Variables of your derived types:
    TYPE(global_settings) :: gs
    TYPE(systems)        :: sys
    TYPE(molecular_dynamics)    :: md

    TYPE(static)            :: stats
    TYPE(dipoles)            :: dips
    TYPE(raman)            :: rams
    CHARACTER(LEN=str_len)          :: input_file_name
    
    !$omp parallel
    num_threads = omp_get_num_threads()
    !$omp end parallel
    
    CALL SYSTEM_CLOCK(count_0, count_rate, count_max) !Starting time
    time_init = count_0*1.0_dp/count_rate
    
    CALL init_global_settings(gs)
    CALL init_systems(sys)
    CALL init_molecular_dynamics(md)
    CALL init_static(stats)
    
    CALL output_config_info()

    !***************************************************************************
    !                                                                        ***
    !                                                                           ***
    !                                                                               ***

    CALL parse_command_line(input_file_name)
    
    CALL parse_input(gs, sys, md, stats, dips, rams, input_file_name)
    PRINT *, REPEAT('-', 30)
    CALL check_input(gs, sys, md, stats, dips, rams)
    PRINT *, REPEAT('-', 30)
    !write(*,*) "input_file_name", input_file_name
    !write(*,*) "temperature", gs%temp
    !write(*,*) "laserin", gs%laser_in
    !write(*,*) "read_function", gs%spectral_type%read_function
    !write(*,*) "type_input ", gs%spectral_type%type_input
    !write(*,*) "type_static ", gs%spectral_type%type_static
    !write(*,*) "type_dipole ", gs%spectral_type%type_dipole
    !                                                                                  ***
    !                                                                           ***
    !                                                                        ***
    !***************************************************************************




    
!    CALL read_input(filename, static_pol_file, static_dip_free_file, static_dip_x_file, static_dip_y_file, static_dip_z_file, &
!                    normal_freq_file, normal_displ_file, read_function, system, length, box_all, box_x, box_y, box_z, dt, &
!                    type_input, wannier_free, wannier_x, wannier_y, wannier_z, input_mass, periodic, direction, averaging, &
!                    type_dipole, cell_type, rtp_dipole_x, rtp_dipole_y, rtp_dipole_z, framecount_rtp, dt_rtp, &
!                    frag_type, type_static, force_file, laser_in, check_pade, dx, framecount_rtp_pade)
    
    CALL conversion(md%dt, md%dom, rams%RR%dt_rtp, rams%RR%dom_rtp, md%freq_range, md%sinc_const)
!
!    ! TEMP setup
!    gs%spectral_type%read_function = read_function
!    gs%spectral_type%type_dipole = type_dipole
!    gs%spectral_type%type_static = type_static
!    gs%spectral_type%type_input = type_input
!    gs%laser_in = laser_in
!
!    sys%filename = filename
!    sys%system = system
!    sys%cell%cell_type = cell_type
!    sys%periodic = periodic
!    sys%frag_type = frag_type
!    sys%input_mass = input_mass
!    sys%cell%box_all = box_all
!    sys%cell%box_x = box_x
!    sys%cell%box_y = box_y
!    sys%cell%box_z = box_z
!    
!    stats%force_file = force_file
!    stats%normal_freq_file = normal_freq_file
!    stats%normal_displ_file = normal_displ_file
!    stats%dx = dx
!    
!    dips%static_dip_file = static_dip_free_file
!    
!    rams%static_pol_file = static_pol_file
!    rams%wannier_free = wannier_free
!    rams%wannier_x = wannier_x
!    rams%wannier_y = wannier_y
!    rams%wannier_z = wannier_z
!    rams%averaging = averaging
!    rams%direction = direction
!
!    rams%RR%dt_rtp = dt_rtp
!    rams%RR%dom_rtp = dom_rtp
!    rams%RR%check_pade = check_pade
!    rams%RR%framecount_rtp_pade = framecount_rtp_pade
!    rams%RR%framecount_rtp = framecount_rtp
!
!    
!    md%freq_range = freq_range
!    md%dt = dt
!    md%dom = dom
!    md%sinc_const = sinc_const
!
!!
!!    !***************************************************************************
    IF (gs%spectral_type%read_function=='P') THEN
        CALL read_coord(gs, sys)
        CALL masses_charges(gs, sys)
             print*,'test'  
        CALL spec_power(gs, sys, md)
!        !***************************************************************************
!        !***************************************************************************
!    ELSEIF (gs%spectral_type%read_function=='MD-IR') THEN
!        CALL read_coord(gs, sys)
!        IF (sys%system=='1' .OR. sys%system=='2' .AND. gs%spectral_type%type_dipole=='1') THEN !!fragment approach or whole supercell
!            IF (sys%cell%cell_type=='1' .OR. sys%cell%cell_type=='2') THEN !!KP or SC
!                CALL masses_charges(gs, sys)
!            END IF
!        END IF
!
!        CALL spec_ir(gs, sys, md, dips)
!        !***************************************************************************
!
!        !***************************************************************************
!    ELSEIF (gs%spectral_type%read_function=='MD-R') THEN
!        sys%filename = wannier_free! <----  MUST BE ADJUSTED
!        CALL read_coord(gs, sys)
!        CALL masses_charges(gs, sys)
!
!        CALL spec_raman(gs, sys, md, dips, rams)
!        !***************************************************************************
!
!        !***************************************************************************
    ELSEIF (gs%spectral_type%read_function=='NMA') THEN
        CALL read_coord(gs, sys)
        CALL masses_charges(gs, sys)
        CALL read_normal_modes(gs, sys, stats)

        CALL normal_mode_analysis(sys, stats)
!        !***************************************************************************
!
!        !***************************************************************************
!    ELSEIF (gs%spectral_type%read_function=='IR') THEN
!        CALL read_coord(gs, sys)
!        CALL masses_charges(gs, sys)
!        CALL read_normal_modes(gs, sys, stats)
!        CALL read_static(dips%static_dip_file, dips%static_dip, gs, sys, rams)
!        IF (type_static=='1') THEN
!            CALL normal_mode_analysis(sys, stats)
!        END IF
!        CALL finite_diff_static(gs, sys, stats, dips, rams)
!
!        CALL spec_static_ir(sys, stats, dips)
!        !***************************************************************************
!
!        !***************************************************************************
!    ELSEIF (gs%spectral_type%read_function=='R') THEN
!        CALL read_coord(gs, sys)
!        CALL masses_charges(gs, sys)
!        CALL read_normal_modes(gs, sys, stats)
!        CALL read_static(dips%static_dip_file, dips%static_dip, gs, sys, rams)
!        IF (type_dipole=='2') THEN
!            CALL read_static(static_dip_x_file, static_dip_x, gs, sys, rams)
!            CALL read_static(static_dip_y_file, static_dip_y, gs, sys, rams)
!            CALL read_static(static_dip_z_file, static_dip_z, gs, sys, rams)
!        END IF
!        IF (type_static=='1') THEN
!            CALL normal_mode_analysis(sys, stats)
!        END IF
!        CALL finite_diff_static(gs, sys, stats, dips, rams)
!
!        CALL spec_static_raman(gs, sys, stats, dips, rams)
!        !***************************************************************************
!
!        !***************************************************************************
!    ELSEIF (gs%spectral_type%read_function=='ABS') THEN
!        CALL read_coord(gs, sys)
!        CALL read_static_resraman(static_dip_x_file, rams%RR%static_dip_x_rtp, sys, rams)
!        CALL read_static_resraman(static_dip_y_file, rams%RR%static_dip_y_rtp, sys, rams)
!        CALL read_static_resraman(static_dip_z_file, rams%RR%static_dip_z_rtp, sys, rams)
!        CALL finite_diff_static_resraman(rams%RR%static_dip_x_rtp, rams%RR%static_dip_y_rtp, rams%RR%static_dip_z_rtp, sys, rams) !<-- CHANGE ?
!
!        CALL spec_abs(gs, sys, rams)
!        !***************************************************************************
!
!        !***************************************************************************
!    ELSEIF (gs%spectral_type%read_function=='RR') THEN
!        CALL read_coord(gs, sys, rams)
!        CALL masses_charges(gs, sys)
!        CALL read_normal_modes(gs, sys, stats)
!        CALL read_static_resraman(static_dip_x_file, rams%RR%static_dip_x_rtp, sys, rams)
!        CALL read_static_resraman(static_dip_y_file, rams%RR%static_dip_y_rtp, sys, rams)
!        CALL read_static_resraman(static_dip_z_file, rams%RR%static_dip_z_rtp, sys, rams)
!        IF (type_static=='1') THEN
!            CALL normal_mode_analysis(sys, stats)
!        END IF
!
!        CALL finite_diff_static_resraman(rams%RR%static_dip_x_rtp, rams%RR%static_dip_y_rtp, rams%RR%static_dip_z_rtp, sys, rams)
!        CALL spec_abs(gs, sys, rams)
!
!        CALL spec_static_resraman(gs, sys, stats, rams)
!        !***************************************************************************
!        !***************************************************************************
!    ELSEIF (gs%spectral_type%read_function=='MD-RR') THEN
!        !sys%filename = rtp_dipole_x ! <----  MUST BE ADJUSTED
!        !CALL read_coord(gs, sys)
!        !natom = sys%natom
!        !framecount = sys%framecount
!        !mol_num = sys%mol_num
!        !CALL spec_resraman(natom,framecount,element,rtp_dipole_x,rtp_dipole_y,rtp_dipole_z,type_input,mol_num,system,&
!        !     read_function,dt,z_iso_resraman,z_aniso_resraman,dom,dom_rtp,laser_in_resraman,y_out)
    END IF
!
    CALL deallocate_types(gs, sys, md, stats, rams, dips)

    CALL SYSTEM_CLOCK(count_1, count_rate, count_max) !Ending time
    time_final = count_1*1.0_dp/count_rate
    elapsed_time = time_final - time_init !Elapsed time

    WRITE (*, 1003) INT(elapsed_time), elapsed_time - INT(elapsed_time) !Write elapsed time
1003 FORMAT('  Wall Clock = ', i0, F0.9)

END PROGRAM vib2d
