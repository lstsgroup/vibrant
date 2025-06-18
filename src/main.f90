PROGRAM vib2d

    USE, INTRINSIC           :: ISO_C_BINDING
    USE kinds, ONLY: dp
    USE constants, ONLY: speed_light, const_planck, const_permit, pi, const_charge, const_boltz, damping_constant, joule_unit, &
                         debye, ev_unit, action_unit, bohr2ang, hartreebohr2evang, at_u, ang, fs2s, reccm2ev, t_cor, temp, &
                         hessian_factor
    USE vib_types, ONLY: global_settings, systems, md, static, dipoles, raman
    USE setup, ONLY: read_input, masses_charges, conversion, pbc_orthorombic, pbc_hexagonal
    USE read_traj, ONLY: read_coord, read_coord_frame, read_normal_modes, read_static, read_static_resraman
    USE dipole_calc, ONLY: center_mass, wannier, wannier_frag, solv_frag_index
    USE vel_cor, ONLY: cvv, cvv_iso, cvv_aniso, cvv_only_x, cvv_resraman
    USE fin_diff, ONLY: central_diff, forward_diff, finite_diff_static, finite_diff_static_resraman
    USE calc_spectra, ONLY: spec_power, spec_ir, spec_raman, normal_mode_analysis, spec_static_ir, spec_static_raman, &
                            spec_abs, spec_static_resraman, spec_resraman
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
    TYPE(static)            :: stats
    TYPE(dipoles)            :: dips
    TYPE(raman)            :: rams
!$omp parallel
    num_threads = omp_get_num_threads()
!$omp end parallel

    CALL SYSTEM_CLOCK(count_0, count_rate, count_max) !Starting time
    time_init = count_0*1.0_dp/count_rate

    CALL output_config_info()

    CALL read_input(filename, static_pol_file, static_dip_free_file, static_dip_x_file, static_dip_y_file, static_dip_z_file, &
                    normal_freq_file, normal_displ_file, read_function, system, length, box_all, box_x, box_y, box_z, dt, &
                    type_input, wannier_free, wannier_x, wannier_y, wannier_z, input_mass, periodic, direction, averaging, &
                    type_dipole, cell_type, rtp_dipole_x, rtp_dipole_y, rtp_dipole_z, framecount_rtp, dt_rtp, &
                    frag_type, type_static, force_file, laser_in, check_pade, dx, framecount_rtp_pade)

    CALL conversion(dt, dom, dt_rtp, dom_rtp, freq_range, sinc_const)

    ! TEMP setup
    gs%spectral_type%read_function = read_function
    gs%spectral_type%type_dipole = type_dipole
    gs%spectral_type%type_static = type_static
    gs%laser_in = laser_in

    sys%filename = filename
    sys%system = system
    sys%framecount_rtp = framecount_rtp
    sys%periodic = periodic

    stats%force_file = force_file
    stats%normal_freq_file = normal_freq_file
    stats%normal_displ_file = normal_displ_file
    stats%dx = dx

    dips%static_dip_file = static_dip_free_file

    rams%static_pol_file = static_pol_file
!
!    !***************************************************************************
!    IF (read_function=='P') THEN
!        CALL read_coord(natom, framecount, element, coord, filename, periodic, mol_num, system, read_function, framecount_rtp, type_dipole)
!        CALL masses_charges(natom, mass_atom, atom_mass_inv_sqrt, mass_mat, element, mass_tot, charge)
!        CALL spec_power(z, zhat, type_input, freq_range, natom, framecount, dt, element, filename, coord_v, v, &
!                        input_mass, dom, mass_atom, read_function, mol_num, mass_tot, coord, system, frag_type)
!        DEALLOCATE (element, z, zhat, coord_v, mass_atom, coord)
!        !***************************************************************************
!
!        !***************************************************************************
!    ELSEIF (read_function=='MD-IR') THEN
!        IF (system=='1' .OR. system=='2' .AND. type_dipole=='1') THEN !!fragment approach or whole supercell
!            IF (cell_type=='1' .OR. cell_type=='2') THEN !!KP or SC
!                CALL read_coord(natom, framecount, element, coord, filename, periodic, mol_num, system, read_function, framecount_rtp, type_dipole)
!                CALL masses_charges(natom, mass_atom, atom_mass_inv_sqrt, mass_mat, element, mass_tot, charge)
!                CALL read_coord_frame(natom, framecount, element, filename, coord_v)
!                !  CALL wannier(element,filename,natom,box_all,box_x,box_y,box_z,vec,vec_pbc,debye,mol_num,&
!                !        periodic,mass_tot,framecount,mass_atom,coord_v,dip)
!                !    CALL frag_index(natom,filename,element,coord_v,fragment,natom_frag,framecount)
!                CALL center_mass(natom_frag, natom, refpoint, coord_v, filename, element, box_all, box_x, box_y, &
!                                 box_z, vec, vec_pbc, fragment, mass_atom, framecount, cell_type, mass_tot_frag, frag_type, mol_num, &
!                                 nfrag, type_dipole, system, mass_tot_cell)
!                CALL wannier_frag(natom_frag, filename, natom, element, coord_v, box_all, box_x, box_y, box_z, vec, vec_pbc, dipole, &
!                                  refpoint, fragment, framecount, mass_tot_frag, mol_num, system, type_dipole, charge, mass_tot_cell)
!                CALL spec_ir(z, zhat, freq_range, natom, framecount, dt, element, filename, coord_v, v, input_mass, dom, &
!                             mol_num, box_all, box_x, box_y, box_z, vec, vec_pbc, periodic, mass_tot, mass_atom, type_input, dip, &
!                             read_function, coord, type_dipole, dipole, system, mass_tot_frag, sinc_const, nfrag, frag_type)
!            ELSEIF (cell_type=='3') THEN !!SC with solvent
!                CALL read_coord(natom, framecount, element, coord, filename, periodic, mol_num, system, read_function, framecount_rtp, type_dipole)
!                CALL read_coord_frame(natom, framecount, element, filename, coord_v)
!                CALL solv_frag_index(natom, coord_v, filename, element, box_all, vec, vec_pbc, &
!                                     box_x, box_y, box_z, mass_atom, framecount, cell_type, refpoint, natom_frag, fragment, mass_tot_frag)
!                CALL wannier_frag(natom_frag, filename, natom, element, coord_v, box_all, box_x, box_y, box_z, vec, vec_pbc, dipole, &
!                                  refpoint, fragment, framecount, mass_tot_frag, mol_num, system, type_dipole, charge, mass_tot_cell)
!                CALL spec_ir(z, zhat, freq_range, natom, framecount, dt, element, filename, coord_v, v, input_mass, dom, &
!                             mol_num, box_all, box_x, box_y, box_z, vec, vec_pbc, periodic, mass_tot, mass_atom, type_input, dip, &
!                             read_function, coord, type_dipole, dipole, system, mass_tot_frag, sinc_const, nfrag, frag_type)
!            END IF
!        ELSEIF (system=='2') THEN !!molecular approach (Berry phase)
!            CALL read_coord(natom, framecount, element, coord, filename, periodic, mol_num, system, read_function, framecount_rtp, type_dipole)
!            CALL spec_ir(z, zhat, freq_range, natom, framecount, dt, element, filename, coord_v, v, input_mass, dom, &
!                         mol_num, box_all, box_x, box_y, box_z, vec, vec_pbc, periodic, mass_tot, mass_atom, type_input, dip, &
!                         read_function, coord, type_dipole, dipole, system, mass_tot_frag, sinc_const, nfrag, frag_type)
!        END IF
!        !***************************************************************************
!
!        !***************************************************************************
!    ELSEIF (read_function=='MD-R') THEN
!        CALL read_coord(natom, framecount, element, coord, wannier_free, periodic, mol_num, system, read_function, framecount_rtp, type_dipole)
!        CALL masses_charges(natom, mass_atom, atom_mass_inv_sqrt, mass_mat, element, mass_tot, charge)
!        CALL spec_raman(natom, framecount, element, coord, wannier_free, wannier_x, wannier_y, wannier_z, mass_atom, mass_tot, periodic, &
!                        mol_num, dt, dom, coord_v, v, type_input, box_all, box_x, box_y, box_z, vec, vec_pbc, read_function, &
!                        z_iso, z_aniso, z_ortho, z_para, laser_in, filename, averaging, direction, &
!                        type_dipole, system, natom_frag, fragment, refpoint, dipole, cell_type, mass_tot_frag, frag_type, nfrag, charge, mass_tot_cell)
!        DEALLOCATE (coord, coord_v, mass_atom, element)
!        !***************************************************************************
!
!        !***************************************************************************
    IF (read_function=='NMA') THEN
        !CALL read_coord(natom, framecount, element, coord, filename, periodic, mol_num, system, read_function, &
        !                framecount_rtp, type_dipole)
        !CALL masses_charges(natom, mass_atom, atom_mass_inv_sqrt, mass_mat, element, mass_tot, charge)
        !CALL read_normal_modes(natom, element, normal_freq_file, normal_displ_file, freq, disp, nmodes, &
        !                         read_function, type_static, force_file, force)
        !CALL normal_mode_analysis(natom, force, dx, mass_mat, nmodes, freq, disp)

        CALL read_coord(gs, sys)
        CALL masses_charges(gs, sys)
        CALL read_normal_modes(gs, sys, stats)
        CALL normal_mode_analysis(sys, stats)

        !***************************************************************************

        !***************************************************************************
    ELSEIF (read_function=='IR') THEN
        !static_pol_file = sys%static_pol_file
        !pol = sys%pol
        !static_dip_free_file = sys%static_dip_free_file
        !type_dipole = sys%type_dipole
        !static_dip_free = dips%static_dip_free
        !type_static = sys%type_static
        !natom = sys%natom
        !!element = sys%element
        !!coord = sys%coord
        !!mass_atom = sys%mass_atom
        !nmodes = stats%nmodes
        !disp = stats%disp
        !atom_mass_inv_sqrt = sys%atom_mass_inv_sqrt
        !dx = stats%dx
        !pol = rams%pol
        !static_dip_free = dips%static_dip
        !freq = stats%freq
        !dips%dip_dq= dip_dq
        CALL read_coord(gs, sys)
        CALL masses_charges(gs, sys)
        CALL read_normal_modes(gs, sys, stats)
        CALL read_static(gs, sys, dips, rams)
        IF (type_static=='1') THEN
            !CALL normal_mode_analysis(natom, force, dx, mass_mat, nmodes, freq, disp)
            CALL normal_mode_analysis(sys, stats)
        END IF
        CALL finite_diff_static(gs, sys, stats, dips, rams)
        CALL spec_static_ir(sys, stats, dips)

        DEALLOCATE (sys%element, sys%coord, sys%mass_atom)
        !***************************************************************************

        !***************************************************************************
    ELSEIF (read_function=='R') THEN
        !CALL read_coord(natom, framecount, element, coord, filename, periodic, mol_num, system, read_function, framecount_rtp, type_dipole)
        !CALL masses_charges(natom, mass_atom, atom_mass_inv_sqrt, mass_mat, element, mass_tot, charge)
        !CALL read_normal_modes(natom, element, normal_freq_file, normal_displ_file, freq, disp, nmodes, &
        !                       read_function, type_static, force_file, force)
        !CALL read_static(natom, element, static_pol_file, pol, static_dip_free_file, type_dipole, static_dip_free, type_static)
        CALL read_coord(gs, sys)
        CALL masses_charges(gs, sys)
        CALL read_normal_modes(gs, sys, stats)
        CALL read_static(gs, sys, dips, rams)
        IF (type_dipole=='2') THEN
            CALL read_static(gs, sys, dips, rams) ! <---- THIS SHOULD and DOES NOT WORK and MUST BE ADJUSTED
            CALL read_static(gs, sys, dips, rams) ! <---- THIS SHOULD and DOES NOT WORK and MUST BE ADJUSTED
            CALL read_static(gs, sys, dips, rams) ! <---- THIS SHOULD and DOES NOT WORK and MUST BE ADJUSTED
            !CALL read_static(natom, element, static_pol_file, pol, static_dip_x_file, type_dipole, static_dip_x, type_static)
            !CALL read_static(natom, element, static_pol_file, pol, static_dip_y_file, type_dipole, static_dip_y, type_static)
            !CALL read_static(natom, element, static_pol_file, pol, static_dip_z_file, type_dipole, static_dip_z, type_static)
        END IF
        IF (type_static=='1') THEN
            !CALL normal_mode_analysis(natom, force, dx, mass_mat, nmodes, freq, disp)
            CALL normal_mode_analysis(sys, stats)
        END IF
        !CALL finite_diff_static(natom, nmodes, pol, pol_dq, disp, atom_mass_inv_sqrt, dx, static_dip_free, static_dip_x, &
        !                        static_dip_y, static_dip_z, type_dipole, read_function, dip_dq)
        CALL finite_diff_static(gs, sys, stats, dips, rams)
        !static_pol_file = sys%static_pol_file
        !pol = sys%pol
        !static_dip_free_file = sys%static_dip_free_file
        !type_dipole = sys%type_dipole
        !static_dip_free = dips%static_dip_free
        !type_static = sys%type_static
        !natom = sys%natom
        !element = sys%element
        !coord = sys%coord
        !mass_atom = sys%mass_atom
        !nmodes = stats%nmodes
        !disp = stats%disp
        !atom_mass_inv_sqrt = sys%atom_mass_inv_sqrt
        !dx = stats%dx
        !pol = rams%pol
        !static_dip_free = dips%static_dip
        !freq = stats%freq
        !dips%dip_dq= dip_dq
        !pol_dq = ramns%pol_dq
        CALL spec_static_raman(gs, sys, stats, dips, rams) 
        !CALL spec_static_raman(nmodes, pol_dq, laser_in, freq, raman_int, element, coord, disp, natom)
        DEALLOCATE (stats%freq, stats%disp)
        DEALLOCATE (rams%pol_dq)
        DEALLOCATE (sys%element, sys%coord, sys%mass_atom)
        !***************************************************************************
!
!        !***************************************************************************
!    ELSEIF (read_function=='ABS') THEN
!        CALL read_coord(natom, framecount, element, coord, filename, periodic, mol_num, system, read_function, framecount_rtp, type_dipole)
!        CALL read_static_resraman(natom, element, static_dip_x_file, framecount_rtp, static_dip_x_rtp)
!        CALL read_static_resraman(natom, element, static_dip_y_file, framecount_rtp, static_dip_y_rtp)
!        CALL read_static_resraman(natom, element, static_dip_z_file, framecount_rtp, static_dip_z_rtp)
!        CALL finite_diff_static_resraman(natom, pol_rtp, static_dip_x_rtp, static_dip_y_rtp, &
!                                         static_dip_z_rtp, framecount_rtp, dt_rtp)
!        CALL spec_abs(nmodes, natom, pol_rtp, freq, framecount_rtp, framecount_rtp_pade, check_pade, &
!                      dom_rtp, read_function, zhat_pol_rtp)
!        DEALLOCATE (zhat_pol_rtp)
!        !***************************************************************************
!
!        !***************************************************************************
!    ELSEIF (read_function=='RR') THEN
!        CALL read_coord(natom, framecount, element, coord, filename, periodic, mol_num, system, read_function, framecount_rtp, type_dipole)
!        CALL masses_charges(natom, mass_atom, atom_mass_inv_sqrt, mass_mat, element, mass_tot, charge)
!        CALL read_normal_modes(natom, element, normal_freq_file, normal_displ_file, freq, disp, nmodes, &
!                               read_function, type_static, force_file, force)
!        CALL read_static_resraman(natom, element, static_dip_x_file, framecount_rtp, static_dip_x_rtp)
!        CALL read_static_resraman(natom, element, static_dip_y_file, framecount_rtp, static_dip_y_rtp)
!        CALL read_static_resraman(natom, element, static_dip_z_file, framecount_rtp, static_dip_z_rtp)
!        IF (type_static=='1') THEN
!            CALL normal_mode_analysis(natom, force, dx, mass_mat, nmodes, freq, disp)
!        END IF
!        CALL finite_diff_static_resraman(natom, pol_rtp, static_dip_x_rtp, static_dip_y_rtp, &
!                                         static_dip_z_rtp, framecount_rtp, dt_rtp)
!        CALL spec_abs(nmodes, natom, pol_rtp, freq, framecount_rtp, framecount_rtp_pade, check_pade, &
!                      dom_rtp, read_function, zhat_pol_rtp)
!        CALL spec_static_resraman(nmodes, natom, zhat_pol_rtp, laser_in, freq, framecount_rtp, dom_rtp, dx, &
!                                  disp, check_pade, atom_mass_inv_sqrt)
!
!        DEALLOCATE (element, coord, mass_atom)
!        !***************************************************************************
!
!        !***************************************************************************
!    ELSEIF (read_function=='MD-RR') THEN
!        !CALL read_coord(natom,framecount,element,coord,rtp_dipole_x,periodic,mol_num,system,read_function,framecount_rtp,type_dipole)
!        ! CALL spec_resraman(natom,framecount,element,rtp_dipole_x,rtp_dipole_y,rtp_dipole_z,type_input,mol_num,system,&
!        !     read_function,dt,t_cor,pi,z_iso_resraman,z_aniso_resraman,dom,speed_light,const_planck,const_boltz,&
!        !    const_permit,temp,dom_rtp,laser_in_resraman,y_out)
    END IF

    CALL SYSTEM_CLOCK(count_1, count_rate, count_max) !Ending time
    time_final = count_1*1.0_dp/count_rate
    elapsed_time = time_final - time_init !Elapsed time

    WRITE (*, 1003) INT(elapsed_time), elapsed_time - INT(elapsed_time) !Write elapsed time
1003 FORMAT('  Wall Clock = ', i0, F0.9)

END PROGRAM vib2d
