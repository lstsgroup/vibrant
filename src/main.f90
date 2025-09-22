PROGRAM vib2d

    USE, INTRINSIC           :: ISO_C_BINDING
    USE kinds, ONLY: dp, str_len
    USE constants, ONLY: speed_light, const_planck, const_permit, pi, const_charge, const_boltz, joule_unit, &
                         debye, ev_unit, action_unit, bohr2ang, hartreebohr2evang, am_u, at_u, ang, fs2s, reccm2ev, &
                         hessian_factor, au2vm
    USE read_input, ONLY: parse_command_line, parse_input, check_input
    USE vib_types, ONLY: global_settings, systems, molecular_dynamics, static, dipoles, &
                         raman, init_global_settings, init_systems, init_molecular_dynamics, init_static, deallocate_types
    USE setup, ONLY: read_input, masses_charges, conversion, pbc_orthorombic, pbc_hexagonal
    USE read_traj, ONLY: read_coord, read_coord_frame, read_normal_modes, read_static, read_static_resraman
    USE dipole_calc, ONLY: center_mass, wannier_frag, solv_frag_index !wannier
    USE vel_cor, ONLY: cvv, cvv_iso, cvv_aniso, cvv_only_x, cvv_resraman
    USE fin_diff, ONLY: central_diff, forward_diff, finite_diff_static, finite_diff_static_resraman
    USE calc_spectra, ONLY: spec_power, normal_mode_analysis, spec_static_ir, spec_static_raman, &
                            spec_ir, spec_raman, spec_abs, spec_static_resraman!, spec_resraman
    USE omp_lib, ONLY: omp_get_num_threads
    USE timing, ONLY: timings
    USE config_info, ONLY: output_config_info

    IMPLICIT NONE

    INCLUDE 'fftw3.f03'

    INTEGER                                         :: num_threads
    ! Variables of your derived types:
    TYPE(global_settings) :: gs
    TYPE(systems)        :: sys
    TYPE(molecular_dynamics)    :: md

    TYPE(static)            :: stats
    TYPE(dipoles)            :: dips
    TYPE(raman)            :: rams
    CHARACTER(LEN=str_len)          :: input_file_name

    ! start timer for init
    CALL timings%register("initializing")

!$omp parallel
    num_threads = omp_get_num_threads()
!$omp end parallel

    CALL init_global_settings(gs)
    CALL init_systems(sys)
    CALL init_molecular_dynamics(md)
    CALL init_static(stats)

    CALL output_config_info()

    !***************************************************************************
    !                                                                        ***
    CALL parse_command_line(input_file_name)

    CALL parse_input(gs, sys, md, stats, dips, rams, input_file_name)
    PRINT *, REPEAT('-', 30)
    CALL check_input(gs, sys, md, stats, dips, rams)
    PRINT *, REPEAT('-', 30)
    !                                                                        ***
    !***************************************************************************

    CALL conversion(md%dt, md%freq_range, rams%RR%dt_rtp, rams%RR%freq_range_rtp)
!!
!!    !***************************************************************************
    IF (gs%spectral_type%read_function=='P') THEN
        CALL timings%register("reading coordinates")
        CALL read_coord(sys%filename, gs, sys)
        CALL timings%register("calculating charges")
        CALL masses_charges(sys)
        CALL timings%register("calculating power spectrum")
        CALL spec_power(gs, sys, md)
!        !***************************************************************************
!        !***************************************************************************
    ELSEIF (gs%spectral_type%read_function=='MD-IR') THEN
        CALL timings%register("reading coordinates")
        CALL read_coord(dips%dip_file, gs, sys, dips)
        !  IF (sys%system=='1' .OR. sys%system=='2' .AND. dips%type_dipole=='wannier') THEN !!fragment approach or whole supercell
        !      IF (sys%cell%cell_type=='1' .OR. sys%cell%cell_type=='2') THEN !!KP or SC
        !          CALL masses_charges(sys)
        !      END IF
        !  END IF

        CALL timings%register("calculating IR spectrum")
        CALL spec_ir(gs, sys, md, dips)
!        !***************************************************************************
!
!        !***************************************************************************
    ELSEIF (gs%spectral_type%read_function=='MD-R') THEN
        !   sys%filename = wannier_free! <----  MUST BE ADJUSTED
        IF (dips%type_dipole=='berry') THEN
            CALL timings%register("reading coordinates")
            CALL read_coord(dips%dip_file, gs, sys, dips)
        ELSEIF (dips%type_dipole=='dfpt') THEN
            CALL timings%register("reading coordinates")
            CALL read_coord(dips%dip_x_file, gs, sys, dips)
        END IF
        CALL timings%register("calculating charges")
        CALL masses_charges(sys)

        CALL timings%register("calculating raman spectrum")
        CALL spec_raman(gs, sys, md, dips, rams)
        !***************************************************************************

        !***************************************************************************
    ELSEIF (gs%spectral_type%read_function=='NMA') THEN
        CALL timings%register("reading coordinates")
        CALL read_coord(sys%filename, gs, sys, dips)
        CALL timings%register("calculating charges")
        CALL masses_charges(sys)
        CALL timings%register("reading normal modes")
        CALL read_normal_modes(gs, sys, stats)
        CALL timings%register("normal mode analysis")
        CALL normal_mode_analysis(sys, stats)
!        !***************************************************************************
!
!        !***************************************************************************
    ELSEIF (gs%spectral_type%read_function=='IR') THEN
        CALL timings%register("reading coordinates")
        CALL read_coord(sys%filename, gs, sys, dips)
        CALL timings%register("calculating charges")
        CALL masses_charges(sys)
        CALL timings%register("reading normal modes")
        CALL read_normal_modes(gs, sys, stats)
        CALL timings%register("reading static dipoles")
        CALL read_static(sys, dips, rams)
        IF (stats%diag_hessian=='y') THEN
            CALL timings%register("normal mode analysis")
            CALL normal_mode_analysis(sys, stats)
        END IF
        CALL timings%register("finite differences")
        CALL finite_diff_static(gs, sys, stats, dips, rams)

        CALL timings%register("calculating IR spectrum")
        CALL spec_static_ir(gs, stats, dips)
!        !***************************************************************************
!
!        !***************************************************************************
    ELSEIF (gs%spectral_type%read_function=='R') THEN
        CALL timings%register("reading coordinates")
        CALL read_coord(sys%filename, gs, sys, dips)
        CALL timings%register("calculating charges")
        CALL masses_charges(sys)
        CALL timings%register("reading normal modes")
        CALL read_normal_modes(gs, sys, stats)
        CALL timings%register("reading static dipoles")
        CALL read_static(sys, dips, rams)
        !  IF (type_dipole=='2') THEN
        !     CALL read_static(static_dip_x_file, static_dip_x, gs, sys, rams)
        !    CALL read_static(static_dip_y_file, static_dip_y, gs, sys, rams)
        !   CALL read_static(static_dip_z_file, static_dip_z, gs, sys, rams)
        ! END IF
        IF (stats%diag_hessian=='y') THEN
            CALL timings%register("normal mode analysis")
            CALL normal_mode_analysis(sys, stats)
        END IF
        CALL timings%register("finite differences")
        CALL finite_diff_static(gs, sys, stats, dips, rams)

        CALL timings%register("calculate Raman spectrum")
        CALL spec_static_raman(gs, sys, stats, rams)
!        !***************************************************************************
!
!        !***************************************************************************
    ELSEIF (gs%spectral_type%read_function=='ABS') THEN
        CALL timings%register("reading coordinates")
        CALL read_coord(sys%filename, gs, sys, dips)
        CALL timings%register("reading dipoles")
        CALL read_static_resraman(dips%dip_x_file, rams%RR%static_dip_x_rtp, sys, rams)
        CALL read_static_resraman(dips%dip_y_file, rams%RR%static_dip_y_rtp, sys, rams)
        CALL read_static_resraman(dips%dip_z_file, rams%RR%static_dip_z_rtp, sys, rams)
        CALL timings%register("finite differences")
        CALL finite_diff_static_resraman(sys, rams) !<-- CHANGE ?

        CALL timings%register("calculate absorption spectrum")
        CALL spec_abs(gs, sys, dips, rams)
        !***************************************************************************

        !***************************************************************************
    ELSEIF (gs%spectral_type%read_function=='RR') THEN
        CALL timings%register("reading coordinates")
        CALL read_coord(sys%filename, gs, sys, dips)
        CALL timings%register("calculating charges")
        CALL masses_charges(sys)
        CALL timings%register("reading normal modes")
        CALL read_normal_modes(gs, sys, stats)
        CALL timings%register("reading dipoles")
        CALL read_static_resraman(dips%dip_x_file, rams%RR%static_dip_x_rtp, sys, rams)
        CALL read_static_resraman(dips%dip_y_file, rams%RR%static_dip_y_rtp, sys, rams)
        CALL read_static_resraman(dips%dip_z_file, rams%RR%static_dip_z_rtp, sys, rams)

        IF (stats%diag_hessian=='y') THEN
            CALL timings%register("normal mode analysis")
            CALL normal_mode_analysis(sys, stats)
        END IF

        CALL timings%register("finite differences")
        CALL finite_diff_static_resraman(sys, rams)
        CALL timings%register("calculate absorption spectrum")
        CALL spec_abs(gs, sys, dips, rams)

        CALL timings%register("calculate resonance Raman spectrum")
        CALL spec_static_resraman(gs, sys, stats, rams)
        !***************************************************************************
        !***************************************************************************
!    ELSEIF (gs%spectral_type%read_function=='MD-RR') THEN
!        !sys%filename = rtp_dipole_x ! <----  MUST BE ADJUSTED
!        !CALL read_coord(gs, sys, dips)
!        !natom = sys%natom
!        !framecount = sys%framecount
!        !mol_num = sys%mol_num
!        !CALL spec_resraman(natom,framecount,element,rtp_dipole_x,rtp_dipole_y,rtp_dipole_z,type_input,mol_num,system,&
!        !     read_function,dt,z_iso_resraman,z_aniso_resraman,freq_range,freq_range_rtp,laser_in_resraman,y_out)
    END IF

    CALL timings%report_all()

    CALL deallocate_types(sys, md, stats, rams, dips)

END PROGRAM vib2d
