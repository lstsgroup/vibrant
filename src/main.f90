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


!> @brief Main driver program for Vibrant.
!>
!> Controls the execution of all Vibrant workflows for vibrational and electronic
!> spectroscopy simulations, including IR, Raman, resonance Raman, and absorption
!> spectra based on static or molecular dynamics data.
!>
!> Initializes all data structures, parses input files and user options,
!> determines the calculation type from `read_function`, and gives
!> tasks to specialized routines for reading data, performing finite-difference
!> calculations, and generating spectra.
!>
!> After a completed process, all allocated arrays are deallocated and timing statistics are reported.
!>
!> @authors
!>   Ekin E. Winogradow, Johannes Scheffler, Moritz Leucke, Dorothea Golze
!>
!> @license
!>   Apache License 2.0
!>
PROGRAM vib2d

    USE, INTRINSIC           :: ISO_C_BINDING
    USE ISO_FORTRAN_ENV, ONLY: output_unit, error_unit
    USE kinds, ONLY: dp, str_len
    USE constants, ONLY: speed_light, const_planck, const_permit, pi, const_charge, const_boltz, joule_unit, &
                         debye, ev_unit, action_unit, bohr2ang, hartreebohr2evang, am_u, at_u, ang, fs2s, reccm2ev, &
                         hessian_factor, au2vm
    USE read_input, ONLY: parse_command_line, parse_input, check_input
    USE vib_types, ONLY: global_settings, systems, molecular_dynamics, static, dipoles, &
                         raman, init_global_settings, init_systems, init_molecular_dynamics, init_static, init_raman, deallocate_types
    USE setup, ONLY: masses_charges
    USE cell_types, ONLY: build_hmat, pbc, invert3x3, determinant3x3
    USE read_traj, ONLY: read_coord, read_coord_frame, read_normal_modes, read_static, read_static_resraman
    USE dipole_calc, ONLY: compute_dipole
    USE vel_cor, ONLY: cvv, cvv_iso, cvv_aniso, cvv_only_x, cvv_resraman
    USE fin_diff, ONLY: central_diff, forward_diff, finite_diff_static, finite_diff_static_resraman
    USE calc_spectra, ONLY: spec_power, normal_mode_analysis, spec_static_ir, spec_static_raman, &
                            spec_ir, spec_raman, spec_abs, spec_static_resraman, spec_abs_md, spec_resraman
    USE omp_lib, ONLY: omp_get_num_threads
    USE timing, ONLY: timings
    USE config_info, ONLY: output_config_info

    IMPLICIT NONE

    INCLUDE 'fftw3.f03'

    INTEGER                                         :: num_threads
    CHARACTER(LEN=str_len)                          :: input_file_name

    ! Variables of your derived types:
    TYPE(global_settings) :: gs
    TYPE(systems)        :: sys
    TYPE(molecular_dynamics)    :: md
    TYPE(static)            :: stats
    TYPE(dipoles)            :: dips
    TYPE(raman)            :: rams

    ! start timer for init
    CALL timings%register("initializing")

!$omp parallel
    num_threads = omp_get_num_threads()
!$omp end parallel

    !!Initialize
    CALL init_global_settings(gs)
    CALL init_systems(sys)
    CALL init_molecular_dynamics(md)
    CALL init_static(stats)
    CALL init_raman(rams)
    CALL output_config_info()

    !Parsing of the input
    CALL parse_command_line(input_file_name)
    CALL parse_input(gs, sys, md, stats, dips, rams, input_file_name)
    WRITE (*, '(90A, /)') REPEAT("-", 90)
    CALL check_input(gs, sys, md, stats, dips, rams)
    WRITE (*, '(90A)') REPEAT("-", 90)

    !***************************************************************************!
    !***************************************************************************!

    !!Power spectrum
    IF (gs%spectral_type%read_function=='P') THEN
        CALL timings%register("reading coordinates")
        CALL read_coord(sys%filename, gs, sys, dips)
        CALL timings%register("calculating charges")
        CALL masses_charges(sys)
        CALL timings%register("calculating power spectrum")
        CALL spec_power(gs, sys, md)
    !***************************************************************************!
    !***************************************************************************!

    !!MD-based IR spectrum
    ELSEIF (gs%spectral_type%read_function=='MD-IR') THEN
        CALL timings%register("reading coordinates")
        CALL read_coord(dips%dip_file, gs, sys, dips)
        IF (dips%type_dipole=='wannier') THEN !!fragment approach or whole supercell
            CALL masses_charges(sys)
        END IF

        CALL timings%register("calculating IR spectrum")
        CALL spec_ir(gs, sys, md, dips)
    !***************************************************************************!
    !***************************************************************************!

    !!MD-based Raman spectrum
    ELSEIF (gs%spectral_type%read_function=='MD-R') THEN
        !   sys%filename = wannier_free! <----  MUST BE ADJUSTED
        IF (dips%type_dipole=='berry' .OR. dips%type_dipole=='wannier') THEN
            CALL timings%register("reading coordinates")
            CALL read_coord(dips%dip_file, gs, sys, dips)
        ELSEIF (dips%type_dipole=='dfpt') THEN
            CALL timings%register("reading coordinates")
            CALL read_coord(dips%dip_x_file, gs, sys, dips)
        END IF
        CALL timings%register("calculating charges")
        CALL masses_charges(sys) !THIS FUNCTION IS NOT NEEDED HERE ?!

        CALL timings%register("calculating raman spectrum")
        CALL spec_raman(gs, sys, md, dips, rams)

    !***************************************************************************!
    !***************************************************************************!

    !!Normal mode analysis
    ELSEIF (gs%spectral_type%read_function=='NMA') THEN
        CALL timings%register("reading coordinates")
        CALL read_coord(sys%filename, gs, sys, dips)
        CALL timings%register("calculating charges")
        CALL masses_charges(sys)
        CALL timings%register("reading normal modes")
        CALL read_normal_modes(gs, sys, stats)
        CALL timings%register("normal mode analysis")
        CALL normal_mode_analysis(sys, stats, gs)
    !***************************************************************************!
    !***************************************************************************!

    !!Static IR spectrum
    ELSEIF (gs%spectral_type%read_function=='IR') THEN
        CALL timings%register("reading coordinates")
        CALL read_coord(sys%filename, gs, sys, dips)
        CALL timings%register("calculating charges")
        CALL masses_charges(sys)
        CALL timings%register("reading normal modes")
        CALL read_normal_modes(gs, sys, stats)
        CALL timings%register("reading static dipoles")
        CALL read_static(gs, sys, dips, rams)
        IF (stats%diag_hessian=='y') THEN
            CALL timings%register("normal mode analysis")
            CALL normal_mode_analysis(sys, stats, gs)
        END IF
        CALL timings%register("finite differences")
        CALL finite_diff_static(gs, sys, stats, dips, rams)

        CALL timings%register("calculating IR spectrum")
        CALL spec_static_ir(gs, sys, stats, dips)
    !***************************************************************************!
    !***************************************************************************!

    !!Static Raman spectrum
    ELSEIF (gs%spectral_type%read_function=='R') THEN
        CALL timings%register("reading coordinates")
        CALL read_coord(sys%filename, gs, sys, dips)
        CALL timings%register("calculating charges")
        CALL masses_charges(sys)
        CALL timings%register("reading normal modes")
        CALL read_normal_modes(gs, sys, stats)
        CALL timings%register("reading static dipoles")
        CALL read_static(gs, sys, dips, rams)
        IF (stats%diag_hessian=='y') THEN
            CALL timings%register("normal mode analysis")
            CALL normal_mode_analysis(sys, stats, gs)
        END IF
        CALL timings%register("finite differences")
        CALL finite_diff_static(gs, sys, stats, dips, rams)

        CALL timings%register("calculate Raman spectrum")
        CALL spec_static_raman(gs, sys, stats, dips, rams)
    
    !***************************************************************************!
    !***************************************************************************!

    !!Absorption spectrum
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
    !***************************************************************************!
    !***************************************************************************!

    !!Static resonance Raman spectrum
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
            CALL normal_mode_analysis(sys, stats, gs)
        END IF

        CALL timings%register("finite differences")
        CALL finite_diff_static_resraman(sys, rams)
        CALL timings%register("calculate absorption spectrum")
        CALL spec_abs(gs, sys, dips, rams)

        CALL timings%register("calculate resonance Raman spectrum")
        CALL spec_static_resraman(gs, sys, stats, rams)
    !***************************************************************************!
    !***************************************************************************!

    !!MD-averaged absorption spectrum
    ELSEIF (gs%spectral_type%read_function=='MD-ABS') THEN
        CALL timings%register("reading initial coordinates")
        CALL read_coord(dips%dip_x_file, gs, sys, dips, rams)
        CALL timings%register("read coordinates for all frames")
        CALL read_coord_frame(rams%RR%framecount_rtp, dips%dip_x_file, rams%RR%dip_x_rtp, sys)
        CALL read_coord_frame(rams%RR%framecount_rtp, dips%dip_y_file, rams%RR%dip_y_rtp, sys)
        CALL read_coord_frame(rams%RR%framecount_rtp, dips%dip_z_file, rams%RR%dip_z_rtp, sys)
        CALL timings%register("calculate MD-based absorption spectrum")
        CALL spec_abs_md(gs, sys, md, rams, dips)
    !***************************************************************************!
    !***************************************************************************!

    !!MD-based resonance Raman spectrum
    ELSEIF (gs%spectral_type%read_function=='MD-RR') THEN
        CALL timings%register("reading initial coordinates")
        CALL read_coord(dips%dip_x_file, gs, sys, dips, rams)
        CALL timings%register("read coordinates for all frames")
        CALL read_coord_frame(rams%RR%framecount_rtp, dips%dip_x_file, rams%RR%dip_x_rtp, sys)
        CALL read_coord_frame(rams%RR%framecount_rtp, dips%dip_y_file, rams%RR%dip_y_rtp, sys)
        CALL read_coord_frame(rams%RR%framecount_rtp, dips%dip_z_file, rams%RR%dip_z_rtp, sys)
        CALL timings%register("calculate MD-based absorption spectrum")
        CALL spec_abs_md(gs, sys, md, rams, dips)
        CALL timings%register("calculate resonance Raman spectrum")
        CALL spec_resraman(gs, sys, md, rams, dips)
    END IF

    !!Report all timings
    CALL timings%report_all()
    !!Deallocate arrays
    CALL deallocate_types(gs, sys, md, stats, rams, dips)

END PROGRAM vib2d
