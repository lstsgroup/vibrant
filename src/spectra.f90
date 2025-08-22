MODULE calc_spectra

   USE setup, ONLY: read_input, conversion
   USE kinds, ONLY: dp
   USE vib_types, ONLY: global_settings, systems, molecular_dynamics, static, dipoles, raman

   USE constants, ONLY: pi, fs2s, debye, speed_light, const_planck, const_boltz, const_permit, cm2m, a3_to_debye_per_e, &
                        hartreebohr2evang, hessian_factor, bohr2ang, reccm2ev, am_u, debye2cm, avo_num, au2vm, ang, at_u, &
                        speed_light_au, debye, reccm2au
   USE read_traj, ONLY: read_coord_frame
   USE fin_diff, ONLY: central_diff, forward_diff
   USE vel_cor, ONLY: cvv, cvv_iso, cvv_aniso, cvv_only_x, cvv_resraman
   USE dipole_calc, ONLY: center_mass, wannier_frag, wannier
   USE pade, ONLY: interpolate

   USE, INTRINSIC                              :: ISO_C_BINDING
   USE OMP_LIB

   IMPLICIT NONE

   INCLUDE 'fftw3.f03'

   PUBLIC :: spec_power, normal_mode_analysis, spec_static_ir, spec_static_raman, spec_ir, spec_raman, spec_abs, spec_static_resraman!, spec_resraman

CONTAINS
   SUBROUTINE spec_power(gs, sys, md)

      TYPE(global_settings), INTENT(INOUT)        :: gs
      TYPE(systems), INTENT(INOUT)                :: sys
      TYPE(molecular_dynamics), INTENT(INOUT)     :: md

      INTEGER                                                  :: stat, i
      INTEGER(kind=dp)                                          :: plan
      REAL(kind=dp)                               :: freq_range, freq_res, power_const
      REAL(kind=dp), DIMENSION(:), ALLOCATABLE      :: power_int, freq

      ALLOCATE (md%zhat(0:md%t_cor*2 - 1), power_int(0:md%t_cor*2 - 1), freq(0:md%t_cor*2 - 1))

      md%zhat = COMPLEX(0._dp, 0.0_dp)

      CALL read_coord_frame(sys%natom, sys%filename, md%coord_v, sys)

      IF (sys%type_traj == 'pos') THEN   !!If it is from positions, do finite differences first
         CALL central_diff(sys%natom, md%coord_v, md%v, sys, md)

         CALL cvv(sys%natom, md%v, sys, gs, md)

      ELSEIF (sys%type_traj == 'vel') THEN   !!If it is from velocities, compute autocorrelation directly

         CALL cvv(sys%natom, md%coord_v, sys, gs, md)

      END IF

      freq_range = REAL((1.0_dp/(md%dt*fs2s))/speed_light, kind=dp)
      freq_res = REAL(freq_range/(2.0_dp*md%t_cor), kind=dp)

      CALL dfftw_plan_dft_r2c_1d(plan, 2*md%t_cor, md%z, md%zhat, FFTW_ESTIMATE) !!!FFT
      CALL dfftw_execute_dft_r2c(plan, md%z, md%zhat)
      CALL dfftw_destroy_plan(plan)

      md%zhat = REAL(md%zhat, kind=dp)
      !!Unit conversion to K.cm
      power_const = (md%dt*fs2s*am_u*speed_light/const_boltz)*2.0_dp/(sys%natom*3.0_dp)

      OPEN (UNIT=63, FILE='power_spec.txt', STATUS='unknown', IOSTAT=stat) !!write the output
      DO i = 0, 2*md%t_cor - 1
         freq(i) = i*freq_res
         power_int(i) = md%zhat(i)*power_const
         IF (freq(i) .GE. 5000_dp) CYCLE
         WRITE (63, *) freq(i), power_int(i)
      END DO

      CLOSE (63)
      DEALLOCATE (power_int, freq)

   END SUBROUTINE spec_power
!***********************************************************************************************!
!***********************************************************************************************!
   SUBROUTINE spec_ir(gs, sys, md, dips)
      TYPE(global_settings), INTENT(INOUT)        :: gs
      TYPE(systems), INTENT(INOUT)                :: sys
      TYPE(molecular_dynamics), INTENT(INOUT)     :: md
      TYPE(dipoles), INTENT(INOUT)     :: dips

      INTEGER                                                  :: stat, i
      INTEGER(kind=dp)                                          :: plan
      REAL(kind=dp)                                          :: freq_range, freq_res, sinc_const, ir_const
      REAL(kind=dp), DIMENSION(:), ALLOCATABLE               :: ir_int, freq

        IF (dips%type_dipole=='wannier') THEN !!fragment approach or whole supercell
      !      IF (sys%cell%cell_type=='1' .OR. sys%cell%cell_type=='2') THEN !!KP or SC
                CALL read_coord_frame(sys%natom, dips%dip_file, md%coord_v, sys)
                print*, sys%cell%box_x, dips%dip_file
                CALL center_mass(dips%dip_file, sys%fragments%fragment, gs, sys, md, dips)
                CALL wannier_frag(sys%fragments%natom_frag, sys%filename, dips%dipole, sys%fragments%fragment, gs, sys, md, dips)
      !      ELSEIF (sys%cell%cell_type=='3') THEN !!SC with solvent
      !          CALL read_coord_frame(sys%natom, sys%filename, md%coord_v, sys)
      !          CALL solv_frag_index(sys%filename, sys%fragments%natom_frag, sys%fragments%fragment, sys, md, dips)
      !          CALL wannier_frag(sys%fragments%natom_frag, sys%filename, dips%dipole, sys%fragments%fragment, gs, sys, md, dips)
      !       END IF
    END IF

      ALLOCATE (md%zhat(0:2*md%t_cor - 1), ir_int(0:2*md%t_cor - 1), freq(0:2*md%t_cor - 1))

      md%zhat = COMPLEX(0._dp, 0.0_dp)

      !   IF (sys%system=='2' .AND. dips%type_dipole=='wannier') THEN
      !       sys%mol_num = 1
      !   END IF
!    IF (sys%system=='1' .OR. (sys%system=='2' .AND. dips%type_dipole=='wannier')) THEN  !!fragment approach or the whole cell
      !       CALL central_diff(sys%mol_num, dips%dipole, md%v, sys, md)
      !       CALL cvv(sys%fragments%nfrag, md%v, sys, md)
     
   ! IF (dips%type_dipole == 'berry') THEN !!Berry phase dipoles
      !   CALL read_coord_frame(sys%mol_num, dips%dip_file, md%coord_v, sys)
     !    CALL central_diff(sys%mol_num, md%coord_v, md%v, sys, md)
         CALL central_diff(sys%mol_num, dips%dipole, md%v, sys, md)
         CALL cvv(sys%mol_num, md%v, sys, gs, md)
    !  END IF

      CALL dfftw_plan_dft_r2c_1d(plan, 2*md%t_cor, md%z(0:2*md%t_cor - 1), md%zhat(0:2*md%t_cor - 1), FFTW_ESTIMATE) !!FFT!!
      CALL dfftw_execute_dft_r2c(plan, md%z, md%zhat)
      CALL dfftw_destroy_plan(plan)

      freq_range = REAL((1.0_dp/(md%dt*fs2s))/speed_light, kind=dp)
      freq_res = REAL(freq_range/(2.0_dp*md%t_cor), kind=dp)
      !!sinc function
      sinc_const = freq_res*md%dt*fs2s*2.0_dp*pi*speed_light

      md%zhat = REAL(md%zhat, kind=dp)
      !!conversion of debye^2 to (C*m)^2
      md%zhat(:) = md%zhat(:)*debye2cm*debye2cm !!conversion of debye^2 to (C*m)^2
      !!conversion of IR units to K*cm*km*mol^-1
      ir_const = avo_num*md%dt*fs2s*2.0_dp*10.0d0/(12.0_dp*const_permit*speed_light*const_boltz)

      OPEN (UNIT=61, FILE='IR_spectrum.txt', STATUS='unknown', IOSTAT=stat) !!write output
      DO i = 0, 2*md%t_cor - 1
         freq(i) = i*freq_res
         ir_int(i) = md%zhat(i)*ir_const*(sinc_const*(i)/SIN(sinc_const*(i)))**2._dp
         IF (freq(i) .GE. 5000_dp) CYCLE
         ir_int(0) = 0.00_dp
         WRITE (61, *) freq(i), -1.0_dp*REAL(ir_int(i), kind=dp)
      END DO
      CLOSE (61)

      DEALLOCATE (ir_int, freq)

   END SUBROUTINE spec_ir

!!!****************************************************************************************!
!!!****************************************************************************************!
   SUBROUTINE spec_raman(gs, sys, md, dips, rams)
      TYPE(global_settings), INTENT(INOUT)        :: gs
      TYPE(systems), INTENT(INOUT)        :: sys
      TYPE(molecular_dynamics), INTENT(INOUT)        :: md
      TYPE(dipoles), INTENT(INOUT)        :: dips
      TYPE(raman), INTENT(INOUT)        :: rams

      INTEGER                                                  :: stat, i, j
      INTEGER(kind=dp)                                          :: plan
      INTEGER, DIMENSION(:), ALLOCATABLE                         :: natom_frag_x, natom_frag_free
      INTEGER, DIMENSION(:), ALLOCATABLE                         :: natom_frag_y, natom_frag_z
      INTEGER, DIMENSION(:, :, :), ALLOCATABLE                     :: fragment_x, fragment_free
      INTEGER, DIMENSION(:, :, :), ALLOCATABLE                     :: fragment_y, fragment_z
      COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE                 :: zhat_iso, zhat_aniso
      REAL(kind=dp), DIMENSION(:), ALLOCATABLE                    :: zhat_unpol_x, zhat_depol_x, zhat_para_all, zhat_depol
      REAL(kind=dp), DIMENSION(:), ALLOCATABLE                    :: raman_ortho, raman_para, raman_depol, raman_unpol
      REAL(kind=dp), DIMENSION(:), ALLOCATABLE                    :: raman_const, sinc_func, freq
      REAL(kind=dp)                                             :: f, freq_res
      REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE                :: dip_free, dip_x, dip_y, dip_z
      REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE                :: alpha_x, alpha_y, alpha_z
      REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE                :: alpha_diff_x, alpha_diff_y, alpha_diff_z

!!!!ALLOCATION!!!
      ALLOCATE (alpha_x(sys%framecount, sys%mol_num, 3))
      ALLOCATE (alpha_y(sys%framecount, sys%mol_num, 3))
      ALLOCATE (alpha_z(sys%framecount, sys%mol_num, 3))

!    IF (rams%averaging=='1') THEN

!!FIELD_FREE!!!
      !       IF (sys%system=='1' .OR. dips%type_dipole=='wannier') THEN
      !           IF (sys%cell%cell_type.NE.'3') THEN
      !               !CALL read_coord_frame(sys, md)
      !              CALL read_coord_frame(sys%natom, rams%wannier_free, md%coord_v, sys) !<- this needs to be ADJUSTED
      !              ! filename different and fragment_free different
      !              CALL center_mass(rams%wannier_free, fragment_free, gs, sys, md, dips)
      !              CALL wannier_frag(sys%fragments%natom_frag, rams%wannier_free, dip_free, fragment_free, gs, sys, md, dips)
      !          ELSEIF (sys%cell%cell_type=='3') THEN
      !              CALL read_coord_frame(sys%natom, rams%wannier_free, md%coord_v, sys)
      !              CALL solv_frag_index(rams%wannier_free, natom_frag_free, fragment_free, sys, md, dips)
!
      !             CALL wannier_frag(natom_frag_free, rams%wannier_free, dip_free, fragment_free, gs, sys, md, dips)
      !         END IF
      !     ELSEIF (sys%system=='2') THEN
      IF (dips%type_dipole == 'berry') THEN
         CALL read_coord_frame(sys%mol_num, dips%dip_file, dip_free, sys)
      END IF
      !     END IF

!!!X-FIELD!!!
      CALL read_coord_frame(sys%natom, dips%dip_x_file, md%coord_v, sys)
      !     IF (sys%system=='1' .OR. dips%type_dipole=='wannier') THEN
      !         IF (sys%cell%cell_type.NE.'3') THEN
      !             CALL center_mass(rams%wannier_x, fragment_x, gs, sys, md, dips)
      !             CALL wannier_frag(sys%fragments%natom_frag, rams%wannier_x, dip_x, fragment_x, gs, sys, md, dips)

      !          ELSEIF (sys%cell%cell_type=='3') THEN
      !              CALL solv_frag_index(rams%wannier_x, natom_frag_x, fragment_x, sys, md, dips)
      !              CALL wannier_frag(natom_frag_x, rams%wannier_x, dip_x, fragment_x, gs, sys, md, dips)
      !          END IF
      !         IF (sys%system=='1') THEN
      !             CALL forward_diff(sys%mol_num, alpha_x, dip_free, dip_x, gs, sys)
      !          ELSEIF (sys%system=='2' .AND. dips%type_dipole=='wannier') THEN
      !              CALL forward_diff(sys%fragments%nfrag, alpha_x, dip_free, dip_x, gs, sys)
      !          END IF
      !      ELSEIF (sys%system=='2') THEN

      IF (dips%type_dipole == 'berry') THEN
         CALL forward_diff(sys%mol_num, alpha_x, dip_free, md%coord_v, gs, sys, dips)

      ELSEIF (dips%type_dipole == 'dfpt') THEN
         !alpha_x = REAL(md%coord_v*((8.988d+15)/(5.142d+11*3.33564d-30)), kind=dp) !conversion to debye/E
         alpha_x = REAL(md%coord_v*a3_to_debye_per_e, kind=dp) !conversion of A^3 to debye/E
      END IF
      ! END IF

      ! IF (sys%system=='2' .AND. dips%type_dipole=='wannier') THEN
      !     CALL central_diff(sys%fragments%nfrag, alpha_x, alpha_diff_x, sys, md)
      ! ELSE
      CALL central_diff(sys%mol_num, alpha_x, alpha_diff_x, sys, md)
      ! END IF

!!!Y-FIELD!!!

      CALL read_coord_frame(sys%natom, dips%dip_y_file, md%coord_v, sys)
      !  IF (sys%system=='1' .OR. dips%type_dipole=='wannier') THEN
      !      IF (sys%cell%cell_type.NE.'3') THEN

      !         CALL center_mass(rams%wannier_y, fragment_y, gs, sys, md, dips)
      !         CALL wannier_frag(sys%fragments%natom_frag, rams%wannier_y, dip_y, fragment_y, gs, sys, md, dips)
      !     ELSEIF (sys%cell%cell_type=='3') THEN

      !        CALL solv_frag_index(rams%wannier_y, natom_frag_y, fragment_y, sys, md, dips)
      !        CALL wannier_frag(natom_frag_y, rams%wannier_y, dip_y, fragment_y, gs, sys, md, dips)
      !    END IF
      !    IF (sys%system=='1') THEN
      !        CALL forward_diff(sys%mol_num, alpha_y, dip_free, dip_y, gs, sys)
      !    ELSEIF (sys%system=='2' .AND. dips%type_dipole=='wannier') THEN
      !        CALL forward_diff(sys%fragments%nfrag, alpha_y, dip_free, dip_y, gs, sys)
      !    END IF
      ! ELSEIF (sys%system=='2') THEN
      IF (dips%type_dipole == 'berry') THEN
         CALL forward_diff(sys%mol_num, alpha_y, dip_free, md%coord_v, gs, sys, dips)

      ELSEIF (dips%type_dipole == 'dfpt') THEN
         alpha_y = REAL(md%coord_v*a3_to_debye_per_e, kind=dp) !conversion of A^3 to debye/E
         !alpha_y = REAL(md%coord_v*((8.988d+15)/(5.142d+11*3.33564d-30)), kind=dp) !conversion to debye/E
      END IF
      ! END IF

      !  IF (sys%system=='2' .AND. dips%type_dipole=='wannier') THEN
      !      CALL central_diff(sys%fragments%nfrag, alpha_y, alpha_diff_y, sys, md)
      !  ELSE
      CALL central_diff(sys%mol_num, alpha_y, alpha_diff_y, sys, md)
      !  END IF

!!!Z-FIELD!!!

      CALL read_coord_frame(sys%natom, dips%dip_z_file, md%coord_v, sys)
      !   IF (sys%system=='1' .OR. dips%type_dipole=='wannier') THEN
      !       IF (sys%cell%cell_type.NE.'3') THEN

!                CALL center_mass(rams%wannier_z, fragment_z, gs, sys, md, dips)
!                CALL wannier_frag(sys%fragments%natom_frag, rams%wannier_z, dip_z, fragment_z, gs, sys, md, dips)
!            ELSEIF (sys%cell%cell_type=='3') THEN

      !               CALL solv_frag_index(rams%wannier_z, natom_frag_z, fragment_z, sys, md, dips)
      !               CALL wannier_frag(natom_frag_z, rams%wannier_z, dip_z, fragment_z, gs, sys, md, dips)
      !          END IF
      !          IF (sys%system=='1') THEN
      !              CALL forward_diff(sys%mol_num, alpha_z, dip_free, dip_z, gs, sys)
      !          ELSEIF (sys%system=='2' .AND. dips%type_dipole=='wannier') THEN
      !              CALL forward_diff(sys%fragments%nfrag, alpha_z, dip_free, dip_z, gs, sys)
      !          END IF

      !     ELSEIF (sys%system=='2') THEN

      IF (dips%type_dipole == 'berry') THEN
         CALL forward_diff(sys%mol_num, alpha_z, dip_free, md%coord_v, gs, sys, dips)

      ELSEIF (dips%type_dipole == 'dfpt') THEN
         alpha_z = REAL(md%coord_v*a3_to_debye_per_e, kind=dp) !conversion of A^3 to debye/E
      END IF
      ! END IF

      ! IF (sys%system=='2' .AND. dips%type_dipole=='wannier') THEN
      !     CALL central_diff(sys%fragments%nfrag, alpha_z, alpha_diff_z, sys, md)
      ! ELSE
      CALL central_diff(sys%mol_num, alpha_z, alpha_diff_z, sys, md)
      ! END IF

!!!ACF AND FFT CALC!!!
      !   PRINT *, sys%fragments%nfrag, 'sys%fragments%nfrag check'
      ALLOCATE (zhat_iso(0:md%t_cor*2), zhat_aniso(0:md%t_cor*2))
      ALLOCATE (raman_ortho(0:md%t_cor*2), raman_para(0:md%t_cor*2), raman_depol(0:md%t_cor*2), raman_unpol(0:md%t_cor*2))
      ALLOCATE (freq(0:md%t_cor*2), raman_const(0:md%t_cor*2))

      zhat_iso = COMPLEX(0._dp, 0.0_dp)
      zhat_aniso = COMPLEX(0._dp, 0.0_dp)
      ! raman_unpol = COMPLEX(0._dp, 0.0_dp)
      ! raman_depol = COMPLEX(0._dp, 0.0_dp)

      !IF (sys%system=='1' .OR. (sys%system=='2' .AND. dips%type_dipole=='wannier')) THEN
      !    CALL cvv_iso(sys%fragments%nfrag, rams%z_iso, alpha_diff_x, alpha_diff_y, alpha_diff_z, sys, md)
      !  ELSE
      CALL cvv_iso(sys%mol_num, rams%z_iso, alpha_diff_x, alpha_diff_y, alpha_diff_z, sys, md)
      ! END IF

      CALL dfftw_plan_dft_r2c_1d(plan, 2*md%t_cor, rams%z_iso, zhat_iso, FFTW_ESTIMATE)
      CALL dfftw_execute_dft_r2c(plan, rams%z_iso, zhat_iso)
      CALL dfftw_destroy_plan(plan)

      !   IF (sys%system=='1' .OR. (sys%system=='2' .AND. dips%type_dipole=='wannier')) THEN
      !       CALL cvv_aniso(sys%fragments%nfrag, rams%z_aniso, alpha_diff_x, alpha_diff_y, alpha_diff_z, sys, md)
      !  ELSE
      CALL cvv_aniso(sys%mol_num, rams%z_aniso, alpha_diff_x, alpha_diff_y, alpha_diff_z, sys, md)
      !  END IF

      CALL dfftw_plan_dft_r2c_1d(plan, 2*md%t_cor, rams%z_aniso, zhat_aniso, FFTW_ESTIMATE)
      CALL dfftw_execute_dft_r2c(plan, rams%z_aniso, zhat_aniso)
      CALL dfftw_destroy_plan(plan)

      !!!Find the frequency resolution in cm^-1
      freq_res = REAL(md%freq_range/(2.0_dp*md%t_cor), kind=dp)
      f = freq_res*md%dt*1.883652d-4

      zhat_iso = REAL(zhat_iso, kind=dp)
      zhat_aniso = REAL(zhat_aniso, kind=dp)
      raman_ortho = 0.0d0
      raman_para = 0.0d0
      raman_unpol = 0.0d0
      raman_depol = 0.0d0

      !!Unit conversion of Debye^2/(E^2*s^2) into C^4*s^2/kg^2
      zhat_iso(:) = zhat_iso(:)*debye2cm*debye2cm/(au2vm*au2vm)
      zhat_aniso(:) = zhat_aniso(:)*debye2cm*debye2cm/(au2vm*au2vm)

      DO i = 0, 2*md%t_cor - 2
         freq(i) = i*freq_res
         !!conversion of the Raman intensities into m^2*K*cm*10^-30!!
         raman_const(i) = const_planck/(8.0_dp*const_boltz*const_permit*const_permit) &
                          *1.e+30*md%dt*fs2s*((((rams%laser_in/reccm2ev - freq(i))/cm2m)**4)/freq(i))* &
                          (1.0_dp/(1.0_dp - EXP(-1._dp*const_planck*speed_light*cm2m*freq(i)/ &
                                                (const_boltz*gs%temp))))*2.0_dp
      END DO

!!!ORTHOGONAL!!!
      OPEN (UNIT=63, FILE='raman_orthogonal.txt', STATUS='unknown', IOSTAT=stat)
      DO i = 0, 2*md%t_cor - 2
         zhat_aniso(i + 1) = REAL(zhat_aniso(i + 1), kind=dp)*(f*(i + 1)/SIN(f*(i + 1)))**2._dp
         raman_ortho(i) = ((REAL(zhat_aniso(i), kind=dp))/15.0_dp)*raman_const(i)
         raman_ortho(0) = 0.0_dp
         IF (freq(i) .GE. 5000.0_dp) CYCLE
         WRITE (63, *) freq(i), raman_ortho(i)
      END DO
      CLOSE (63)

!!!PARALLEL!!!
      OPEN (UNIT=64, FILE='raman_parallel.txt', STATUS='unknown', IOSTAT=stat)
      DO i = 0, 2*md%t_cor - 2
         zhat_iso(i + 1) = REAL(zhat_iso(i + 1), kind=dp)*(f*(i + 1)/SIN(f*(i + 1)))**2._dp
         raman_para(i) = (zhat_iso(i) + (zhat_aniso(i)*4.0_dp/45.0_dp))*raman_const(i)
         raman_para(0) = 0.0_dp
         IF (freq(i) .GE. 5000.0_dp) CYCLE
         WRITE (64, *) freq(i), raman_para(i)
      END DO
      CLOSE (64)

!!!UNPOL!!!
      OPEN (UNIT=65, FILE='raman_unpolarized.txt', STATUS='unknown', IOSTAT=stat)

      DO i = 0, 2*md%t_cor - 2
         raman_unpol(i) = raman_ortho(i) + raman_para(i)
         raman_unpol(0) = 0.00_dp
         IF (freq(i) .GE. 5000.0_dp) CYCLE
         WRITE (65, *) freq(i), raman_unpol(i)
      END DO
      CLOSE (65)

!!!DEPOL RATIO!!!
      OPEN (UNIT=66, FILE='raman_depolarization_ratio.txt', STATUS='unknown', IOSTAT=stat)

      DO i = 0, 2*md%t_cor - 2
         raman_depol(i) = REAL(raman_ortho(i), kind=dp)/REAL(raman_para(i), kind=dp)
         IF (freq(i) .GE. 5000.0_dp) CYCLE
         WRITE (66, *) freq(i), raman_depol(i)
      END DO

      CLOSE (66)
      DEALLOCATE (rams%z_iso, rams%z_aniso)
      DEALLOCATE (raman_depol, raman_para, raman_unpol, raman_ortho, zhat_aniso, zhat_iso)

!    ELSEIF (rams%averaging=='2') THEN

!        IF (dips%type_dipole=='wannier') THEN
!            CALL read_coord_frame(sys%natom, rams%wannier_free, md%coord_v, sys)
      !sys%filename = rams%wannier_free
      !CALL read_coord_frame(sys, md)
!            CALL wannier(sys%filename, dip_free, sys, md)
!        ELSEIF (dips%type_dipole=='berry') THEN

      !           CALL read_coord_frame(sys%natom, rams%wannier_free, dip_free, sys)
      !       END IF

!!!X-FIELD!!!
      !      IF (rams%direction=='1') THEN
!!!X-FIELD!!!
      !    sys%filename = rams%wannier_x

      !    CALL read_coord_frame(sys%natom, rams%wannier_x, md%coord_v, sys)
      !    IF (dips%type_dipole=='wannier') THEN
      !        CALL wannier(rams%wannier_x, dip_x, sys, md)
      !       CALL forward_diff(sys%mol_num, alpha_x, dip_free, dip_x, gs, sys)
      !   ELSEIF (dips%type_dipole=='berry') THEN
      !       CALL forward_diff(sys%mol_num, alpha_x, dip_free, md%coord_v, gs, sys)
      !   END IF
      !   CALL central_diff(sys%natom, alpha_x, alpha_diff_x, sys, md)

!!!Y-FIELD!!!
      !ELSEIF (rams%direction=='2') THEN

      !   CALL read_coord_frame(sys%natom, rams%wannier_y, md%coord_v, sys)
      !   IF (dips%type_dipole=='wannier') THEN
      !      CALL wannier(rams%wannier_y, dip_y, sys, md)
      !     CALL forward_diff(sys%mol_num, alpha_y, dip_free, dip_y, gs, sys)
      ! ELSEIF (dips%type_dipole=='berry') THEN
      !    CALL forward_diff(sys%mol_num, alpha_y, dip_free, md%coord_v, gs, sys)
      ! END IF
      ! CALL central_diff(sys%natom, alpha_y, alpha_diff_y, sys, md)

      ! ELSEIF (rams%direction=='3') THEN
!!!Z-FIELD!!!
      !  CALL read_coord_frame(sys%natom, rams%wannier_z, md%coord_v, sys)
      !   IF (dips%type_dipole=='wannier') THEN
      !        CALL wannier(rams%wannier_z, dip_z, sys, md)
      !         CALL forward_diff(sys%mol_num, alpha_z, dip_free, dip_z, gs, sys)
      !      ELSEIF (dips%type_dipole=='berry') THEN
      !           CALL forward_diff(sys%mol_num, alpha_z, dip_free, md%coord_v, gs, sys)
      !        END IF
      !        CALL central_diff(sys%natom, alpha_z, alpha_diff_z, sys, md)
      !CALL central_diff(md%dt, sys%natom, sys%framecount, alpha_z, alpha_diff_z, gs%spectral_type%read_function, sys%mol_num, sys%system)

      ! END IF

      !ALLOCATE (zhat_para(0:md%t_cor*2), zhat_unpol_x(0:md%t_cor*2), zhat_ortho(0:md%t_cor*2), zhat_depol_x(0:md%t_cor*2))

      ! zhat_para = COMPLEX(0._dp, 0.0_dp)
      ! zhat_ortho = COMPLEX(0._dp, 0.0_dp)
      ! zhat_unpol_x = COMPLEX(0._dp, 0.0_dp)
      ! zhat_depol_x = COMPLEX(0._dp, 0.0_dp)

!!IF ONLY ISOTROPIC AVERAGING IS CONSIDERED!!
      ! CALL cvv_only_x(sys%mol_num, sys%framecount, rams%z_para, rams%z_ortho, alpha_diff_x, &
      !               alpha_diff_y, alpha_diff_z, rams%direction, md)

      !CALL dfftw_plan_dft_r2c_1d(plan, 2*md%t_cor, rams%z_para, zhat_para, FFTW_ESTIMATE)
      !CALL dfftw_execute_dft_r2c(plan, rams%z_para, zhat_para)

      ! CALL dfftw_plan_dft_r2c_1d(plan, 2*md%t_cor, rams%z_ortho, zhat_ortho, FFTW_ESTIMATE)
      ! CALL dfftw_execute_dft_r2c(plan, rams%z_ortho, zhat_ortho)

!!ORTHOGONAL!!
      ! OPEN (UNIT=68, FILE='result_fft_water_lib_ortho_iso.txt', STATUS='unknown', IOSTAT=stat)
      ! zhat_ortho = REAL(zhat_ortho, kind=dp)
      ! freq_res = REAL(md%freq_range/md%t_cor, kind=dp)
      ! f = freq_res*md%dt*1.883652d-4

      !DO i = 0, 2*md%t_cor - 2
      !  zhat_ortho(i + 1) = REAL(zhat_ortho(i + 1), kind=dp)*(f*(i + 1)/SIN(f*(i + 1)))**2._dp

      !   zhat_ortho(i) = REAL(zhat_ortho(i), kind=dp)*((const_planck)/(8.0_dp*const_boltz*const_permit*const_permit) &
      !                                                 *1d-29*0.421_dp*md%dt &
      !                                                 *(((rams%laser_in - ((i)*freq_res))**4)/((i)*freq_res)) &
      !                                                 *(1.0_dp/(1.0_dp - EXP((-1.438777_dp*((i)*freq_res)) &
      !                                                                         /gs%temp))))*2.0_dp*pi*2.0_dp

      !     zhat_ortho(0) = 0.0_dp
      !     IF ((i*freq_res).GE.5000_dp) CYCLE
      !     WRITE (68, *) i*freq_res, (REAL(zhat_ortho(i), kind=dp))
      ! END DO
      ! CLOSE (68)

!!PARALLEL!!
      !OPEN (UNIT=67, FILE='result_fft_water_lib_para_iso.txt', STATUS='unknown', IOSTAT=stat)
      !zhat_para = REAL(zhat_para, kind=dp)
      !freq_res = REAL(md%freq_range/md%t_cor, kind=dp)
      !f = freq_res*md%dt*1.883652d-4

      !DO i = 0, 2*md%t_cor - 2
      !zhat_para(i + 1) = REAL(zhat_para(i + 1), kind=dp)*(f*(i + 1)/SIN(f*(i + 1)))**2._dp

      ! zhat_para(i) = REAL(zhat_para(i), kind=dp)*((const_planck)/(8.0_dp*const_boltz*const_permit*const_permit) &
      !                                            *1d-29*0.421_dp*md%dt &
      !                                            *(((rams%laser_in - ((i)*freq_res))**4)/((i)*freq_res)) &
      !                                            *(1.0_dp/(1.0_dp - EXP((-1.438777_dp*((i)*freq_res)) &
      !                                                                    /gs%temp))))*2.0_dp*pi*2.0_dp
      !  zhat_para(0) = 0.0_dp
      !   IF ((i*freq_res).GE.5000_dp) CYCLE
      !    WRITE (67, *) (i)*freq_res, (REAL(zhat_para(i), kind=dp))!,REAL(integral(i),kind=dp)
      ! END DO
      ! CLOSE (67)

!!UNPOL!!
      ! OPEN (UNIT=69, FILE='result_fft_water_lib_unpol_iso.txt', STATUS='unknown', IOSTAT=stat)
      ! freq_res = REAL(md%freq_range/md%t_cor, kind=dp)

      !  DO i = 0, 2*md%t_cor - 2
      !    zhat_unpol_x(i) = zhat_para(i) + zhat_ortho(i)
      !     IF ((i*freq_res).GE.5000_dp) CYCLE
      !      WRITE (69, *) i*freq_res, REAL(zhat_unpol_x(i), kind=dp)
      !   END DO
      !   CLOSE (69)

!!DEPOL RATIO!!
      ! OPEN (UNIT=70, FILE='result_fft_water_lib_depol_iso.txt', STATUS='unknown', IOSTAT=stat)

      ! DO i = 0, 2*md%t_cor - 2
      !    zhat_depol_x(i) = REAL(zhat_ortho(i), kind=dp)/REAL(zhat_para(i), kind=dp)
      !    IF ((i*freq_res).GE.5000_dp) CYCLE
      !     WRITE (70, *) i*freq_res, REAL(zhat_depol_x(i), kind=dp)!REAL(zhat_ortho(i),kind=dp)/REAL(zhat_para(i),kind=dp)
      !  END DO
      !   CLOSE (70)

      !    DEALLOCATE (rams%z_para, rams%z_ortho, zhat_unpol_x, zhat_depol_x, zhat_para, zhat_ortho)
      ! END IF

      ! CALL dfftw_destroy_plan(plan)

      !DEALLOCATE(dip_free,dip_x,dip_y,dip_z)
      !  IF (sys%system=='1') THEN
      !      DEALLOCATE (fragment_x, fragment_y, fragment_z, fragment_free)
      !  DEALLOCATE(natom_frag_x,natom_frag_y,natom_frag_z,natom_frag_free)
      !   END IF

!    DEALLOCATE (alpha_x, alpha_y, alpha_z)
!    DEALLOCATE (alpha_diff_x, alpha_diff_y, alpha_diff_z)

   END SUBROUTINE spec_raman

!!....................................................................................................................!
!....................................................................................................................!
   SUBROUTINE normal_mode_analysis(sys, stats)
      TYPE(systems), INTENT(INOUT)        :: sys
      TYPE(static), INTENT(INOUT)        :: stats

      INTEGER                                                     :: stat, i, j, m, n, p, k, info, lwork, lwmax, lda
      REAL(kind=dp)                                                :: factor
      REAL(kind=dp), DIMENSION(:), ALLOCATABLE                       :: w, work, w_new
      REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE                     :: hessian_new, atomic_displacements
      REAL(kind=dp), DIMENSION(:, :, :, :), ALLOCATABLE                 :: hessian
      LOGICAL, DIMENSION(9)                                        :: mk = .TRUE.

      lwmax = MAX(1, 3*sys%natom*64)  ! conservative guess; LAPACK recommends this for DSYEV
      lda = sys%natom*3
      stats%nmodes = 3*sys%natom - 6 !only for non-linear molecules

      ALLOCATE (work(lwmax), w(sys%natom*3), w_new(sys%natom*3))

      factor = REAL(1.0_dp/(2.0_dp*stats%dx), kind=dp)

      ALLOCATE (hessian(0:sys%natom - 1, 0:2, 0:sys%natom - 1, 0:2), hessian_new(0:sys%natom*3 - 1, 0:sys%natom*3 - 1))

!hessian=factor*hessian_factor*(stats%force(2,:,:,:,:)-stats%force(1,:,:,:,:))
      hessian = hartreebohr2evang*factor*hessian_factor*(stats%force(2, :, :, :, :) - stats%force(1, :, :, :, :))

      p = 0
      DO i = 0, sys%natom - 1
         DO m = 0, 2
            k = 0
            DO j = 0, sys%natom - 1
               DO n = 0, 2
                  hessian_new(i + m + p, j + n + k) = sys%mass_mat(i + 1, j + 1)*hessian(i, m, j, n)
               END DO
               k = k + 2
            END DO
         END DO
         p = p + 2
      END DO

!hessian_new(:,:)=RESHAPE(hessian(:,:,:,:), (/3*sys%natom, 3*sys%natom/))
      hessian_new(:, :) = REAL((hessian_new(:, :) + TRANSPOSE(hessian_new(:, :)))/2.0_dp, kind=dp)
      n = SIZE(hessian_new, 1)

      PRINT *, hessian_new(1, 1), "hess"

! work size query
      lwork = -1
      CALL DSYEV('V', 'U', n, hessian_new, lda, w, work, lwork, info)
      lwork = MIN(lwmax, INT(work(1)))

! get eigenvalues and eigenvectors
      CALL dsyev('V', 'U', n, hessian_new, lda, w, work, lwork, info)

      hessian_new = TRANSPOSE(hessian_new)

      w = REAL(w*SQRT(ABS(w))/ABS(w), kind=dp)
      w = REAL(w/(2.0_dp*pi*speed_light), kind=dp)

      ALLOCATE (stats%freq(stats%nmodes), atomic_displacements(stats%nmodes, sys%natom*3), stats%disp(stats%nmodes, sys%natom, 3))

      DO i = 7, sys%natom*3
         stats%freq(i - 6) = w(i)
      END DO

      atomic_displacements(1:stats%nmodes, 1:sys%natom*3) = hessian_new(6:3*sys%natom - 1, :)

      m = 0
      DO j = 0, sys%natom - 1 !sys%natom
         DO k = 0, 2 !dims
            stats%disp(1:stats%nmodes, j + 1, k + 1) = atomic_displacements(1:stats%nmodes, j + k + 1 + m)
         END DO
         m = m + 2
      END DO

      PRINT *, stats%freq(1:3)

      OPEN (UNIT=13, FILE='normal_mode_freq.txt', STATUS='unknown', IOSTAT=stat)
      DO i = 1, stats%nmodes !!atom_num: 1st atom
         WRITE (13, *) stats%freq(i)
      END DO

      OPEN (UNIT=14, FILE='normal_mode_displ.txt', STATUS='unknown', IOSTAT=stat)
      DO i = 1, stats%nmodes !!atom_num: 1st atom
         DO j = 1, sys%natom !!dims: x dimension
            WRITE (14, *) stats%disp(i, j, 1:3)
         END DO
      END DO

   END SUBROUTINE normal_mode_analysis

!....................................................................................................................!
!....................................................................................................................!

   SUBROUTINE spec_static_ir(gs, sys, stats, dips)

      TYPE(global_settings), INTENT(INOUT)        :: gs
      TYPE(systems), INTENT(INOUT)        :: sys
      TYPE(static), INTENT(INOUT)        :: stats
      TYPE(dipoles), INTENT(INOUT)        :: dips

      REAL(kind=dp), DIMENSION(:), ALLOCATABLE        :: ir_int
      INTEGER                                                  :: stat, i, k, x, freq_res
      INTEGER                                                  :: start_freq, end_freq
      REAL(kind=dp), DIMENSION(:), ALLOCATABLE                    :: gamma_sq, data2!,broad
      REAL(kind=dp)                                             :: broad, ir_factor

      ALLOCATE (gamma_sq(stats%nmodes), ir_int(stats%nmodes))

      start_freq = 1
      end_freq = INT(MAXVAL(stats%freq) + 1000.0_dp)
      freq_res = INT(end_freq - start_freq)

      PRINT *, "Max freq: ", MAXVAL(stats%freq)
      PRINT *, "end_freq: ", end_freq
      PRINT *, "freq_res: ", freq_res

      ALLOCATE (data2(freq_res + 1))
      data2 = 0.0_dp

      DO k = 1, stats%nmodes
         gamma_sq(k) = SQRT(DOT_PRODUCT(dips%dip_dq(k, :), dips%dip_dq(k, :)))
      END DO

      !IR factor used in cp2k which is equal to 42.256
      !! first convert debye²angstrom⁻²amu⁻¹ to C^2/kg then to km/mol
      ir_factor = (debye2cm/ang)**2.0_dp*avo_num*1.0e-3_dp/(12.0_dp*const_permit*(speed_light*cm2m)**2.0_dp)/am_u
      ir_int(:) = (gamma_sq(:)**2.0_dp)*ir_factor

    !!!Broadening the spectrum!!
      DO i = start_freq, end_freq
         broad = 0.0_dp
         DO x = 1, stats%nmodes
            broad = broad + (ir_int(x)*(1.0_dp/(gs%fwhm*SQRT(2.0_dp*pi)))*EXP(-0.50_dp*((i - stats%freq(x))/gs%fwhm)**2.0_dp))
         END DO
         data2(i) = data2(i) + broad
      END DO

      OPEN (UNIT=98, FILE='result_static_ir.txt', STATUS='unknown', IOSTAT=stat)
      DO i = start_freq, end_freq
         WRITE (98, *) i, data2(i)
      END DO
      CLOSE (98)

      DEALLOCATE (gamma_sq, data2, ir_int, dips%dip_dq, stats%freq, stats%disp)

   END SUBROUTINE spec_static_ir
!....................................................................................................................!
!....................................................................................................................!
   SUBROUTINE spec_static_raman(gs, sys, stats, dips, rams)

      TYPE(global_settings), INTENT(INOUT)        :: gs
      TYPE(systems), INTENT(INOUT)        :: sys
      TYPE(static), INTENT(INOUT)        :: stats
      TYPE(dipoles), INTENT(INOUT)        :: dips
      TYPE(raman), INTENT(INOUT)        :: rams

      INTEGER                                                  :: stat, i, j, x, freq_res
      INTEGER                                                  :: start_freq, end_freq
      REAL(kind=dp), DIMENSION(:), ALLOCATABLE                    :: iso_sq, aniso_sq, ram_const, data2!,broad
      REAL(kind=dp)                                             :: broad

      ALLOCATE (iso_sq(stats%nmodes), aniso_sq(stats%nmodes))
      ALLOCATE (rams%raman_int(stats%nmodes), ram_const(stats%nmodes))

      start_freq = 1
      end_freq = INT(MAXVAL(stats%freq) + 1000.0_dp)
      freq_res = INT(end_freq - start_freq)
      ALLOCATE (data2(freq_res + 1))
      data2 = 0.0_dp

    !!!Isotropic and anisotropic contributions!!
      iso_sq(:) = REAL((rams%pol_dq(:, 1, 1) + rams%pol_dq(:, 2, 2) + rams%pol_dq(:, 3, 3))/3.0_dp, kind=dp)**2.0_dp

      aniso_sq(:) = (0.50_dp*(((rams%pol_dq(:, 1, 1) - rams%pol_dq(:, 2, 2))**2.0_dp) + &
                              ((rams%pol_dq(:, 2, 2) - rams%pol_dq(:, 3, 3))**2.0_dp) &
                              + ((rams%pol_dq(:, 3, 3) - rams%pol_dq(:, 1, 1))**2.0_dp))) &
                    + (3.0_dp*((rams%pol_dq(:, 1, 2)**2.0_dp) + (rams%pol_dq(:, 2, 3)**2.0_dp) &
                               + (rams%pol_dq(:, 3, 1)**2.0_dp)))

    !!!Conversion from angstrom^4 amu⁻¹ to m^4 kg -1
      iso_sq = iso_sq*(ang**4._dp)/am_u
      aniso_sq = aniso_sq*(ang**4._dp)/am_u

    !!! Conversion of static Raman units into 10^{-30}*cm^2/sr
      ram_const(:) = (const_planck/(8.0_dp*speed_light*cm2m*const_permit*const_permit)*1.e+30* &
                      REAL(((rams%laser_in/reccm2ev - stats%freq(:))**4.0_dp)/(stats%freq(:)*cm2m**3.0_dp), kind=dp)* &
                      (1.0_dp/(1.0_dp - EXP(-1._dp*const_planck*speed_light*cm2m*stats%freq(:)/ &
                                            (const_boltz*gs%temp)))))/(cm2m**2._dp)

    !!!Calculation of the unpolarized Raman intensities!!
      rams%raman_int(:) = REAL(((7.0_dp*aniso_sq(:)) + (45.0_dp*iso_sq(:)))/45.0_dp, kind=dp)*ram_const(:)

    !!!Broadening the spectrum!!
      DO i = start_freq, end_freq
         broad = 0.0_dp
         DO x = 1, stats%nmodes
            broad = broad + (rams%raman_int(x)*(1.0_dp/(gs%fwhm*SQRT(2.0_dp*pi)))* &
                             EXP(-0.50_dp*((i - stats%freq(x))/gs%fwhm)**2.0_dp))
         END DO
         data2(i) = data2(i) + broad
      END DO

      OPEN (UNIT=98, FILE='result_static_raman.txt', STATUS='unknown', IOSTAT=stat)
      DO i = start_freq, end_freq
         WRITE (98, *) i, data2(i)
      END DO
      CLOSE (98)

    !!Write Molden output
      rams%raman_int = REAL(rams%raman_int/MINVAL(rams%raman_int), kind=dp)
      OPEN (UNIT=15, FILE='raman.mol', STATUS='unknown', IOSTAT=stat)
      WRITE (15, *) "[Molden Format]"
      WRITE (15, *) "[GEOMETRIES] XYZ"
      WRITE (15, *) sys%natom
      WRITE (15, *)
      DO i = 1, sys%natom
         WRITE (15, *) sys%element(i), sys%coord(i, 1), sys%coord(i, 2), sys%coord(i, 3)
      END DO
      WRITE (15, *) "[stats%freq]"
      DO i = 1, stats%nmodes
         WRITE (15, *) stats%freq(i)
      END DO
      WRITE (15, *) "[INT]"
      DO i = 1, stats%nmodes
         WRITE (15, *) rams%raman_int(i)
      END DO
      WRITE (15, *) "[FR-sys%coord]"
      WRITE (15, *) sys%natom
      WRITE (15, *)
      DO i = 1, sys%natom
         WRITE (15, *) sys%element(i), sys%coord(i, 1)/bohr2ang, sys%coord(i, 2)/bohr2ang, sys%coord(i, 3)/bohr2ang
      END DO
      WRITE (15, *) "[FR-NORM-sys%coord]"
      DO i = 1, stats%nmodes
         WRITE (15, *) "vibration", i
         DO j = 1, sys%natom
            WRITE (15, *) stats%disp(i, j, 1)/bohr2ang, stats%disp(i, j, 2)/bohr2ang, stats%disp(i, j, 3)/bohr2ang
         END DO
      END DO
      CLOSE (15)

      DEALLOCATE (iso_sq, aniso_sq, data2, ram_const)

   END SUBROUTINE spec_static_raman
!....................................................................................................................!
!....................................................................................................................!

   SUBROUTINE spec_abs(gs, sys, dips, rams)
      TYPE(global_settings), INTENT(INOUT)        :: gs
      TYPE(systems), INTENT(INOUT)        :: sys
      TYPE(dipoles), INTENT(INOUT)        :: dips
      TYPE(raman), INTENT(INOUT)        :: rams
      CHARACTER(len=256) :: filename
      INTEGER                                                       :: stat, i, j, k, m, x, o, dims, dir
      INTEGER(kind=dp)                                               :: plan
      REAL(kind=dp)                                                  :: rtp_freq_res, freq_au
      REAL(kind=dp), DIMENSION(:, :, :, :), ALLOCATABLE                   :: trace, abs_intens
      COMPLEX(kind=dp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE            :: y_out

      IF (gs%spectral_type%read_function == 'RR') THEN
         dims = 3
         dir = 2
      ELSEIF (gs%spectral_type%read_function == 'ABS') THEN
         dims = 1
         dir = 1
         sys%natom = 1
      END IF

      ALLOCATE (rams%RR%zhat_pol_rtp(sys%natom, dims, dir, 3, 3, rams%RR%framecount_rtp))

      rams%RR%zhat_pol_rtp = COMPLEX(0._dp, 0.0_dp)

        !!!FFT of the RTP polarizabilities
      DO j = 1, sys%natom
         DO i = 1, dims
            DO k = 1, dir
               DO m = 1, 3
                  DO o = 1, 3
                     CALL dfftw_plan_dft_r2c_1d(plan, rams%RR%framecount_rtp, &
                                                rams%RR%pol_rtp(j, i, k, m, o, 1:rams%RR%framecount_rtp), &
                                                rams%RR%zhat_pol_rtp(j, i, k, m, o, 1:rams%RR%framecount_rtp), FFTW_ESTIMATE)
                     CALL dfftw_execute_dft_r2c(plan, rams%RR%pol_rtp(j, i, k, m, o, 1:rams%RR%framecount_rtp), &
                                                rams%RR%zhat_pol_rtp(j, i, k, m, o, 1:rams%RR%framecount_rtp))
                     CALL dfftw_destroy_plan(plan)
                  END DO
               END DO
            END DO
         END DO
      END DO

      IF (rams%RR%check_pade == 'y') THEN

         ALLOCATE (y_out(sys%natom, dims, dir, 3, 3, rams%RR%framecount_rtp_pade))
        !!Call Pade
!$OMP PARALLEL DO COLLAPSE(5)
         DO j = 1, sys%natom
            DO i = 1, dims
               DO k = 1, dir
                  DO m = 1, 3
                     DO o = 1, 3
                        CALL interpolate(rams%RR%framecount_rtp, rams%RR%zhat_pol_rtp(j, i, k, m, o, 1:rams%RR%framecount_rtp), &
                                         rams%RR%framecount_rtp_pade, y_out(j, i, k, m, o, :))
                     END DO
                  END DO
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO
         rams%RR%framecount_rtp = rams%RR%framecount_rtp_pade
         rams%RR%zhat_pol_rtp = y_out
         DEALLOCATE (y_out)
      END IF

!!!Dividing by electric field and multiplying by rams%RR%dt_rtp which is coming from FFT
      !rams%RR%zhat_pol_rtp = rams%RR%zhat_pol_rtp/dips%e_field
      rams%RR%zhat_pol_rtp = rams%RR%zhat_pol_rtp*(rams%RR%dt_rtp*fs2s)/dips%e_field

!!!Finding frequency range
      rtp_freq_res = REAL(rams%RR%freq_range_rtp/rams%RR%framecount_rtp, kind=dp)

!!!Calculate absorption spectra
      ALLOCATE (trace(sys%natom, dims, dir, rams%RR%framecount_rtp))
      ALLOCATE (abs_intens(sys%natom, dims, dir, rams%RR%framecount_rtp))

      trace = 0.0_dp
      trace(:, :, :, :) = DIMAG(rams%RR%zhat_pol_rtp(:, :, :, 1, 1, :)) + DIMAG(rams%RR%zhat_pol_rtp(:, :, :, 2, 2, :)) &
                          + DIMAG(rams%RR%zhat_pol_rtp(:, :, :, 3, 3, :))

      !!Conversion of absorption spectrum units into a.u.
      abs_intens(:, :, :, :) = (4.0_dp*pi*debye*trace(:, :, :, :))/(3.0_dp*speed_light_au*at_u)

      ! DO j = 1, sys%natom      !! shifted atom index
      ! DO i = 1, dims          !! displacement direction
      ! DO k = 1, dir      !! + / - shift
      !    Create a unique file name for this shifted structure
      !  WRITE (filename, '("absorption_spectrum_",I0,"_",I0,"_",I0,".txt")') j, i, k
      !  OPEN (UNIT=13, FILE=filename, STATUS='unknown', IOSTAT=stat)
      !  IF (stat /= 0) THEN
      !     PRINT *, "Error opening file: ", TRIM(filename)
      !      STOP
      !    END IF

      !   Loop over time steps to write the full spectrum
      !     DO o = 1, rams%RR%framecount_rtp
      !         WRITE (13, *) o*rtp_freq_res*reccm2ev, &
      !             abs_intens(j, i, k, o)*o*rtp_freq_res*(-1.0_dp)
      !        END DO

      !         CLOSE (13)  ! close after finishing one structure
      !       END DO
      !        END DO
!      END DO

      !! Conversion from cm-1 to a.u.
      freq_au = rtp_freq_res*(-1.0_dp)*reccm2au

      DO j = 1, 1 !!atom_num: 1st atom
         DO i = 1, 1 !!dims: x dimension
            DO k = 1, 1 !! + direction
               OPEN (UNIT=13, FILE='absorption_spectrum.txt', STATUS='unknown', IOSTAT=stat)
               DO o = 1, rams%RR%framecount_rtp
                  WRITE (13, *) o*rtp_freq_res*reccm2ev, abs_intens(j, i, k, o)*o*freq_au
               END DO
            END DO
         END DO
      END DO
      CLOSE (13)
      DEALLOCATE (rams%RR%pol_rtp, trace, abs_intens)

   END SUBROUTINE spec_abs
!....................................................................................................................!
!....................................................................................................................!
   SUBROUTINE spec_static_resraman(gs, sys, stats, rams)

      TYPE(global_settings), INTENT(INOUT)        :: gs
      TYPE(systems), INTENT(INOUT)        :: sys
      TYPE(static), INTENT(INOUT)        :: stats
      TYPE(raman), INTENT(INOUT)        :: rams

      INTEGER                                                       :: stat, i, j, k, m, x, freq_res, l, o, n, r
      INTEGER                                                       :: start_freq, end_freq, rtp_point
      INTEGER(kind=dp)                                               :: plan
      REAL(kind=dp)                                                  :: broad, factor
      REAL(kind=dp)                                                  :: rtp_freq_res, pade_freq_res
      REAL(kind=dp), DIMENSION(:), ALLOCATABLE                         :: data2, ram_const
      REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE                       :: iso_sq, aniso_sq, raman_int
      REAL(kind=dp), DIMENSION(:, :, :, :), ALLOCATABLE                   :: zhat_pol_dq_rtp
      REAL(kind=dp), DIMENSION(:, :, :, :, :), ALLOCATABLE                 :: zhat_pol_dxyz_rtp

      broad = 0.0_dp
      factor = 1._dp/(2.0_dp*stats%dx)
      start_freq = 1.0_dp
      end_freq = INT(MAXVAL(stats%freq) + 1000.0_dp)
      freq_res = INT(end_freq - start_freq)

      ALLOCATE (data2(freq_res*stats%nmodes))
      ALLOCATE (zhat_pol_dxyz_rtp(sys%natom, 3, 3, 3, rams%RR%framecount_rtp))
      ALLOCATE (zhat_pol_dq_rtp(stats%nmodes, 3, 3, rams%RR%framecount_rtp))
      ALLOCATE (iso_sq(stats%nmodes, rams%RR%framecount_rtp), aniso_sq(stats%nmodes, rams%RR%framecount_rtp))
      ALLOCATE (raman_int(stats%nmodes, rams%RR%framecount_rtp), ram_const(stats%nmodes))

      zhat_pol_dq_rtp = 0.0_dp
      data2 = 0.0_dp

!!!Finding laser frequency
      rtp_freq_res = REAL(rams%RR%freq_range_rtp/rams%RR%framecount_rtp, kind=dp)
      rtp_point = ANINT(rams%laser_in/(rtp_freq_res*reccm2ev), kind=dp)

      PRINT *, rams%laser_in, "rams%laser_in", rtp_freq_res, "rtp_freq_res", &
         rams%RR%freq_range_rtp, "rams%RR%freq_range_rtp", rtp_point, 'rtp_point', rams%RR%framecount_rtp

!!!Finite differences
      zhat_pol_dxyz_rtp(:, :, :, :, :) = REAL((REAL(rams%RR%zhat_pol_rtp(:, :, 2, :, :, :), kind=dp) &
                                               - REAL(rams%RR%zhat_pol_rtp(:, :, 1, :, :, :), kind=dp))*factor, kind=dp)

!!!Derivatives w.r.t. mass weighted normal coordinates
      DO i = 1, stats%nmodes
         DO j = 1, sys%natom
            DO o = 1, rams%RR%framecount_rtp
               zhat_pol_dq_rtp(i, :, :, o) = zhat_pol_dq_rtp(i, :, :, o) &
                                             + (zhat_pol_dxyz_rtp(j, 1, :, :, o)*stats%disp(i, j, 1)*sys%atom_mass_inv_sqrt(j)) &
                                             + (zhat_pol_dxyz_rtp(j, 2, :, :, o)*stats%disp(i, j, 2)*sys%atom_mass_inv_sqrt(j)) &
                                             + (zhat_pol_dxyz_rtp(j, 3, :, :, o)*stats%disp(i, j, 3)*sys%atom_mass_inv_sqrt(j))
            END DO
         END DO
      END DO

!!!Isotropic and anisotropic contributions!!
      iso_sq(:, :) = REAL((zhat_pol_dq_rtp(:, 1, 1, :) + zhat_pol_dq_rtp(:, 2, 2, :) &
                           + zhat_pol_dq_rtp(:, 3, 3, :))/3.0_dp, kind=dp)**2.0_dp

      aniso_sq(:, :) = (0.50_dp*(((zhat_pol_dq_rtp(:, 1, 1, :) - zhat_pol_dq_rtp(:, 2, 2, :))**2.0_dp) + &
                                 ((zhat_pol_dq_rtp(:, 2, 2, :) - zhat_pol_dq_rtp(:, 3, 3, :))**2.0_dp) &
                                 + ((zhat_pol_dq_rtp(:, 3, 3, :) - zhat_pol_dq_rtp(:, 1, 1, :))**2.0_dp))) &
                       + (3.0_dp*((zhat_pol_dq_rtp(:, 1, 2, :)**2.0_dp) &
                                  + (zhat_pol_dq_rtp(:, 2, 3, :)**2.0_dp) &
                                  + (zhat_pol_dq_rtp(:, 3, 1, :)**2.0_dp)))

!!!Conversion from (debye/E)^2 angstrom^-2 amu⁻¹ to angstrom^6 angstrom^-2 amu⁻¹
      iso_sq = iso_sq/(a3_to_debye_per_e*a3_to_debye_per_e)
      aniso_sq = aniso_sq/(a3_to_debye_per_e*a3_to_debye_per_e)

!!!Conversion from angstrom^4 amu⁻¹ to m^4 kg^-1
      iso_sq = iso_sq*(ang**4._dp)/am_u
      aniso_sq = aniso_sq*(ang**4._dp)/am_u

!!! Conversion of static resonance Raman units into 10^{-30}*cm^2/sr
      ram_const(:) = (const_planck/(8.0_dp*speed_light*cm2m*const_permit*const_permit)*1.e+30* &
                      REAL(((rams%laser_in/reccm2ev - stats%freq(:))**4.0_dp)/(stats%freq(:)*cm2m**3.0_dp), kind=dp)* &
                      (1.0_dp/(1.0_dp - EXP(-1._dp*const_planck*speed_light*cm2m*stats%freq(:)/ &
                                            (const_boltz*gs%temp)))))/(cm2m**2._dp)

!!!Calculation of the unpolarized resonance Raman intensities!!
      raman_int(:, rtp_point) = REAL(((7.0_dp*aniso_sq(:, rtp_point)) + (45.0_dp*iso_sq(:, rtp_point)))/45.0_dp, kind=dp)* &
                                ram_const(:)

!!!Broadening the spectrum!!
      DO x = start_freq, end_freq
         broad = 0.0_dp
         DO i = 1, stats%nmodes
            broad = broad + (raman_int(i, rtp_point)*(1.0_dp/(gs%fwhm*SQRT(2.0_dp*pi))) &
                             *EXP(-0.50_dp*((x - stats%freq(i))/gs%fwhm)**2.0_dp))
         END DO
         data2(x) = data2(x) + broad
      END DO

      OPEN (UNIT=98, FILE='result_static_resraman.txt', STATUS='unknown', IOSTAT=stat)
      DO i = start_freq, end_freq
         WRITE (98, *) i, data2(i)
      END DO
      CLOSE (98)

      DEALLOCATE (iso_sq, aniso_sq, data2, ram_const, raman_int, stats%disp, stats%freq)
      DEALLOCATE (rams%RR%zhat_pol_rtp, zhat_pol_dxyz_rtp, zhat_pol_dq_rtp)

   END SUBROUTINE spec_static_resraman

!!....................................................................................................................!
!!....................................................................................................................!
!
   SUBROUTINE spec_resraman(natom, framecount, element, rtp_dipole_x, rtp_dipole_y, rtp_dipole_z, type_input, &
                            mol_num, system, read_function, dt, z_iso_resraman, z_aniso_resraman, freq_range, &
                            freq_range_rtp, laser_in_resraman, y_out)

      CHARACTER(LEN=40), INTENT(INOUT)                          :: read_function, system
      CHARACTER(LEN=40), INTENT(INOUT)                          :: type_input
      CHARACTER(LEN=40), INTENT(INOUT)                          :: rtp_dipole_x, rtp_dipole_y, rtp_dipole_z
      CHARACTER(LEN=2), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: element
      INTEGER, INTENT(INOUT)                                    :: natom, framecount, mol_num
      REAL(kind=dp), INTENT(INOUT)                               :: dt, freq_range, freq_range_rtp
      REAL(kind=dp), INTENT(INOUT)                               :: laser_in_resraman
      COMPLEX(kind=dp), DIMENSION(:, :), ALLOCATABLE, INTENT(INOUT) :: z_iso_resraman, z_aniso_resraman
      COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)   :: y_out

      TYPE(global_settings) :: gs
      TYPE(systems)         :: sys
      TYPE(molecular_dynamics)   :: md
      TYPE(dipoles)   :: dips
      CHARACTER(LEN=40)                                        :: chara
      INTEGER                                                  :: stat, i, j, k, m, t0, t1
      INTEGER(kind=dp)                                          :: plan
      COMPLEX(kind=dp), DIMENSION(:, :, :), ALLOCATABLE             :: yx_out, yy_out, yz_out
      COMPLEX(kind=dp), DIMENSION(:, :), ALLOCATABLE               :: zhat_iso_resraman, zhat_aniso_resraman
      COMPLEX(kind=dp), DIMENSION(:, :, :), ALLOCATABLE             :: zhat_resraman_x, zhat_resraman_y, zhat_resraman_z
      REAL(kind=dp)                                             :: f, freq_res, rtp_freq_res, pade_freq_res, laser_in
      REAL(kind=dp), DIMENSION(:), ALLOCATABLE                    :: trace, abs_intens, trace_pade, abs_intens_pade
      REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE                  :: zhat_unpol_resraman
      REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE                :: alpha_resraman_x, alpha_resraman_y, alpha_resraman_z
      REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE                :: alpha_resraman_x_diff_re, alpha_resraman_y_diff_re
      REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE                :: alpha_resraman_z_diff_re
      REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE                :: alpha_resraman_x_diff_im, alpha_resraman_y_diff_im
      REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE                :: alpha_resraman_z_diff_im
      REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE                :: alpha_resraman_x_im, alpha_resraman_y_im
      REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE                :: alpha_resraman_z_im
      REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE                :: alpha_resraman_x_re, alpha_resraman_y_re
      REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE                :: alpha_resraman_z_re
      REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE                :: alpha_x, alpha_y, alpha_z
      REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE                :: dip_x, dip_y, dip_z

      ALLOCATE (zhat_resraman_x(framecount, natom, 3), zhat_resraman_y(framecount, natom, 3), &
                zhat_resraman_z(framecount, natom, 3))
      ALLOCATE (alpha_resraman_x(framecount, natom, 3), alpha_resraman_y(framecount, natom, 3), &
                alpha_resraman_z(framecount, natom, 3))
      ALLOCATE (alpha_resraman_x_diff_re(framecount - 2, natom, 3), alpha_resraman_y_diff_re(framecount - 2, natom, 3), &
                alpha_resraman_z_diff_re(framecount - 2, natom, 3))
      ALLOCATE (alpha_resraman_x_diff_im(framecount - 2, natom, 3), alpha_resraman_y_diff_im(framecount - 2, natom, 3), &
                alpha_resraman_z_diff_im(framecount - 2, natom, 3))
      ALLOCATE (alpha_x(framecount, natom, 3), alpha_y(framecount, natom, 3), alpha_z(framecount, natom, 3))
      ALLOCATE (yx_out(framecount, 10000, 3), yy_out(framecount, 10000, 3), yz_out(framecount, 10000, 3))

!!X-Field!!
      CALL read_coord_frame(mol_num, rtp_dipole_x, dip_x, sys)
      CALL forward_diff(mol_num, alpha_x, dip_x, dip_x, gs, sys, dips)

      zhat_resraman_x = COMPLEX(0._dp, 0.0_dp)
      zhat_resraman_y = COMPLEX(0._dp, 0.0_dp)
      zhat_resraman_z = COMPLEX(0._dp, 0.0_dp)

      DO i = 1, framecount
         DO j = 1, 3
            CALL dfftw_plan_dft_r2c_1d(plan, natom, alpha_x(i, 1:natom, j), zhat_resraman_x(i, 1:natom, j), FFTW_ESTIMATE)
            CALL dfftw_execute_dft_r2c(plan, alpha_x(i, 1:natom, j), zhat_resraman_x(i, 1:natom, j)) !!!important to specify arrays!!
            CALL dfftw_destroy_plan(plan)
         END DO
      END DO

      alpha_resraman_x_re = REAL(zhat_resraman_x, kind=dp)
      alpha_resraman_x_im = AIMAG(zhat_resraman_x)

!!Call Pade
      DO i = 1, 1 !!framecount
         DO j = 1, 3
            CALL interpolate(natom - 1, zhat_resraman_x(i, 1:natom, j), 10000, yx_out(i, :, j))
         END DO
      END DO

!OPEN(UNIT=40,FILE='y_out.txt',STATUS='unknown',IOSTAT=stat)
!DO i=1,10000
      ! WRITE(40,*) i/10000._dp,REAL(yx_out(i),kind=dp),AIMAG(yx_out(i))
!ENDDO
!CLOSE(40)

!OPEN(UNIT=40,FILE='zhat_out.txt',STATUS='unknown',IOSTAT=stat)
!DO i=1,natom-1
      ! WRITE(40,*) i/(natom-1._dp),REAL(zhat_resraman_x(1,i,1),kind=dp),AIMAG(zhat_resraman_x(1,i,1))
!ENDDO
!CLOSE(40)

      CALL central_diff(natom, alpha_resraman_x_re, alpha_resraman_x_diff_re, sys, md)
      CALL central_diff(natom, alpha_resraman_x_im, alpha_resraman_x_diff_im, sys, md)

!!Y-Field!!

      CALL read_coord_frame(natom, rtp_dipole_y, dip_y, sys)
      CALL forward_diff(mol_num, alpha_y, dip_y, dip_y, gs, sys, dips)

      DO i = 1, framecount
         DO j = 1, 3
            CALL dfftw_plan_dft_r2c_1d(plan, natom, alpha_y(i, 1:natom, j), zhat_resraman_y(i, 1:natom, j), FFTW_ESTIMATE)
            CALL dfftw_execute_dft_r2c(plan, alpha_y(i, 1:natom, j), zhat_resraman_y(i, 1:natom, j)) !!!important to specify arrays!!
            CALL dfftw_destroy_plan(plan)
         END DO
      END DO

      alpha_resraman_y_re = REAL(zhat_resraman_x, kind=dp)
      alpha_resraman_y_im = AIMAG(zhat_resraman_x)

!!Call Pade
      DO i = 1, 1 !!framecount
         DO j = 1, 3
            CALL interpolate(natom, zhat_resraman_y(i, 1:natom, j), 10000, yy_out(i, :, j))
         END DO
      END DO

!OPEN(UNIT=40,FILE='yy_out.txt',STATUS='unknown',IOSTAT=stat)
!DO i=1,10000
!  WRITE(40,*) i/10000._dp,REAL(yy_out(i),kind=dp),AIMAG(yy_out(i))
!ENDDO
!CLOSE(40)

!pade_interpolation.f90OPEN(UNIT=40,FILE='zhaty_out.txt',STATUS='unknown',IOSTAT=stat)
!DO i=1,natom-1
!  WRITE(40,*) i/(natom-1._dp),REAL(zhat_resraman_y(1,i,2),kind=dp),AIMAG(zhat_resraman_y(1,i,2))
!ENDDO
!CLOSE(40)

      CALL central_diff(natom, alpha_resraman_y_re, alpha_resraman_y_diff_re, sys, md)
      CALL central_diff(natom, alpha_resraman_y_im, alpha_resraman_y_diff_im, sys, md)

!!Z-Field!!

      CALL read_coord_frame(natom, rtp_dipole_z, dip_z, sys)
      CALL forward_diff(mol_num, alpha_z, dip_z, dip_z, gs, sys, dips)

      DO i = 1, framecount
         DO j = 1, 3
            CALL dfftw_plan_dft_r2c_1d(plan, natom, alpha_z(i, 1:natom - 1, j), zhat_resraman_z(i, 1:natom, j), FFTW_ESTIMATE)
            CALL dfftw_execute_dft_r2c(plan, alpha_z(i, 1:natom, j), zhat_resraman_z(i, 1:natom, j)) !!!important to specify arrays!!
            CALL dfftw_destroy_plan(plan)
         END DO
      END DO

      alpha_resraman_z_re = REAL(zhat_resraman_x, kind=dp)
      alpha_resraman_z_im = AIMAG(zhat_resraman_x)

!!Call Pade
      DO i = 1, 1 !!framecount
         DO j = 1, 3
            CALL interpolate(natom, zhat_resraman_z(i, 1:natom, j), 10000, yz_out(i, :, j))
         END DO
      END DO

!OPEN(UNIT=40,FILE='yz_out.txt',STATUS='unknown',IOSTAT=stat)
!DO i=1,10000
!  WRITE(40,*) i/10000._dp,REAL(yz_out(i),kind=dp),AIMAG(yz_out(i))
!ENDDO
!CLOSE(40)

!OPEN(UNIT=40,FILE='zhatz_out.txt',STATUS='unknown',IOSTAT=stat)
!DO i=1,natom-1
!  WRITE(40,*) i/(natom-1._dp),REAL(zhat_resraman_z(1,i,3),kind=dp),AIMAG(zhat_resraman_z(1,i,3))
!ENDDO
!CLOSE(40)

      CALL central_diff(natom, alpha_resraman_z_re, alpha_resraman_z_diff_re, sys, md)
      CALL central_diff(natom, alpha_resraman_z_im, alpha_resraman_z_diff_im, sys, md)

!!!Calculate absorption spectra

      rtp_freq_res = REAL(freq_range_rtp/(natom), kind=dp)
      pade_freq_res = REAL(freq_range_rtp/(10000), kind=dp)

      ALLOCATE (abs_intens(natom), trace(natom))
      ALLOCATE (abs_intens_pade(10000), trace_pade(10000))

      DO i = 1, 1
         DO j = 1, natom
            trace(j) = DIMAG(zhat_resraman_x(i, j, 1)) + DIMAG(zhat_resraman_y(i, j, 2)) + DIMAG(zhat_resraman_z(i, j, 3))
            abs_intens(j) = (4.0_dp*pi*trace(j))/(3.0_dp*speed_light)
         END DO
      END DO

      OPEN (UNIT=41, FILE='absorption_spectra.txt', STATUS='unknown', IOSTAT=stat)
      DO i = 1, natom
         WRITE (41, *) i*rtp_freq_res*1.23984198e-4, abs_intens(i)*i*rtp_freq_res
      END DO
      CLOSE (41)

      DO i = 1, 1
         DO j = 1, 10000
            trace_pade(j) = DIMAG(yx_out(i, j, 1)) + DIMAG(yy_out(i, j, 2)) + DIMAG(yz_out(i, j, 3))
            abs_intens_pade(j) = (4.0_dp*pi*trace_pade(j))/(3.0_dp*speed_light)
         END DO
      END DO

      OPEN (UNIT=42, FILE='absorption_spectra_pade.txt', STATUS='unknown', IOSTAT=stat)
      DO i = 1, 10000
         WRITE (42, *) i*pade_freq_res*1.23984198e-4, abs_intens_pade(i)*i*pade_freq_res
      END DO
      CLOSE (42)

!!Generate the spectrum!!

      CALL cvv_resraman(framecount, natom, dt, alpha_resraman_x_diff_re, alpha_resraman_y_diff_re, &
                        alpha_resraman_z_diff_re, alpha_resraman_x_diff_im, alpha_resraman_y_diff_im, alpha_resraman_z_diff_im, &
                        z_iso_resraman, z_aniso_resraman, md)

      ALLOCATE (zhat_iso_resraman(0:md%t_cor*2, natom), zhat_aniso_resraman(0:md%t_cor*2, natom))
      ALLOCATE (zhat_unpol_resraman(0:md%t_cor*2, natom))

      zhat_iso_resraman = COMPLEX(0._dp, 0.0_dp)
      zhat_aniso_resraman = COMPLEX(0._dp, 0.0_dp)

      DO j = 1, natom
         CALL dfftw_plan_dft_1d(plan, 2*md%t_cor, z_iso_resraman(0:md%t_cor*2, j), zhat_iso_resraman(0:md%t_cor*2, j), &
                                FFTW_FORWARD, FFTW_ESTIMATE)
         CALL dfftw_execute_dft(plan, z_iso_resraman(0:md%t_cor*2, j), zhat_iso_resraman(0:md%t_cor*2, j)) !!!important to specify arrays!!
         CALL dfftw_destroy_plan(plan)

         CALL dfftw_plan_dft_1d(plan, 2*md%t_cor, z_aniso_resraman(0:md%t_cor*2, j), zhat_aniso_resraman(0:md%t_cor*2, j), &
                                FFTW_FORWARD, FFTW_ESTIMATE)
         CALL dfftw_execute_dft(plan, z_aniso_resraman(0:md%t_cor*2, j), zhat_aniso_resraman(0:md%t_cor*2, j)) !!!important to specify arrays!!
         CALL dfftw_destroy_plan(plan)
      END DO

      freq_res = REAL(freq_range/(2*md%t_cor), kind=dp)
      j = ANINT(laser_in_resraman/rtp_freq_res, kind=dp)

!!!!UNPOLARIZED!!!!

!zhat_iso_resraman=AIMAG(zhat_iso_resraman)
!zhat_aniso_resraman=AIMAG(zhat_aniso_resraman)

!OPEN(UNIT=30,FILE='zhat_aimag.txt',STATUS='unknown',IOSTAT=stat)
!DO i=0,2*md%t_cor-2
! WRITE(30,*) REAL(zhat_aniso_resraman(i,j),kind=dp),AIMAG(zhat_aniso_resraman(i,j))
!ENDDO
!CLOSE(30)

      f = freq_res*dt*1.883652d-4
      OPEN (UNIT=73, FILE='o-NP_resraman.txt', STATUS='unknown', IOSTAT=stat)
      DO i = 0, 2*md%t_cor - 2
         ! j=22
         !zhat_iso_resraman(i+1,j),AIMAG(zhat_iso_resraman(i+1,j),kind=dp)*(f*(i+1)/SIN(f*(i+1)))**2._dp
         zhat_iso_resraman(i + 1, j) = (zhat_iso_resraman(i + 1, j))*(f*(i + 1)/SIN(f*(i + 1)))**2._dp
         zhat_aniso_resraman(i + 1, j) = (zhat_aniso_resraman(i + 1, j))*(f*(i + 1)/SIN(f*(i + 1)))**2._dp

         zhat_unpol_resraman(i, j) = (zhat_iso_resraman(i, j)) + (zhat_aniso_resraman(i, j)*7.0_dp/45.0_dp)* &
                                     ((const_planck)/(8.0_dp*const_boltz*const_permit*const_permit) &
                                      *1d-29*0.421_dp*dt*((((laser_in*j) - ((i)*freq_res))**4)/((i)*freq_res)) &
                                      *(1.0_dp/(1.0_dp - EXP((-1.438777_dp*((i)*freq_res))/gs%temp))))*2.0_dp*2.0_dp*pi
         zhat_unpol_resraman(0, j) = 0.0_dp
         IF ((i*freq_res) .GE. 5000.0_dp) CYCLE
         !WRITE(73,*) i*freq_res,REAL(zhat_unpol_resraman(i,j),kind=dp),j
         WRITE (73, *) i*freq_res, zhat_unpol_resraman(i, j), j
! ENDDO
      END DO
      CLOSE (73)

   END SUBROUTINE spec_resraman
END MODULE calc_spectra
