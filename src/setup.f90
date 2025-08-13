MODULE setup

   USE kinds, ONLY: dp
   USE constants, ONLY: speed_light
   USE vib_types, ONLY: global_settings, systems

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: read_input, masses_charges, conversion, pbc_orthorombic, pbc_hexagonal, pbc_hexagonal_old, pbc_orthorombic_old !constants,

CONTAINS

!*************************************************************************************************
!*************************************************************************************************

   SUBROUTINE read_input(filename, static_pol_file, static_dip_free_file, static_dip_x_file, static_dip_y_file, &
                         static_dip_z_file, normal_freq_file, normal_displ_file, read_function, system, &
                         length, box_all, box_x, box_y, box_z, dt, type_input, wannier_free, wannier_x, wannier_y, wannier_z, &
                         input_mass, periodic, direction, averaging, type_dipole, cell_type, rtp_dipole_x, rtp_dipole_y, &
                         rtp_dipole_z, framecount_rtp, dt_rtp, frag_type, type_static, force_file, &
                         laser_in, check_pade, dx, framecount_rtp_pade)

      CHARACTER(LEN=40), INTENT(OUT)              :: filename, static_pol_file, read_function, length, system, type_input, periodic
      CHARACTER(LEN=40), INTENT(OUT)              :: wannier_free, wannier_x, wannier_y, wannier_z, input_mass, type_static
      CHARACTER(LEN=40), INTENT(OUT)              :: direction, averaging, type_dipole, cell_type, force_file
      CHARACTER(LEN=40), INTENT(OUT)              :: static_dip_free_file, static_dip_x_file, static_dip_y_file, static_dip_z_file
      CHARACTER(LEN=40), INTENT(OUT)              :: normal_freq_file, normal_displ_file, frag_type
      CHARACTER(LEN=40), INTENT(OUT)              :: rtp_dipole_x, rtp_dipole_y, rtp_dipole_z, check_pade
      REAL(kind=dp), INTENT(OUT)                  :: dt, dt_rtp, box_all, box_x, box_y, box_z, laser_in, dx
      INTEGER, INTENT(OUT)                        :: framecount_rtp, framecount_rtp_pade

      DO
         WRITE (*, *) 'Enter which function you want to calculate (type "P" for Power spectrum , "MD-IR" for MD-based IR spectrum, &
        &         "MD-R" for MD-based Raman spectrum, "MD-RR" for MD-based resonance Raman, "NMA" for normal mode analysis, &
        &         "IR" for static IR spectrum, "R" for static Raman spectrum, "ABS" for absorption spectrum,"RR" for static &
        &         resonance Raman spectrum)'
         READ (*, *) read_function
      IF (read_function .NE. 'P' .AND. read_function .NE. 'MD-IR' .AND. read_function .NE. 'MD-R' .AND. read_function .NE. 'MD-RR' &
    .AND. read_function .NE. 'NMA' .AND. read_function .NE. 'IR' .AND. read_function .NE. 'R' .AND. read_function .NE. 'ABS' .AND. &
             read_function .NE. 'RR') THEN
            WRITE (*, *) 'Please type P, MD-IR, MD-R, MD-RR, NMA, IR, R, ABS or RR!'
            CYCLE
         END IF
         EXIT
      END DO

!read_function='RR'

      DO
         IF (read_function == 'P') THEN
            WRITE (*, *) 'Enter the type of the input file (type 1 for positions, 2 for velocities)'
            READ (*, *) type_input
            IF (type_input .NE. '1' .AND. type_input .NE. '2') THEN
               WRITE (*, *) 'Please type 1 or 2!'
               CYCLE
            END IF
         END IF
         EXIT
      END DO

      DO
         IF (read_function == 'P') THEN
            WRITE (*, *) 'Do you want to apply mass weighting (y/n)?'
            READ (*, *) input_mass
            IF (input_mass .NE. 'y' .AND. input_mass .NE. 'n') THEN
               WRITE (*, *) 'Please type y or n!'
               CYCLE
            END IF
         END IF
         EXIT
      END DO

      DO
         IF (read_function == 'IR' .OR. read_function == 'R' .OR. read_function == 'RR') THEN
            WRITE (*, *) 'Do you want the normal modes to be calculated (type "1") or read from an external file (type "2")?'
            READ (*, *) type_static
            IF (type_static .NE. '1' .AND. type_static .NE. '2') THEN
               WRITE (*, *) 'Please type 1 or 2!'
               CYCLE
            END IF
         END IF
         EXIT
      END DO

      DO
         IF (read_function == 'MD-R' .OR. read_function == 'R' .OR. read_function == 'IR' .OR. read_function == 'RR' .OR. &
             read_function == 'ABS') THEN
            WRITE (*, *) 'Which one do you want to use: Wannier centers (1), Berry phase dipole moments (2)&
    &                    or DFPT polarizabilities (3)?'
            READ (*, *) type_dipole
            IF (type_dipole .NE. '1' .AND. type_dipole .NE. '2' .AND. type_dipole .NE. '3') THEN
               WRITE (*, *) 'Please type 1, 2 or 3!!'
               CYCLE
            END IF
         END IF
         EXIT
      END DO

      DO
         IF (read_function == 'MD-IR') THEN
            WRITE (*, *) 'Which one do you want to use: Wannier centers (1) or Berry phase dipole moments (2)?'
            READ (*, *) type_dipole
            IF (type_dipole .NE. '1' .AND. type_dipole .NE. '2') THEN
               WRITE (*, *) 'Please type 1 or 2!!'
               CYCLE
            END IF
         END IF
         EXIT
      END DO

      DO
         IF (read_function == 'MD-IR' .OR. read_function == 'MD-R') THEN
            WRITE (*, *) 'Do you want to apply the fragment approach (1) or the molecular approach? (2)'
            READ (*, *) system
            IF (system .NE. '1' .AND. system .NE. '2') THEN
               WRITE (*, *) 'Please type 1 or 2!'
               CYCLE
            END IF
         END IF
         EXIT
      END DO

      DO
         IF (read_function == 'MD-IR' .OR. read_function == 'MD-R') THEN
            WRITE (*, *) 'Is it the k-point trajectory (1), supercell trajectory (2) or solvent trajectory (3)? '
            READ (*, *) cell_type
            IF (cell_type .NE. '1' .AND. cell_type .NE. '2' .AND. cell_type .NE. '3') THEN
               WRITE (*, *) 'Please type 1, 2 or 3!'
               CYCLE
            END IF
         END IF
         EXIT
      END DO

      DO
         IF (read_function == 'MD-R') THEN
            IF (type_dipole == '3') THEN
               IF (cell_type == '1') THEN
                  wannier_free = 'alpha_x_kp_5000.xyz'
                  wannier_x = 'alpha_x_kp_5000.xyz'
                  wannier_y = 'alpha_y_kp_5000.xyz'
                  wannier_z = 'alpha_z_kp_5000.xyz'
               ELSEIF (cell_type == '2') THEN
                  ! wannier_free='alpha_x_ismail.xyz'
                  ! wannier_x='alpha_x_ismail.xyz'
                  ! wannier_y='alpha_y_ismail.xyz'
                  ! wannier_z='alpha_z_ismail.xyz'
                  !  wannier_free='alpha_x_ismail_400.xyz'
                  !  wannier_x='alpha_x_ismail_400.xyz'
                  !  wannier_y='alpha_y_ismail_400.xyz'
                  !  wannier_z='alpha_z_ismail_400.xyz'
                  !  wannier_free='alpha_x.xyz'
                  ! wannier_x='alpha_x.xyz'
                  ! wannier_y='alpha_y.xyz'
                  ! wannier_z='alpha_z.xyz'
                  wannier_free = 'alpha_x_solv_5000_new.xyz'
                  wannier_x = 'alpha_x_solv_5000_new.xyz'
                  wannier_y = 'alpha_y_solv_5000_new.xyz'
                  wannier_z = 'alpha_z_solv_5000_new.xyz'
               END IF
            ELSEIF (type_dipole == '2') THEN
               IF (cell_type == '1') THEN
                  wannier_free = 'berry_dipole_free_5000_kp.xyz'
                  wannier_x = 'berry_dipole_x_5000_kp_large.xyz'
                  wannier_y = 'berry_dipole_y_5000_kp_large.xyz'
                  wannier_z = 'berry_dipole_z_5000_kp_large.xyz'
               ELSEIF (cell_type == '2') THEN
                  !  wannier_free='dipole_o-NP_free.xyz'
                  ! wannier_x='dipole_o-NP_X.xyz'
                  ! wannier_y='dipole_o-NP_Y.xyz'
                  ! wannier_z='dipole_o-NP_Z.xyz'
                  !  wannier_x='dipole_o-NP_X_smallfield.xyz'
                  !  wannier_y='dipole_o-NP_Y_smallfield.xyz'
                  !  wannier_z='dipole_o-NP_Z_smallfield.xyz'
                  wannier_free = 'berry_dipole_free_5000_sc.xyz'
                  wannier_x = 'berry_dipole_X_5000_sc.xyz'
                  wannier_y = 'berry_dipole_Y_5000_sc.xyz'
                  wannier_z = 'berry_dipole_Z_5000_sc.xyz'
               END IF
            ELSEIF (type_dipole == '1') THEN
               IF (cell_type == '1') THEN
                  wannier_free = 'wannier_free_COF-1_kp_5000.xyz'
                  wannier_x = 'wannier_X_COF-1_kp_large_5000.xyz'
                  wannier_y = 'wannier_Y_COF-1_kp_large_5000.xyz'
                  wannier_z = 'wannier_Z_COF-1_kp_large_5000.xyz'
               ELSEIF (cell_type == '2') THEN
                  wannier_free = 'wannier_free_COF-1_sc_5000.xyz'
                  wannier_x = 'wannier_X_COF-1_sc_5000.xyz'
                  wannier_y = 'wannier_Y_COF-1_sc_5000.xyz'
                  wannier_z = 'wannier_Z_COF-1_sc_5000.xyz'
               ELSEIF (cell_type == '3') THEN
                  wannier_free = 'COF-1_solv_wannier_free.xyz'
                  wannier_x = 'COF-1_solv_wannier_X.xyz'
                  wannier_y = 'COF-1_solv_wannier_Y.xyz'
                  wannier_z = 'COF-1_solv_wannier_Z.xyz'
               END IF
               !      wannier_free='wannier_free.xyz'
               !     wannier_x='wannier_x.xyz'
               !    wannier_y='wannier_y.xyz'
               !   wannier_z='wannier_z.xyz'

               !  wannier_free='new_wan.xyz'
               ! wannier_x='polarizability_xp-wannier.xyz'
               ! wannier_y='polarizability_yp-wannier.xyz'
               ! wannier_z='polarizability_zp-wannier.xyz'

               ! WRITE(*,*) 'Enter the field free dipole moment data'
               ! READ(*,*) wannier_free
               ! WRITE(*,*) 'Enter the X-field dipole moment data'
               ! READ(*,*) wannier_x
               ! WRITE(*,*) 'Enter the Y-field dipole moment data'
               ! READ(*,*) wannier_y
               ! WRITE(*,*) 'Enter the Z-field dipole moment data'
               ! READ(*,*) wannier_z
            END IF
         END IF
         EXIT
      END DO

      DO
         IF (read_function == 'MD-R') THEN
            WRITE (*, *) 'What is the type of averaging you want to apply? (type 1 for orientational, 2 for isotropic)'
            READ (*, *) averaging
            IF (averaging .NE. '1' .AND. averaging .NE. '2') THEN
               WRITE (*, *) 'Please type 1 or 2!'
               CYCLE
            END IF
         END IF
         EXIT
      END DO

      DO
         IF (read_function == 'MD-R' .AND. averaging == '2') THEN
            WRITE (*, *) 'What is the direction of the applied electric field (type 1 for x, 2 for y, 3 for z)?'
            READ (*, *) direction
            IF (direction .NE. '1' .AND. direction .NE. '2' .AND. direction .NE. '3') THEN
               WRITE (*, *) 'Please type 1, 2 or 3!'
               CYCLE
            END IF
         END IF
         EXIT
      END DO

      DO
         IF (read_function == 'MD-R') THEN
            !IF (read_function=='MD-RR' .OR. read_function=='RR') THEN
            WRITE (*, *) 'What is the wavenumber (cm^-1) of the incident laser?'
            READ (*, *) laser_in
         END IF
         EXIT
      END DO

      DO
         IF (read_function == 'MD-IR' .OR. read_function == 'MD-R') THEN
            WRITE (*, *) 'Does the system contain more than one molecule? (y/n)'
            READ (*, *) periodic
            IF (periodic .NE. 'y' .AND. periodic .NE. 'n') THEN
               WRITE (*, *) 'Please type y or n!'
               CYCLE
            END IF
         END IF
         EXIT
      END DO

      DO
         IF (system == '1') THEN
            WRITE (*, *) 'Which fragments do you want to calculate? (BO (1), Ph (2) or B-C bonds (3)?'
            READ (*, *) frag_type
            IF (frag_type .NE. '1' .AND. frag_type .NE. '2' .AND. frag_type .NE. '3') THEN
               WRITE (*, *) 'Please type 1, 2 or 3!'
               CYCLE
            END IF
         END IF
         EXIT
      END DO

      DO
         IF (read_function == 'NMA') THEN
            WRITE (*, *) 'Enter the name of the geometry file:'
            READ (*, *) filename
            WRITE (*, *) 'Enter the name of the force file:'
            READ (*, *) force_file
            WRITE (*, *) 'What is the normal mode displacement in Angstrom?'
            READ (*, *) dx
            !force_file='water-force.data3'
            !dx=0.001d0
         END IF
         EXIT
      END DO

      DO
         IF (read_function == 'IR' .OR. read_function == 'R' .OR. read_function == 'RR' .OR. read_function == 'ABS') THEN
            WRITE (*, *) 'Enter the name of the geometry file:'
            READ (*, *) filename
            IF (read_function .NE. 'ABS') THEN
               IF (type_static == '1') THEN
                  WRITE (*, *) 'Enter the name of the force file:'
                  READ (*, *) force_file
               ELSEIF (type_static == '2') THEN
                  WRITE (*, *) 'Enter the name of the normal mode frequency file:'
                  READ (*, *) normal_freq_file
                  WRITE (*, *) 'Enter the name of the normal mode coordinates file:'
                  READ (*, *) normal_displ_file
               END IF
               WRITE (*, *) 'What is the normal mode displacement in Angstrom?'
               READ (*, *) dx
            END IF
            IF (type_dipole == '3') THEN
               WRITE (*, *) 'Enter the name of the file that contains polarizabilities:'
               READ (*, *) static_pol_file
            ELSEIF (type_dipole == '2') THEN
               IF (read_function == 'IR' .OR. read_function == 'R') THEN
                  WRITE (*, *) 'Enter the name of the static dipoles file (Field-free):'
                  READ (*, *) static_dip_free_file
               END IF
               IF (read_function .NE. 'IR') THEN
                  WRITE (*, *) 'Enter the name of the static dipoles file (X-Field):'
                  READ (*, *) static_dip_x_file
                  WRITE (*, *) 'Enter the name of the static dipoles file (Y-Field):'
                  READ (*, *) static_dip_y_file
                  WRITE (*, *) 'Enter the name of the static dipoles file (Z-Field):'
                  READ (*, *) static_dip_z_file
               END IF
            END IF
            IF (read_function == 'RR' .OR. read_function == 'ABS') THEN
               WRITE (*, *) 'Enter the RT-TDDFT time step:'
               READ (*, *) dt_rtp
               WRITE (*, *) 'Enter the total number of RT-TDDFT steps:'
               READ (*, *) framecount_rtp
               WRITE (*, *) 'Do you want to apply Pade approximants? (y/n)'
               READ (*, *) check_pade
               IF (check_pade == 'y') THEN
                  WRITE (*, *) 'Enter the final number of RT-TDDFT steps after the application of Pade interpolation:'
                  READ (*, *) framecount_rtp_pade
               END IF
            END IF
            IF (read_function == 'R' .OR. read_function == 'RR') THEN
               WRITE (*, *) 'What is the wavenumber (cm^-1) of the incident laser?'
               READ (*, *) laser_in
            END IF
         END IF
         EXIT
      END DO

      DO
         IF (read_function == 'MD-RR') THEN
            rtp_dipole_x = 'o-NP_RTP_dipoles_X_256.xyz'
            rtp_dipole_y = 'o-NP_RTP_dipoles_Y_256.xyz'
            rtp_dipole_z = 'o-NP_RTP_dipoles_Z_256.xyz'
            ! rtp_dipole_x='o-NP_RTP_dipoles_X.xyz'
            ! rtp_dipole_y='o-NP_RTP_dipoles_Y.xyz'
            rtp_dipole_z = 'o-NP_RTP_dipoles_Z.xyz'
            ! WRITE(*,*)'What is the number of RTP frames?'
            ! READ(*,*) framecount_rtp
            ! framecount_rtp=1280
            framecount_rtp = 256
            ! dt_rtp=0.0125_dp
            dt_rtp = 0.0625_dp
            ! WRITE (*, *) 'What is the wavenumber of the incident laser (cm^-1)?'
            ! READ (*, *) laser_in_resraman
            !laser_in_resraman=15797.788309636651_dp
         END IF
         EXIT
      END DO

      ! laser_in_resraman=15797.788309636651_dp
      !laser_in_resraman=15808.596424_dp !r-met NR
      !laser_in_resraman = 57346.490087_dp !r-met RR
      ! laser_in_resraman=41860.518081_dpi
      ! check_pade='n'
      ! framecount_rtp=80000
      ! check_pade='n'

      DO
         IF (read_function == 'P' .OR. read_function == 'MD-IR') THEN
            WRITE (*, *) 'Enter the name of the trajectory'
            READ (*, *) filename
         END IF
         IF (read_function == 'P' .OR. read_function == 'MD-IR' .OR. read_function == 'MD-R') THEN
            WRITE (*, *) 'Enter the time step (fs)'
            READ (*, *) dt
            IF (read_function .NE. 'P') THEN
               WRITE (*, *) 'Are the 3 cell vectors the same length? (y/n)'
               READ (*, *) length
               IF (length == 'y') THEN
                  WRITE (*, *) 'Enter the cell vector (in Armstrong)'
                  READ (*, *) box_all
               ELSEIF (length == 'n') THEN
                  IF (cell_type == '1') THEN
                     box_x = 15.100_dp
                     box_y = 15.100_dp
                     box_z = 13.448_dp
                  ELSEIF (cell_type == '2') THEN
                     box_x = 15.100_dp
                     box_y = 15.101_dp
                     box_z = 13.457_dp
                  ELSEIF (cell_type == '3') THEN
                     box_x = 15.100_dp
                     box_y = 15.100_dp
                     box_z = 20.172_dp
                  END IF
                  !         WRITE(*,*)'Enter the cell vector for X direction (in Armstrong)'
                  !         READ(*,*) box_x
                  !        WRITE(*,*)'Enter the cell vector for Y direction (in Armstrong)'
                  !       READ(*,*) box_y
                  !      WRITE(*,*)'Enter the cell vector for Z direction (in Armstrong)'
                  !     READ(*,*) box_z
               END IF
            END IF
         END IF
         EXIT
      END DO

      IF (length == 'y') THEN
         box_x = box_all
         box_y = box_all
         box_z = box_all
      END IF

   END SUBROUTINE read_input

!*********************************************************************************************
!*********************************************************************************************
   SUBROUTINE masses_charges(gs, sys)!natom, mass_atom, atom_mass_inv_sqrt, mass_mat, element, mass_tot, charge)

      TYPE(global_settings), INTENT(INOUT)   :: gs
      TYPE(systems), INTENT(INOUT)        :: sys
      !CHARACTER(LEN=2), DIMENSION(:), ALLOCATABLE, INTENT(IN)   :: element
      !INTEGER, INTENT(INOUT)                                  :: natom
      !REAL(kind=dp), INTENT(INOUT)                             :: mass_tot
      !REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(OUT)      :: atom_mass_inv_sqrt, mass_atom, charge
      !REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE, INTENT(OUT)    :: mass_mat

      INTEGER                                                :: i, j, stat
      REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE                :: mat1, mat2

   ALLOCATE (sys%atom_mass_inv_sqrt(sys%natom), sys%mass_mat(sys%natom, sys%natom), sys%mass_atom(sys%natom), sys%charge(sys%natom))
      ALLOCATE (mat1(sys%natom, 1), mat2(1, sys%natom))

      sys%mass_atom = 0.0_dp
      sys%mass_tot = 0.0_dp
      DO i = 1, sys%natom
         IF (sys%element(i) == 'O') THEN
            sys%mass_atom(i) = 15.999_dp
            sys%charge(i) = 6.0_dp
         ELSEIF (sys%element(i) == 'H') THEN
            sys%mass_atom(i) = 1.00784_dp
            sys%charge(i) = 1.0_dp
         ELSEIF (sys%element(i) == 'C') THEN
            sys%mass_atom(i) = 12.011_dp
            sys%charge(i) = 4.0_dp
         ELSEIF (sys%element(i) == 'B') THEN
            sys%mass_atom(i) = 10.811_dp
            sys%charge(i) = 3.0_dp
         ELSEIF (sys%element(i) == 'N') THEN
            sys%mass_atom(i) = 14.0067_dp
            sys%charge(i) = 5.0_dp
         ELSEIF (sys%element(i) == 'Fe') THEN
            sys%mass_atom(i) = 55.845_dp
            sys%charge(i) = 2.0_dp
         ELSEIF (sys%element(i) == 'X') THEN
            sys%mass_atom(i) = 0.00_dp
            sys%charge(i) = -2.0_dp
         END IF
         sys%mass_tot = sys%mass_tot + sys%mass_atom(i)
      END DO

      sys%atom_mass_inv_sqrt(:) = SQRT(REAL(1.0_dp/sys%mass_atom(:), kind=dp))

      mat1(:, :) = RESHAPE(sys%atom_mass_inv_sqrt(:), (/sys%natom, 1/))
      mat2(:, :) = RESHAPE(sys%atom_mass_inv_sqrt(:), (/1, sys%natom/))

      sys%mass_mat = MATMUL(mat1, mat2)

      DEALLOCATE (mat1, mat2)

   END SUBROUTINE masses_charges

!!*********************************************************************************************
!!*********************************************************************************************

   SUBROUTINE conversion(dt, dom, dt_rtp, dom_rtp, freq_range, sinc_const)

      REAL(kind=dp), INTENT(OUT)              :: dom, dom_rtp, freq_range, sinc_const
      REAL(kind=dp), INTENT(IN)               ::  dt, dt_rtp

      INTEGER                               :: stat   ! error status of OPEN statements
      INTEGER                               :: i, j, k

      dom = REAL((1.0_dp/(dt*1e-15))/speed_light, kind=dp)
      dom_rtp = REAL((1.0_dp/(dt_rtp*1e-15))/speed_light, kind=dp)

      ! freq_range = REAL(dom/(2.0_dp*md%t_cor), kind=dp)
      !sinc_const = freq_range*dt*1.883652d-4 !!for sinc function

   END SUBROUTINE conversion

!*********************************************************************************************
!*********************************************************************************************

   SUBROUTINE pbc_orthorombic(coord2, coord1, sys)

      TYPE(systems), INTENT(INOUT)                :: sys
      REAL(kind=dp), DIMENSION(3), INTENT(INOUT)                ::coord2, coord1

      sys%cell%vec(:) = coord2(:) - coord1(:)

      sys%cell%vec_pbc(1) = sys%cell%vec(1) - sys%cell%box_x*ANINT((1./sys%cell%box_x)*sys%cell%vec(1))
      sys%cell%vec_pbc(2) = sys%cell%vec(2) - sys%cell%box_y*ANINT((1./sys%cell%box_y)*sys%cell%vec(2))
      sys%cell%vec_pbc(3) = sys%cell%vec(3) - sys%cell%box_z*ANINT((1./sys%cell%box_z)*sys%cell%vec(3))

   END SUBROUTINE pbc_orthorombic

!********************************************************************************************
!********************************************************************************************

   SUBROUTINE pbc_hexagonal(coord2, coord1, sys)

      TYPE(systems), INTENT(INOUT)                :: sys
      REAL(kind=dp), DIMENSION(3), INTENT(INOUT)                :: coord2, coord1
      REAL(kind=dp)                                           :: h_inv(3, 3), a, s(3), hmat(3, 3)
      REAL(kind=dp)                                           :: acosa, asina, sqrt3, det_a

      sqrt3 = 1.73205080756887729352744634_dp

      a = 0.5_dp*(sys%cell%box_x + sys%cell%box_y)
      acosa = 0.5_dp*a
      asina = sqrt3*acosa
      hmat(1, 1) = a; hmat(1, 2) = acosa; hmat(1, 3) = 0.0_dp
      hmat(2, 1) = 0.0_dp; hmat(2, 2) = asina; hmat(2, 3) = 0.0_dp
      hmat(3, 1) = 0.0_dp; hmat(3, 2) = 0.0_dp; hmat(3, 3) = sys%cell%box_z

      det_a = hmat(1, 1)*(hmat(2, 2)*hmat(3, 3) - hmat(2, 3)*hmat(3, 2)) - &
              hmat(1, 2)*(hmat(2, 3)*hmat(3, 1) - hmat(2, 1)*hmat(3, 3)) + &
              hmat(1, 3)*(hmat(2, 1)*hmat(3, 2) - hmat(2, 2)*hmat(3, 1))

      det_a = 1./det_a

      h_inv(1, 1) = (hmat(2, 2)*hmat(3, 3) - hmat(3, 2)*hmat(2, 3))*det_a
      h_inv(2, 1) = (hmat(2, 3)*hmat(3, 1) - hmat(3, 3)*hmat(2, 1))*det_a
      h_inv(3, 1) = (hmat(2, 1)*hmat(3, 2) - hmat(3, 1)*hmat(2, 2))*det_a

      h_inv(1, 2) = (hmat(1, 3)*hmat(3, 2) - hmat(3, 3)*hmat(1, 2))*det_a
      h_inv(2, 2) = (hmat(1, 1)*hmat(3, 3) - hmat(3, 1)*hmat(1, 3))*det_a
      h_inv(3, 2) = (hmat(1, 2)*hmat(3, 1) - hmat(3, 2)*hmat(1, 1))*det_a

      h_inv(1, 3) = (hmat(1, 2)*hmat(2, 3) - hmat(2, 2)*hmat(1, 3))*det_a
      h_inv(2, 3) = (hmat(1, 3)*hmat(2, 1) - hmat(2, 3)*hmat(1, 1))*det_a
      h_inv(3, 3) = (hmat(1, 1)*hmat(2, 2) - hmat(2, 1)*hmat(1, 2))*det_a

      sys%cell%vec(:) = coord2(:) - coord1(:)

      s(1) = h_inv(1, 1)*sys%cell%vec(1) + h_inv(1, 2)*sys%cell%vec(2) + h_inv(1, 3)*sys%cell%vec(3)
      s(2) = h_inv(2, 1)*sys%cell%vec(1) + h_inv(2, 2)*sys%cell%vec(2) + h_inv(2, 3)*sys%cell%vec(3)
      s(3) = h_inv(3, 1)*sys%cell%vec(1) + h_inv(3, 2)*sys%cell%vec(2) + h_inv(3, 3)*sys%cell%vec(3)

      s(1) = s(1) - ANINT(s(1))
      s(2) = s(2) - ANINT(s(2))
      s(3) = s(3) - ANINT(s(3))

      sys%cell%vec_pbc(1) = hmat(1, 1)*s(1) + hmat(1, 2)*s(2) + hmat(1, 3)*s(3)
      sys%cell%vec_pbc(2) = hmat(2, 1)*s(1) + hmat(2, 2)*s(2) + hmat(2, 3)*s(3)
      sys%cell%vec_pbc(3) = hmat(3, 1)*s(1) + hmat(3, 2)*s(2) + hmat(3, 3)*s(3)

   END SUBROUTINE pbc_hexagonal

   SUBROUTINE pbc_hexagonal_old(coord2, coord1, vec, vec_pbc, box_all, box_x, box_y, box_z)

      REAL(kind=dp), DIMENSION(3), INTENT(INOUT)                :: vec, vec_pbc, coord2, coord1
      REAL(kind=dp), INTENT(IN)                                :: box_all, box_x, box_y, box_z

      REAL(kind=dp)                                           :: h_inv(3, 3), a, s(3), hmat(3, 3)
      REAL(kind=dp)                                           :: acosa, asina, sqrt3, det_a

      sqrt3 = 1.73205080756887729352744634_dp

      a = 0.5_dp*(box_x + box_y)
      acosa = 0.5_dp*a
      asina = sqrt3*acosa
      hmat(1, 1) = a; hmat(1, 2) = acosa; hmat(1, 3) = 0.0_dp
      hmat(2, 1) = 0.0_dp; hmat(2, 2) = asina; hmat(2, 3) = 0.0_dp
      hmat(3, 1) = 0.0_dp; hmat(3, 2) = 0.0_dp; hmat(3, 3) = box_z

      det_a = hmat(1, 1)*(hmat(2, 2)*hmat(3, 3) - hmat(2, 3)*hmat(3, 2)) - &
              hmat(1, 2)*(hmat(2, 3)*hmat(3, 1) - hmat(2, 1)*hmat(3, 3)) + &
              hmat(1, 3)*(hmat(2, 1)*hmat(3, 2) - hmat(2, 2)*hmat(3, 1))

      det_a = 1./det_a

      h_inv(1, 1) = (hmat(2, 2)*hmat(3, 3) - hmat(3, 2)*hmat(2, 3))*det_a
      h_inv(2, 1) = (hmat(2, 3)*hmat(3, 1) - hmat(3, 3)*hmat(2, 1))*det_a
      h_inv(3, 1) = (hmat(2, 1)*hmat(3, 2) - hmat(3, 1)*hmat(2, 2))*det_a

      h_inv(1, 2) = (hmat(1, 3)*hmat(3, 2) - hmat(3, 3)*hmat(1, 2))*det_a
      h_inv(2, 2) = (hmat(1, 1)*hmat(3, 3) - hmat(3, 1)*hmat(1, 3))*det_a
      h_inv(3, 2) = (hmat(1, 2)*hmat(3, 1) - hmat(3, 2)*hmat(1, 1))*det_a

      h_inv(1, 3) = (hmat(1, 2)*hmat(2, 3) - hmat(2, 2)*hmat(1, 3))*det_a
      h_inv(2, 3) = (hmat(1, 3)*hmat(2, 1) - hmat(2, 3)*hmat(1, 1))*det_a
      h_inv(3, 3) = (hmat(1, 1)*hmat(2, 2) - hmat(2, 1)*hmat(1, 2))*det_a

      vec(:) = coord2(:) - coord1(:)

      s(1) = h_inv(1, 1)*vec(1) + h_inv(1, 2)*vec(2) + h_inv(1, 3)*vec(3)
      s(2) = h_inv(2, 1)*vec(1) + h_inv(2, 2)*vec(2) + h_inv(2, 3)*vec(3)
      s(3) = h_inv(3, 1)*vec(1) + h_inv(3, 2)*vec(2) + h_inv(3, 3)*vec(3)

      s(1) = s(1) - ANINT(s(1))
      s(2) = s(2) - ANINT(s(2))
      s(3) = s(3) - ANINT(s(3))

      vec_pbc(1) = hmat(1, 1)*s(1) + hmat(1, 2)*s(2) + hmat(1, 3)*s(3)
      vec_pbc(2) = hmat(2, 1)*s(1) + hmat(2, 2)*s(2) + hmat(2, 3)*s(3)
      vec_pbc(3) = hmat(3, 1)*s(1) + hmat(3, 2)*s(2) + hmat(3, 3)*s(3)

   END SUBROUTINE pbc_hexagonal_old
   SUBROUTINE pbc_orthorombic_old(coord2, coord1, vec, vec_pbc, box_all, box_x, box_y, box_z)

      REAL(kind=dp), DIMENSION(3), INTENT(INOUT)                :: vec, vec_pbc, coord2, coord1
      REAL(kind=dp), INTENT(IN)                                :: box_all, box_x, box_y, box_z

      vec(:) = coord2(:) - coord1(:)

      vec_pbc(1) = vec(1) - box_x*ANINT((1./box_x)*vec(1))
      vec_pbc(2) = vec(2) - box_y*ANINT((1./box_y)*vec(2))
      vec_pbc(3) = vec(3) - box_z*ANINT((1./box_z)*vec(3))

   END SUBROUTINE pbc_orthorombic_old
END MODULE setup

