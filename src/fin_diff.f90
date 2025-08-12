MODULE fin_diff
    USE kinds, ONLY: dp
    USE constants, ONLY: bohr2ang, speed_light, fs2s, damping_constant, joule_unit, ev_unit, action_unit
    USE vib_types, ONLY: global_settings, systems, molecular_dynamics, static, dipoles, raman
    IMPLICIT NONE
    PUBLIC :: central_diff, forward_diff, finite_diff_static, finite_diff_static_resraman

CONTAINS
    SUBROUTINE central_diff(natom, shifts, diff, sys, md)

        TYPE(systems), INTENT(INOUT)        :: sys
        TYPE(molecular_dynamics), INTENT(INOUT)        :: md
        INTEGER, INTENT(INOUT)                                    :: natom
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE, INTENT(INOUT)  ::  shifts
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE, INTENT(OUT)    ::  diff

        INTEGER                                                  :: stat, i, j, k, m

        ALLOCATE (diff(sys%framecount, natom, 3))
        !ALLOCATE( md%v(sys%sys%sys%framecount,1:44,3)) !change for fragments

        DO j = 1, sys%framecount - 2
            DO i = 1, natom
                DO k = 1, 3
                    diff(j, i, k) = (shifts(j + 2, i, k) - shifts(j, i, k))/(2.0_dp*md%dt)
                END DO
            END DO
        END DO

    END SUBROUTINE central_diff
!**************************************************************************************************************!
!**************************************************************************************************************!
    SUBROUTINE forward_diff(mol_num, alpha, dip_free, dip_x, gs, sys, dips)
        TYPE(global_settings), INTENT(INOUT) :: gs
        TYPE(systems), INTENT(INOUT)        :: sys
        TYPE(dipoles), INTENT(INOUT)        :: dips
        INTEGER, INTENT(INOUT)                                    :: mol_num
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE, INTENT(INOUT)  :: dip_free, dip_x
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE, INTENT(OUT)    :: alpha

        INTEGER                                                  :: stat, i, j, k, m

        ALLOCATE (alpha(sys%framecount, mol_num, 3))
        !ALLOCATE(alpha(sys%framecount,8:37,3))

        IF (gs%spectral_type%read_function.NE.'MD-RR') THEN
            DO j = 1, sys%framecount
                DO i = 1, mol_num  !!! change to mol_num later
                    DO k = 1, 3
                        alpha(j, i, k) = REAL((dip_x(j, i, k) - dip_free(j, i, k))/dips%e_field, kind=dp)
                        !alpha_x(j,i,k)=(dip_x(j,i,k)-dip_free(j,i,k))
                    END DO
                END DO
            END DO

        ELSEIF (gs%spectral_type%read_function=='MD-RR') THEN
            DO j = 1, sys%framecount
                DO i = 2, mol_num
                    DO k = 1, 3
                        alpha(j, i - 1, k) = REAL((dip_x(j, i, k) - dip_x(j, 1, k)) &
                                                  *(EXP(-7.0_dp*1.0_dp*REAL(i/(mol_num - 1.0_dp), kind=dp))**2.0_dp)/0.005_dp, &
                                                  kind=dp)
                    END DO
                END DO
            END DO

            DO j = 1, sys%framecount
                DO i = 2, mol_num
                    DO k = 1, 3
                        alpha(j, i - 1, k) = alpha(j, i - 1, k)*0.5_dp*(1 + COS(2.0_dp*3.14_dp*i/(2.0_dp*(mol_num - 1))))
                    END DO
                END DO
            END DO
        END IF

    END SUBROUTINE forward_diff

!**************************************************************************************************************!
!**************************************************************************************************************!
    SUBROUTINE finite_diff_static(gs, sys, stats, dips, rams)

        TYPE(global_settings), INTENT(INOUT) :: gs
        TYPE(systems), INTENT(INOUT)   :: sys
        TYPE(static), INTENT(INOUT):: stats
        TYPE(dipoles), INTENT(INOUT)    ::  dips
        TYPE(raman), INTENT(INOUT)   :: rams

        INTEGER                                                      :: stat, i, j, k, m
        REAL(kind=dp)                                                 :: factor
        REAL(kind=dp), DIMENSION(:, :, :, :), ALLOCATABLE                  :: pol_dxyz
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE                  :: dip_dxyz

        factor = REAL(1.0_dp/(2.0_dp*stats%dx), kind=dp)

        IF (gs%spectral_type%read_function=='R') THEN
            ALLOCATE (pol_dxyz(sys%natom, 3, 3, 3))
            ALLOCATE (rams%pol_dq(stats%nmodes, 3, 3))

            rams%pol_dq = 0.0_dp

            IF (dips%type_dipole=='berry') THEN
                rams%pol(:, :, :, 1, :) = REAL((rams%static_dip_x(:, :, :, :) - dips%static_dip(:, :, :, :))/(5.338d-5*1.313d-26), kind=dp)
                rams%pol(:, :, :, 2, :) = REAL((rams%static_dip_y(:, :, :, :) - dips%static_dip(:, :, :, :))/(5.338d-5*1.313d-26), kind=dp)
                rams%pol(:, :, :, 3, :) = REAL((rams%static_dip_z(:, :, :, :) - dips%static_dip(:, :, :, :))/(5.338d-5*1.313d-26), kind=dp)

                DEALLOCATE (dips%static_dip, rams%static_dip_x, rams%static_dip_y, rams%static_dip_z)
            END IF

            pol_dxyz(:, :, :, :) = REAL((rams%pol(:, :, 1, :, :) - rams%pol(:, :, 2, :, :))*factor, kind=dp)

        ELSEIF (gs%spectral_type%read_function=='IR') THEN

            ALLOCATE (dip_dxyz(sys%natom, 3, 3))
            ALLOCATE (dips%dip_dq(stats%nmodes, 3))

            dips%dip_dq = 0.0_dp

            dip_dxyz(:, :, :) = REAL((dips%static_dip(:, :, 1, :) - dips%static_dip(:, :, 2, :))*factor, kind=dp)

            DEALLOCATE (dips%static_dip)
        END IF

        DO j = 1, stats%nmodes
            DO k = 1, sys%natom
                IF (gs%spectral_type%read_function=='R') THEN
                    rams%pol_dq(j, :, :) = rams%pol_dq(j, :, :) + (pol_dxyz(k, 1, :, :)*stats%disp(j, k, 1)*sys%atom_mass_inv_sqrt(k)) &
                                           + (pol_dxyz(k, 2, :, :)*stats%disp(j, k, 2)*sys%atom_mass_inv_sqrt(k)) + (pol_dxyz(k, 3, :, :)* &
                                                                                                                     stats%disp(j, k, 3)*sys%atom_mass_inv_sqrt(k))
                ELSEIF (gs%spectral_type%read_function=='IR') THEN
                    dips%dip_dq(j, :) = dips%dip_dq(j, :) + (dip_dxyz(k, 1, :)*stats%disp(j, k, 1)*sys%atom_mass_inv_sqrt(k)) &
                                        + (dip_dxyz(k, 2, :)*stats%disp(j, k, 2)*sys%atom_mass_inv_sqrt(k)) + (dip_dxyz(k, 3, :)* &
                                                                                                               stats%disp(j, k, 3)*sys%atom_mass_inv_sqrt(k))
                END IF
            END DO
        END DO

        IF (gs%spectral_type%read_function=='R') THEN
            DEALLOCATE (pol_dxyz, rams%pol)
        ELSEIF (gs%spectral_type%read_function=='IR') THEN
            DEALLOCATE (dip_dxyz)
        END IF

    END SUBROUTINE finite_diff_static

!**************************************************************************************************************!
!**************************************************************************************************************!
    SUBROUTINE finite_diff_static_resraman(static_dipole_x_rtp, static_dipole_y_rtp, static_dipole_z_rtp, sys, rams)

        TYPE(systems), INTENT(INOUT)        :: sys
        TYPE(raman), INTENT(INOUT)        :: rams
        REAL(kind=dp), DIMENSION(:, :, :, :, :), ALLOCATABLE, INTENT(INOUT)  :: static_dipole_x_rtp, static_dipole_y_rtp, static_dipole_z_rtp

        INTEGER                                                      ::  i, j, k, m, l
        REAL(kind=dp)                                                 :: damping_factor, conv_unit

        ALLOCATE (rams%RR%pol_rtp(sys%natom, 3, 2, 3, 3, rams%RR%framecount_rtp))

        conv_unit = damping_constant*joule_unit/ev_unit !! J
        damping_factor = conv_unit/action_unit*rams%RR%dt_rtp*fs2s !! s-1
                            
        DO l = 2, rams%RR%framecount_rtp + 1
            rams%RR%pol_rtp(:, :, :, 1, :, l - 1) = static_dipole_x_rtp(:, :, :, :, l) - static_dipole_x_rtp(:, :, :, :, 1)
            rams%RR%pol_rtp(:, :, :, 2, :, l - 1) = static_dipole_y_rtp(:, :, :, :, l) - static_dipole_y_rtp(:, :, :, :, 1)
            rams%RR%pol_rtp(:, :, :, 3, :, l - 1) = static_dipole_z_rtp(:, :, :, :, l) - static_dipole_z_rtp(:, :, :, :, 1)
        ENDDO

      DO l=1,rams%RR%framecount_rtp
        rams%RR%pol_rtp(:,:,:,:,:,l)= rams%RR%pol_rtp(:,:,:,:,:,l)*(EXP(-1.0d0*damping_factor*l))
      ENDDO

       ! DO i = 1, framecount_rtp
        !    pol_rtp(:,:,:,:,:,i) = pol_rtp(:,:,:,:,:,i)*((COS(i/(framecount_rtp - 1.0_dp)/2.0_dp*3.14_dp))**2) !!Hann Window function
        !END DO
    END SUBROUTINE finite_diff_static_resraman
END MODULE fin_diff

