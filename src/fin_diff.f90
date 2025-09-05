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
    SUBROUTINE forward_diff(mol_num, alpha, dip_free, dip_x, gs, sys)
        TYPE(global_settings), INTENT(INOUT) :: gs
        TYPE(systems), INTENT(INOUT)        :: sys
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
                        alpha(j, i, k) = REAL((dip_x(j, i, k) - dip_free(j, i, k))/0.005_dp, kind=dp)
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
        IMPLICIT NONE
        TYPE(global_settings), INTENT(INOUT) :: gs
        TYPE(systems), INTENT(INOUT) :: sys
        TYPE(static), INTENT(INOUT) :: stats
        TYPE(dipoles), INTENT(INOUT) :: dips
        TYPE(raman), INTENT(INOUT) :: rams

        INTEGER :: i_pol, j_pol, k, l, j, xyz
        REAL(kind=dp) :: factor, delta, coeff

        factor = 1.0_dp/(2.0_dp*stats%dx)

        IF (gs%spectral_type%read_function=='R') THEN
            ALLOCATE (rams%pol_dq(stats%nmodes, 3, 3))
            rams%pol_dq = 0.0_dp

            DO i_pol = 1, 3
                DO j_pol = 1, 3
                    DO k = 1, sys%natom
                        DO l = 1, 3
                            coeff = factor*sys%atom_mass_inv_sqrt(k)
                            delta = rams%pol(i_pol, j_pol)%atom(k)%displacement(1)%XYZ(l)%frame(1) &
                                    - rams%pol(i_pol, j_pol)%atom(k)%displacement(2)%XYZ(l)%frame(1)

                            DO j = 1, stats%nmodes
                                rams%pol_dq(j, i_pol, j_pol) = rams%pol_dq(j, i_pol, j_pol) &
                                                               + coeff*stats%disp(j, k, l)*delta
                            END DO
                        END DO
                    END DO
                END DO
            END DO

        ELSEIF (gs%spectral_type%read_function=='IR') THEN
            ALLOCATE (dips%dip_dq(stats%nmodes, 3))
            dips%dip_dq = 0.0_dp

            DO xyz = 1, 3
                DO k = 1, sys%natom
                    DO l = 1, 3
                        coeff = factor*sys%atom_mass_inv_sqrt(k)
                        delta = dips%static_dip(xyz)%atom(k)%displacement(2)%XYZ(l)%frame(1) &
                                - dips%static_dip(xyz)%atom(k)%displacement(1)%XYZ(l)%frame(1)

                        DO j = 1, stats%nmodes
                            dips%dip_dq(j, xyz) = dips%dip_dq(j, xyz) + coeff*stats%disp(j, k, l)*delta
                        END DO
                    END DO
                END DO
            END DO
        END IF

    END SUBROUTINE finite_diff_static

!**************************************************************************************************************!
!**************************************************************************************************************!
    SUBROUTINE finite_diff_static_resraman(sys, rams)

        TYPE(raman), INTENT(INOUT)        :: rams
        TYPE(systems), INTENT(INOUT)        :: sys

        INTEGER                                                      ::  x, y, i, j, k, m, l, i_pol, j_pol, stat
        REAL(kind=dp)                                                 :: damping_factor, conv_unit

        !ALLOCATE
        CALL rams%RR%init_rr_pol(sys%natom, rams%RR%framecount_rtp)

        conv_unit = damping_constant*joule_unit/ev_unit !! J
        damping_factor = conv_unit/action_unit*rams%RR%dt_rtp*fs2s !! s-1

        DO k = 1, 2
            DO i = 1, sys%natom
                DO j = 1, 3
                    DO j_pol = 1, 3
                        DO l = 2, rams%RR%framecount_rtp + 1
                            rams%RR%pol_rtp(1, j_pol)%atom(i)%displacement(k)%XYZ(j)%frame(l-1) = REAL( (rams%RR%static_dip_x_rtp(j_pol)%atom(i)%displacement(k)%XYZ(j)%frame(l)- rams%RR%static_dip_x_rtp(j_pol)%atom(i)%displacement(k)%XYZ(j)%frame(1))  * (EXP(-1.0_dp*damping_factor*(l - 1)))/0.001_dp, kind=dp)
                            rams%RR%pol_rtp(2, j_pol)%atom(i)%displacement(k)%XYZ(j)%frame(l-1) = REAL( (rams%RR%static_dip_y_rtp(j_pol)%atom(i)%displacement(k)%XYZ(j)%frame(l)- rams%RR%static_dip_y_rtp(j_pol)%atom(i)%displacement(k)%XYZ(j)%frame(1))  * (EXP(-1.0_dp*damping_factor*(l-1)))/0.001_dp , kind=dp)
                            rams%RR%pol_rtp(3, j_pol)%atom(i)%displacement(k)%XYZ(j)%frame(l-1) = REAL( (rams%RR%static_dip_z_rtp(j_pol)%atom(i)%displacement(k)%XYZ(j)%frame(l)- rams%RR%static_dip_z_rtp(j_pol)%atom(i)%displacement(k)%XYZ(j)%frame(1))  * (EXP(-1.0_dp*damping_factor*(l-1)))/0.001_dp , kind=dp)
                        END DO
                    END DO
                END DO
            END DO
        END DO

        ! DO i = 1, framecount_rtp
        !    pol_rtp(:,:,:,:,:,i) = pol_rtp(:,:,:,:,:,i)*((COS(i/(framecount_rtp - 1.0_dp)/2.0_dp*3.14_dp))**2) !!Hann Window function
        !END DO
    END SUBROUTINE finite_diff_static_resraman
END MODULE fin_diff
