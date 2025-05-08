MODULE fin_diff
    USE kinds, ONLY: dp
    IMPLICIT NONE
    PUBLIC :: central_diff, forward_diff, finite_diff_static, finite_diff_static_resraman

CONTAINS
    SUBROUTINE central_diff(dt, natom, framecount, coord_v, v, read_function, mol_num, system)

        CHARACTER(LEN=40), INTENT(INOUT)                          :: read_function, system
        INTEGER, INTENT(INOUT)                                    :: natom, framecount, mol_num
        REAL(kind=dp), INTENT(INOUT)                               :: dt
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE, INTENT(INOUT)  :: coord_v
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE, INTENT(OUT)    :: v

        INTEGER                                                  :: stat, i, j, k, m

        ALLOCATE (v(framecount, natom, 3))
!ALLOCATE(v(framecount,1:44,3)) !change for fragments

        DO j = 1, framecount - 2
            DO i = 1, natom
                DO k = 1, 3
                    v(j, i, k) = (coord_v(j + 2, i, k) - coord_v(j, i, k))/REAL(2.0_dp*dt, kind=dp)
                END DO
            END DO
        END DO

    END SUBROUTINE central_diff

!**************************************************************************************************************!
!**************************************************************************************************************!
    SUBROUTINE forward_diff(mol_num, framecount, alpha, dip_free, dip_x, system, read_function)

        CHARACTER(LEN=40), INTENT(IN)                             :: read_function
        CHARACTER(LEN=40), INTENT(INOUT)                          :: system
        INTEGER, INTENT(INOUT)                                    :: mol_num
        INTEGER, INTENT(INOUT)                                    :: framecount
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE, INTENT(INOUT)  :: dip_free, dip_x
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE, INTENT(OUT)    :: alpha

        INTEGER                                                  :: stat, i, j, k, m

        ALLOCATE (alpha(framecount, mol_num, 3))
!ALLOCATE(alpha(framecount,8:37,3))

        IF (read_function.NE.'MD-RR') THEN
            DO j = 1, framecount
                DO i = 1, mol_num  !!! change to mol_num later
                    DO k = 1, 3
                        alpha(j, i, k) = REAL((dip_x(j, i, k) - dip_free(j, i, k))/0.005_dp, kind=dp)
                        !alpha_x(j,i,k)=(dip_x(j,i,k)-dip_free(j,i,k))
                    END DO
                END DO
            END DO

        ELSEIF (read_function=='MD-RR') THEN
            DO j = 1, framecount
                DO i = 2, mol_num
                    DO k = 1, 3
                        alpha(j, i - 1, k) = REAL((dip_x(j, i, k) - dip_x(j, 1, k)) &
                                                  *(EXP(-7.0_dp*1.0_dp*REAL(i/(mol_num - 1.0_dp), kind=dp))**2.0_dp)/0.005_dp, &
                                                  kind=dp)
                    END DO
                END DO
            END DO

            DO j = 1, framecount
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
    SUBROUTINE finite_diff_static(natom, nmodes, pol, pol_dq, disp, atom_mass_inv_sqrt, dx, bohr2ang, static_dip_free, &
                                  static_dip_x, static_dip_y, static_dip_z, type_dipole, read_function, dip_dq)

        INTEGER, INTENT(INOUT)                                         :: natom, nmodes
        CHARACTER(LEN=40), INTENT(IN)                                  :: read_function, type_dipole
        REAL(kind=dp), DIMENSION(:, :, :, :, :), ALLOCATABLE, INTENT(INOUT)  :: pol
        REAL(kind=dp), DIMENSION(:, :, :, :), ALLOCATABLE, INTENT(INOUT)    :: static_dip_free, static_dip_x
        REAL(kind=dp), DIMENSION(:, :, :, :), ALLOCATABLE, INTENT(INOUT)    :: static_dip_y, static_dip_z
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE, INTENT(OUT)        :: pol_dq
        REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE, INTENT(OUT)        :: dip_dq
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE, INTENT(IN)         :: disp
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(IN)           :: atom_mass_inv_sqrt
        REAL(kind=dp), INTENT(IN)                                      :: dx, bohr2ang

        INTEGER                                                      :: stat, i, j, k, m
        REAL(kind=dp)                                                 :: factor
        REAL(kind=dp), DIMENSION(:, :, :, :), ALLOCATABLE                  :: pol_dxyz
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE                  :: dip_dxyz

        factor = REAL(1.0_dp/(2.0_dp*dx), kind=dp)
 
        IF (read_function=='R') THEN 
            ALLOCATE (pol_dxyz(natom, 3, 3, 3))
            ALLOCATE (pol_dq(nmodes, 3, 3))
            
            pol_dq = 0.0_dp
            
            IF (type_dipole=='2') THEN
                pol(:, :, :, 1, :) = REAL((static_dip_x(:, :, :, :) - static_dip_free(:, :, :, :)) / (5.338D-5 * 1.313D-26), kind=dp)
                pol(:, :, :, 2, :) = REAL((static_dip_y(:, :, :, :) - static_dip_free(:, :, :, :)) / (5.338D-5 * 1.313D-26), kind=dp)
                pol(:, :, :, 3, :) = REAL((static_dip_z(:, :, :, :) - static_dip_free(:, :, :, :)) / (5.338D-5 * 1.313D-26), kind=dp)
                
                DEALLOCATE(static_dip_free,static_dip_x,static_dip_y,static_dip_z)
            END IF

            pol_dxyz(:, :, :, :) = REAL((pol(:, :, 1, :, :) - pol(:, :, 2, :, :))*factor, kind=dp)
        
        
        ELSEIF (read_function=='IR') THEN
        
            ALLOCATE (dip_dxyz(natom, 3, 3))
            ALLOCATE (dip_dq(nmodes, 3))
            
            dip_dq = 0.0_dp
       
            dip_dxyz(:, :, :) = REAL((static_dip_free(:, :, 1, :) - static_dip_free(:, :, 2, :))*factor, kind=dp)
        
            DEALLOCATE(static_dip_free)
        END IF
        
        DO j = 1, nmodes
            DO k = 1, natom
                IF (read_function=='R') THEN
                    pol_dq(j, :, :) = pol_dq(j, :, :) + (pol_dxyz(k, 1, :, :)*disp(j, k, 1)*atom_mass_inv_sqrt(k)) &
                                      + (pol_dxyz(k, 2, :, :)*disp(j, k, 2)*atom_mass_inv_sqrt(k)) + (pol_dxyz(k, 3, :, :)* &
                                      disp(j, k, 3)*atom_mass_inv_sqrt(k))
                ELSEIF (read_function=='IR') THEN
                    dip_dq(j, :) = dip_dq(j, :) + (dip_dxyz(k, 1, :)*disp(j, k, 1)*atom_mass_inv_sqrt(k)) &
                                      + (dip_dxyz(k, 2, :)*disp(j, k, 2)*atom_mass_inv_sqrt(k)) + (dip_dxyz(k, 3, :)* &
                                      disp(j, k, 3)*atom_mass_inv_sqrt(k))
                END IF
            END DO
        END DO
        
        IF (read_function=='R') THEN 
             DEALLOCATE (pol_dxyz,pol)
        ELSEIF (read_function=='IR') THEN 
             DEALLOCATE (dip_dxyz)
        ENDIF
    
    END SUBROUTINE finite_diff_static

!**************************************************************************************************************!
!**************************************************************************************************************!
    SUBROUTINE finite_diff_static_resraman(natom, pol_rtp, static_dipole_x_rtp, static_dipole_y_rtp, static_dipole_z_rtp, &
                                           framecount_rtp, speed_light, fs2s, damping_constant, joule_unit, &
                                           ev_unit, action_unit, dt_rtp)

        INTEGER, INTENT(IN)                                           :: natom, framecount_rtp
        REAL(kind=dp), INTENT(IN)                                      :: speed_light, fs2s, damping_constant
        REAL(kind=dp), INTENT(IN)                                      :: joule_unit, ev_unit, action_unit, dt_rtp
        REAL(kind=dp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE, INTENT(OUT)  :: pol_rtp
        REAL(kind=dp), DIMENSION(:, :, :, :, :), ALLOCATABLE, INTENT(INOUT)  :: static_dipole_x_rtp, static_dipole_y_rtp
        REAL(kind=dp), DIMENSION(:, :, :, :, :), ALLOCATABLE, INTENT(INOUT)  :: static_dipole_z_rtp

        INTEGER                                                      :: stat, i, j, k, m, l
        REAL(kind=dp)                                                 :: damping_factor, conv_unit

        ALLOCATE (pol_rtp(natom, 3, 2, 3, 3, framecount_rtp))

        conv_unit = damping_constant*joule_unit/ev_unit !! J
        damping_factor = conv_unit/action_unit*dt_rtp*fs2s !! s-1
                            
        DO l = 2, framecount_rtp + 1
            pol_rtp(:, :, :, 1, :, l - 1) = static_dipole_x_rtp(:, :, :, :, l) - static_dipole_x_rtp(:, :, :, :, 1)
            pol_rtp(:, :, :, 2, :, l - 1) = static_dipole_y_rtp(:, :, :, :, l) - static_dipole_y_rtp(:, :, :, :, 1)
            pol_rtp(:, :, :, 3, :, l - 1) = static_dipole_z_rtp(:, :, :, :, l) - static_dipole_z_rtp(:, :, :, :, 1)
        ENDDO

      DO l=1,framecount_rtp
        pol_rtp(:,:,:,:,:,l)=pol_rtp(:,:,:,:,:,l)*(EXP(-1.0d0*damping_factor*l))
      ENDDO


       ! DO i = 1, framecount_rtp
        !    pol_rtp(:,:,:,:,:,i) = pol_rtp(:,:,:,:,:,i)*((COS(i/(framecount_rtp - 1.0_dp)/2.0_dp*3.14_dp))**2) !!Hann Window function
        !END DO
    END SUBROUTINE finite_diff_static_resraman
END MODULE fin_diff

