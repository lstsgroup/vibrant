MODULE fin_diff
        
IMPLICIT NONE
PUBLIC :: central_diff, forward_diff, finite_diff_static,finite_diff_static_resraman

CONTAINS
SUBROUTINE central_diff(dt,natom,framecount,coord_v,v,read_function,mol_num,system)

CHARACTER(LEN=40),INTENT(INOUT)                          :: read_function,system
INTEGER,INTENT(INOUT)                                    :: natom,framecount,mol_num
REAL(KIND=8),INTENT(INOUT)                               :: dt
REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE,INTENT(INOUT)  :: coord_v
REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE,INTENT(OUT)    :: v

INTEGER                                                  :: stat,i,j,k,m

ALLOCATE(v(framecount,natom,3))
!ALLOCATE(v(framecount,1:44,3)) !change for fragments

DO j=1,framecount-2
    DO i=1,natom
        DO k=1,3
            v(j,i,k)=(coord_v(j+2,i,k)-coord_v(j,i,k))/REAL(2.0d0*dt,KIND=8)
        ENDDO
    ENDDO
ENDDO

END SUBROUTINE central_diff

!**************************************************************************************************************!
!**************************************************************************************************************!
SUBROUTINE forward_diff(mol_num,framecount,alpha,dip_free,dip_x,system,read_function)

CHARACTER(LEN=40),INTENT(IN)                             :: read_function
CHARACTER(LEN=40),INTENT(INOUT)                          :: system
INTEGER,INTENT(INOUT)                                    :: mol_num
INTEGER,INTENT(INOUT)                                    :: framecount
REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE,INTENT(INOUT)  :: dip_free,dip_x
REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE,INTENT(OUT)    :: alpha

INTEGER                                                  :: stat,i,j,k,m

ALLOCATE(alpha(framecount,mol_num,3))
!ALLOCATE(alpha(framecount,8:37,3))

IF (read_function.NE.'MD-RR') THEN
    DO j=1,framecount
        DO i=1,mol_num  !!! change to mol_num later
            DO k=1,3
                alpha(j,i,k)=REAL((dip_x(j,i,k)-dip_free(j,i,k))/0.005d0,KIND=8)
                !alpha_x(j,i,k)=(dip_x(j,i,k)-dip_free(j,i,k))
            ENDDO
        ENDDO
    ENDDO

ELSEIF (read_function=='MD-RR') THEN
    DO j=1,framecount
        DO i=2,mol_num 
            DO k=1,3
                alpha(j,i-1,k)=REAL((dip_x(j,i,k)-dip_x(j,1,k))*(EXP(-7.0d0*1.0d0*REAL(i/(mol_num-1.0d0),KIND=8))**2.0d0)&
                               /0.005d0,KIND=8)
            ENDDO
        ENDDO
    ENDDO

    DO j=1,framecount
        DO i=2,mol_num
            DO k=1,3
                alpha(j,i-1,k)= alpha(j,i-1,k)*0.5d0*(1+COS(2.0d0*3.14d0*i/(2.0d0*(mol_num-1))))
            ENDDO
        ENDDO
    ENDDO
ENDIF

END SUBROUTINE forward_diff

!**************************************************************************************************************!
!**************************************************************************************************************!
SUBROUTINE finite_diff_static(natom,nmodes,pol,pol_dq,disp,mass_atom,dx,bohr2ang,static_dipole_free,&
          static_dipole_x,static_dipole_y,static_dipole_z,type_dipole)

INTEGER,INTENT(INOUT)                                        :: natom,nmodes
CHARACTER(LEN=40),INTENT(IN)                                 :: type_dipole
REAL(KIND=8),DIMENSION(:,:,:,:,:),ALLOCATABLE,INTENT(INOUT)  :: pol
REAL(KIND=8),DIMENSION(:,:,:,:),ALLOCATABLE,INTENT(INOUT)    :: static_dipole_free,static_dipole_x
REAL(KIND=8),DIMENSION(:,:,:,:),ALLOCATABLE,INTENT(INOUT)    :: static_dipole_y,static_dipole_z
REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE,INTENT(OUT)        :: pol_dq
REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE,INTENT(INOUT)      :: disp
REAL(KIND=8),DIMENSION(:),ALLOCATABLE,INTENT(INOUT)          :: mass_atom
REAL(KIND=8),INTENT(IN)                                      :: dx,bohr2ang

INTEGER                                                      :: stat,i,j,k,m
REAL(KIND=8)                                                 :: factor
REAL(KIND=8),DIMENSION(:),ALLOCATABLE                        :: mass_inv_sqrt
REAL(KIND=8),DIMENSION(:,:,:,:),ALLOCATABLE                  :: pol_dxyz

ALLOCATE(pol_dxyz(natom,3,3,3))
ALLOCATE(pol_dq(nmodes,3,3))
ALLOCATE(mass_inv_sqrt(natom))

print*,static_dipole_free(1,1,1,1),'free',static_dipole_x(1,1,1,1),'x',"polarizabilities1"
pol_dq=0.0d0

factor=REAL(1.0d0/(2.0d0*dx),KIND=8)

mass_inv_sqrt(:)=REAL(1.0d0/SQRT(mass_atom(:)),KIND=8)

IF (type_dipole=='2') THEN
    DO j=1,natom
        DO i=1,3
            DO k=1,2
                DO m=1,3
                    pol(j,i,k,1,m)=REAL((static_dipole_x(j,i,k,m)-static_dipole_free(j,i,k,m))/0.005d0,KIND=8)
                    pol(j,i,k,2,m)=REAL((static_dipole_y(j,i,k,m)-static_dipole_free(j,i,k,m))/0.005d0,KIND=8)
                    pol(j,i,k,3,m)=REAL((static_dipole_z(j,i,k,m)-static_dipole_free(j,i,k,m))/0.005d0,KIND=8)
                ENDDO
            ENDDO
        ENDDO
    ENDDO
ENDIF

print*,static_dipole_free(1,1,1,1),'free',static_dipole_x(1,1,1,1),'x',"polarizabilities1"
print*,pol(1,1,1,1,1),"polarizabilities"

DO j=1,natom
    DO i=1,3
        pol_dxyz(j,i,:,:)=REAL((pol(j,i,1,:,:)-pol(j,i,2,:,:))*factor,KIND=8)
    ENDDO
ENDDO

DO j=1,nmodes
    DO k=1,natom
        pol_dq(j,:,:)=pol_dq(j,:,:)+(pol_dxyz(k,1,:,:)*disp(k,j,1)*mass_inv_sqrt(k))&
                      +(pol_dxyz(k,2,:,:)*disp(k,j,2)*mass_inv_sqrt(k))+(pol_dxyz(k,3,:,:)*&
                      disp(k,j,3)*mass_inv_sqrt(k))
    ENDDO
ENDDO

DEALLOCATE(pol_dxyz,mass_inv_sqrt)

END SUBROUTINE finite_diff_static

!**************************************************************************************************************!
!**************************************************************************************************************!
SUBROUTINE finite_diff_static_resraman(natom,pol_rtp,static_dipole_x_rtp,static_dipole_y_rtp,static_dipole_z_rtp,&
         framecount_rtp,speed_light,fs2s,damping_constant,joule_unit,ev_unit,action_unit,dt_rtp)

INTEGER,INTENT(IN)                                           :: natom,framecount_rtp
REAL(KIND=8),INTENT(IN)                                      :: speed_light,fs2s,damping_constant
REAL(KIND=8),INTENT(IN)                                      :: joule_unit,ev_unit,action_unit,dt_rtp
REAL(KIND=8),DIMENSION(:,:,:,:,:,:),ALLOCATABLE,INTENT(OUT)  :: pol_rtp
REAL(KIND=8),DIMENSION(:,:,:,:,:),ALLOCATABLE,INTENT(INOUT)  :: static_dipole_x_rtp,static_dipole_y_rtp
REAL(KIND=8),DIMENSION(:,:,:,:,:),ALLOCATABLE,INTENT(INOUT)  :: static_dipole_z_rtp

INTEGER                                                      :: stat,i,j,k,m,l
REAL(KIND=8)                                                 :: damping_factor,conv_unit

ALLOCATE(pol_rtp(natom,3,2,3,3,framecount_rtp))

conv_unit=damping_constant*joule_unit/ev_unit !! J
damping_factor=conv_unit/action_unit*dt_rtp*fs2s !! s-1 

DO j=1,natom
  DO i=1,3
   DO k=1,2
    DO m=1,3
     DO l=2,framecount_rtp+1
       pol_rtp(j,i,k,1,m,l-1)=REAL((static_dipole_x_rtp(j,i,k,m,l)-static_dipole_x_rtp(j,i,k,m,1))*&
       (EXP(-1.0d0*damping_factor*(l-1)))/0.001d0,KIND=8)
       pol_rtp(j,i,k,2,m,l-1)=REAL((static_dipole_y_rtp(j,i,k,m,l)-static_dipole_y_rtp(j,i,k,m,1))*&
       (EXP(-1.0d0*damping_factor*(l-1)))/0.001d0,KIND=8)
       pol_rtp(j,i,k,3,m,l-1)=REAL((static_dipole_z_rtp(j,i,k,m,l)-static_dipole_z_rtp(j,i,k,m,1))*&
       (EXP(-1.0d0*damping_factor*(l-1)))/0.001d0,KIND=8)
     ENDDO
    ENDDO
   ENDDO
 ENDDO
ENDDO

END SUBROUTINE finite_diff_static_resraman
END MODULE fin_diff

