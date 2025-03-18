MODULE vel_cor

USE dipole_calc,          ONLY: center_mass
USE kinds,              ONLY: dp
        
IMPLICIT NONE
PUBLIC :: cvv, cvv_iso, cvv_aniso, cvv_only_x, cvv_resraman

CONTAINS
SUBROUTINE cvv(natom,framecount,t_cor,coord_v,z,type_input,dt,input_mass,mass_atom,mass_tot,pi,&
                mol_num,read_function,system,frag_type)

CHARACTER(LEN=40),INTENT(INOUT)                          :: system,type_input,input_mass,read_function,frag_type
INTEGER,INTENT(INOUT)                                    :: natom,framecount,mol_num,t_cor
REAL(kind=dp),INTENT(IN)                                  :: dt,pi
REAL(kind=dp),INTENT(INOUT)                               :: mass_tot
REAL(kind=dp),DIMENSION(:,:,:),ALLOCATABLE,INTENT(INOUT)  :: coord_v
REAL(kind=dp),DIMENSION(:),ALLOCATABLE,INTENT(OUT)        :: z
REAL(kind=dp),DIMENSION(:),ALLOCATABLE,INTENT(IN)         :: mass_atom

CHARACTER(LEN=40)                                        :: chara
INTEGER                                                  :: stat,i,j,k,m,t0,t1,l
INTEGER,DIMENSION(:),ALLOCATABLE                         :: norm
REAL(kind=dp),DIMENSION(:,:),ALLOCATABLE                  :: coord

ALLOCATE(z(0:2*t_cor-1),norm(0:2*t_cor-1))

norm=0
z=0.0_dp
k=0
j=0

DO t0=1,framecount-2
    t1=MIN(framecount-2,t0+t_cor)
    DO j=1,natom
        IF (frag_type=='2') THEN 
            k=j+8
        ELSEIF (frag_type=='3') THEN 
            k=j+20
        ELSE
            k=j
        ENDIF
        IF (input_mass=='y') THEN
            z(0:t1-t0)=z(0:t1-t0)+(coord_v(t0,k,1)*coord_v(t0:t1,j,1)+coord_v(t0,j,2)*coord_v(t0:t1,j,2)+&
            coord_v(t0,j,3)*coord_v(t0:t1,j,3))*mass_atom(j)
        ELSE
              !  print*,k,"k",j,"j"
            DO m=1,3
                z(0:t1-t0)=z(0:t1-t0)+coord_v(t0,k,m)*coord_v(t0:t1,k,m)
            ENDDO
        ENDIF
    ENDDO
    norm(0:t1-t0)=norm(0:t1-t0)+1
ENDDO

!IF (type_input=='2') THEN        
 !  z(:)=z(:)*4.785992e12
!ELSEIF (type_input=='1') THEN 
!   z(:)=z(:)*1e10
!ENDIF

z(:)=z(:)/norm(:) !!Normalization

z(t_cor)=0.0_dp
DO i=1,t_cor-1                                                                                                                    
 z(t_cor+i)=z(t_cor-i) !!Data mirroring
ENDDO

DO i=0,2*t_cor-1
 z(i)=z(i)*((COS(i/(t_cor-1.0_dp)/ 2.0_dp*3.14_dp))**2) !!Hann Window function
ENDDO

OPEN(UNIT=61,FILE='result_cvv.txt',STATUS='unknown',IOSTAT=stat) !!Write output
DO i=0,2*t_cor-1
 WRITE(61,*) z(i)
ENDDO
CLOSE(61)

DEALLOCATE(norm)

END SUBROUTINE cvv

!***********************************************************************************
!***********************************************************************************

SUBROUTINE cvv_iso(mol_num,framecount,t_cor,z_iso,alpha_diff_x,alpha_diff_y,alpha_diff_z,dt,&
           pi,frag_type)

CHARACTER(LEN=40),INTENT(INOUT)                          :: frag_type
INTEGER,INTENT(INOUT)                                    :: framecount,mol_num,t_cor
REAL(kind=dp),INTENT(IN)                                  :: dt,pi
REAL(kind=dp),DIMENSION(:),ALLOCATABLE,INTENT(OUT)        :: z_iso
REAL(kind=dp),DIMENSION(:,:,:),ALLOCATABLE,INTENT(INOUT)  :: alpha_diff_x,alpha_diff_y,alpha_diff_z

INTEGER                                                  :: stat,i,j,k,m,t0,t1
INTEGER,DIMENSION(:),ALLOCATABLE                         :: norm

ALLOCATE(z_iso(0:2*t_cor-1),norm(0:2*t_cor-1))

 norm=0
 z_iso=0.0_dp
 DO t0=1,framecount-2
  t1=MIN(framecount-2,t0+t_cor)
  DO j=1,mol_num
      IF (frag_type=='2') THEN 
          k=j+8
      ELSEIF (frag_type=='3') THEN 
          k=j+20
      ELSE 
          k=j
      ENDIF
     z_iso(0:t1-t0)=z_iso(0:t1-t0)+(alpha_diff_x(t0,k,1)+alpha_diff_y(t0,k,2)+alpha_diff_z(t0,k,3))*&
            (alpha_diff_x(t0:t1,k,1)+alpha_diff_y(t0:t1,k,2)+alpha_diff_z(t0:t1,k,3))
  ENDDO
     norm(0:t1-t0)=norm(0:t1-t0)+1
 ENDDO

print*,norm(0:3),'iso norm'
print*,z_iso(0:3),'iso z'
z_iso(:)=z_iso(:)/norm(:)
z_iso(:)=z_iso(:)/9._dp
z_iso(:)=z_iso(:)/(2.0_dp*pi)
!z_iso(:)=z_iso(:)/mol_num

DO i=0,t_cor-1
       z_iso(i)=z_iso(i)*((COS(i/(t_cor-1.0_dp)/ 2.0_dp*3.14_dp))**2)
       !z_iso(i)=z_iso(i)*0.5_dp*(1+COS(2.0_dp*3.14_dp*i/(2.0_dp*(t_cor-1))))
ENDDO

z_iso(t_cor)=0.0_dp

DO i=1,t_cor-1                                                                                                                    
        z_iso(t_cor+i)=z_iso(t_cor-i)
ENDDO

OPEN(UNIT=61,FILE='result_cvv_iso.txt',STATUS='unknown',IOSTAT=stat)

DO i=0,2*t_cor-1
        WRITE(61,*) z_iso(i)
ENDDO
CLOSE(61)
DEALLOCATE(norm)

END SUBROUTINE cvv_iso

!***********************************************************************************
!***********************************************************************************

SUBROUTINE cvv_aniso(mol_num,natom,framecount,t_cor,z_aniso,alpha_diff_x,alpha_diff_y,&
           alpha_diff_z,dt,pi,frag_type)

CHARACTER(LEN=40),INTENT(INOUT)                          :: frag_type
INTEGER,INTENT(INOUT)                                    :: natom,framecount,mol_num
INTEGER,INTENT(INOUT)                                    :: t_cor
REAL(kind=dp),INTENT(IN)                                  :: dt,pi
REAL(kind=dp),DIMENSION(:),ALLOCATABLE,INTENT(OUT)        :: z_aniso
REAL(kind=dp),DIMENSION(:,:,:),ALLOCATABLE,INTENT(INOUT)  :: alpha_diff_x,alpha_diff_y,alpha_diff_z

CHARACTER(LEN=40)                                        :: chara
CHARACTER(LEN=2),DIMENSION(:),ALLOCATABLE                :: element
INTEGER                                                  :: stat,i,j,k,m,t0,t1
REAL(kind=dp),DIMENSION(:,:),ALLOCATABLE                  :: coord
INTEGER,DIMENSION(:),ALLOCATABLE                         :: norm

ALLOCATE(z_aniso(0:2*t_cor-1),norm(0:2*t_cor-1))

norm=0
z_aniso=0.0_dp
 
DO t0=1,framecount-2
    t1=MIN(framecount-2,t0+t_cor)
    DO j=1,mol_num
        IF (frag_type=='2') THEN 
            k=j+8
        ELSEIF (frag_type=='3') THEN 
            k=j+20
        ELSE 
            k=j
        ENDIF
        z_aniso(0:t1-t0)=z_aniso(0:t1-t0)+(alpha_diff_x(t0,k,1)-alpha_diff_y(t0,k,2))*(alpha_diff_x(t0:t1,k,1)&
             -alpha_diff_y(t0:t1,k,2))/2.0_dp
        z_aniso(0:t1-t0)=z_aniso(0:t1-t0)+(alpha_diff_y(t0,k,2)-alpha_diff_z(t0,k,3))*(alpha_diff_y(t0:t1,k,2)&
             -alpha_diff_z(t0:t1,k,3))/2.0_dp
        z_aniso(0:t1-t0)=z_aniso(0:t1-t0)+(alpha_diff_z(t0,k,3)-alpha_diff_x(t0,k,1))*(alpha_diff_z(t0:t1,k,3)&
             -alpha_diff_x(t0:t1,k,1))/2.0_dp
        z_aniso(0:t1-t0)=z_aniso(0:t1-t0)+(alpha_diff_x(t0,k,2)*0.50_dp+alpha_diff_y(t0,k,1)*0.50_dp)&
             *(alpha_diff_x(t0:t1,k,2)*0.50_dp+alpha_diff_y(t0:t1,k,1)*0.50_dp)*3.0_dp
        z_aniso(0:t1-t0)=z_aniso(0:t1-t0)+(alpha_diff_y(t0,k,3)*0.50_dp+alpha_diff_z(t0,k,2)*0.50_dp)&
             *(alpha_diff_y(t0:t1,k,3)*0.50_dp+alpha_diff_z(t0:t1,k,2)*0.50_dp)*3.0_dp
        z_aniso(0:t1-t0)=z_aniso(0:t1-t0)+(alpha_diff_z(t0,k,1)*0.50_dp+alpha_diff_x(t0,k,3)*0.50_dp)&
             *(alpha_diff_z(t0:t1,k,1)*0.50_dp+alpha_diff_x(t0:t1,k,3)*0.50_dp)*3.0_dp
     ENDDO
     norm(0:t1-t0)=norm(0:t1-t0)+1
 ENDDO


print*,norm(0:3),'aniso norm'
print*,z_aniso(0:3),'aniso z'
z_aniso(:)=z_aniso(:)/norm(:)
z_aniso(:)=z_aniso(:)/(2.0_dp*pi)
!z_aniso(:)=REAL(z_aniso(:)/mol_num,kind=dp)

DO i=0,t_cor-1
       z_aniso(i)=z_aniso(i)*((COS(i/(t_cor-1.0_dp)/ 2.0_dp*3.14_dp))**2)
ENDDO

z_aniso(t_cor)=0.0_dp

DO i=1,t_cor-1                                                                                                                    
        z_aniso(t_cor+i)=z_aniso(t_cor-i)
ENDDO

OPEN(UNIT=61,FILE='result_cvv_aniso.txt',STATUS='unknown',IOSTAT=stat)

DO i=0,2*t_cor-1
        WRITE(61,*) z_aniso(i)
ENDDO
CLOSE(61) 

DEALLOCATE(norm)

END SUBROUTINE cvv_aniso

!***********************************************************************************
!***********************************************************************************

SUBROUTINE cvv_resraman(framecount,natom,t_cor,dt,pi,alpha_resraman_x_diff_re,alpha_resraman_y_diff_re,&
           alpha_resraman_z_diff_re,alpha_resraman_x_diff_im,alpha_resraman_y_diff_im,alpha_resraman_z_diff_im,&
           z_iso_resraman,z_aniso_resraman)

INTEGER,INTENT(INOUT)                                    :: framecount,natom,t_cor
REAL(kind=dp),INTENT(IN)                                  :: dt,pi
REAL(kind=dp),DIMENSION(:,:,:),ALLOCATABLE,INTENT(INOUT)  :: alpha_resraman_x_diff_re,alpha_resraman_y_diff_re
REAL(kind=dp),DIMENSION(:,:,:),ALLOCATABLE,INTENT(INOUT)  :: alpha_resraman_z_diff_re
REAL(kind=dp),DIMENSION(:,:,:),ALLOCATABLE,INTENT(INOUT)  :: alpha_resraman_x_diff_im,alpha_resraman_y_diff_im
REAL(kind=dp),DIMENSION(:,:,:),ALLOCATABLE,INTENT(INOUT)  :: alpha_resraman_z_diff_im
COMPLEX(kind=dp),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT)   :: z_iso_resraman,z_aniso_resraman

INTEGER                                                  :: stat,i,j,k,m,t0,t1
INTEGER,DIMENSION(:,:),ALLOCATABLE                       :: norm_iso,norm_aniso
COMPLEX(kind=dp)                                          :: im_unit

!!!ISOTROPIC!!!

ALLOCATE(z_iso_resraman(0:2*t_cor,natom-1),norm_iso(0:2*t_cor,natom-1))
ALLOCATE(z_aniso_resraman(0:2*t_cor,natom-1),norm_aniso(0:2*t_cor,natom-1))

im_unit = (0.0_dp, 1.0_dp)
z_iso_resraman = (0.0_dp, 0.0_dp)
z_aniso_resraman = (0.0_dp, 0.0_dp)

framecount=framecount-2

norm_iso=0.0_dp
DO t0=2,framecount
 t1=MIN(framecount,t0+t_cor)
 DO k=1,natom-1
  !!RE*RE
    z_iso_resraman(0:t1-t0,k)=z_iso_resraman(0:t1-t0,k)+(alpha_resraman_x_diff_re(t0,k,1)+alpha_resraman_y_diff_re(t0,k,2)+&
            alpha_resraman_z_diff_re(t0,k,3))*(alpha_resraman_x_diff_re(t0:t1,k,1)+alpha_resraman_y_diff_re(t0:t1,k,2)+&
            alpha_resraman_z_diff_re(t0:t1,k,3))
  !!IM*IM
    z_iso_resraman(0:t1-t0,k)=z_iso_resraman(0:t1-t0,k)+(alpha_resraman_x_diff_im(t0,k,1)+alpha_resraman_y_diff_im(t0,k,2)+&
            alpha_resraman_z_diff_im(t0,k,3))*(alpha_resraman_x_diff_im(t0:t1,k,1)+alpha_resraman_y_diff_im(t0:t1,k,2)+&
            alpha_resraman_z_diff_im(t0:t1,k,3))
  !!RE*IM
    z_iso_resraman(0:t1-t0,k)=z_iso_resraman(0:t1-t0,k)+((alpha_resraman_x_diff_re(t0,k,1)+alpha_resraman_y_diff_re(t0,k,2)+&
            alpha_resraman_z_diff_re(t0,k,3))*(alpha_resraman_x_diff_im(t0:t1,k,1)+alpha_resraman_y_diff_im(t0:t1,k,2)+&
            alpha_resraman_z_diff_im(t0:t1,k,3)))*im_unit
  !!IM*RE (SUBTRACT)
    z_iso_resraman(0:t1-t0,k)=z_iso_resraman(0:t1-t0,k)-((alpha_resraman_x_diff_im(t0,k,1)+alpha_resraman_y_diff_im(t0,k,2)+&
            alpha_resraman_z_diff_im(t0,k,3))*(alpha_resraman_x_diff_re(t0:t1,k,1)+alpha_resraman_y_diff_re(t0:t1,k,2)+&
            alpha_resraman_z_diff_re(t0:t1,k,3)))*im_unit
   
    norm_iso(0:t1-t0,k)=norm_iso(0:t1-t0,k)+1.0_dp 
 ENDDO
ENDDO

z_iso_resraman(:,:)=z_iso_resraman(:,:)/norm_iso(:,:)
z_iso_resraman(:,:)=z_iso_resraman(:,:)/9._dp
z_iso_resraman(:,:)=z_iso_resraman(:,:)/(2.0_dp*pi)

DO i=0,t_cor-1
   DO j=1,natom-1
       z_iso_resraman(i,j)=z_iso_resraman(i,j)*((COS(i/(t_cor-1.0_dp)/ 2.0_dp*3.14_dp))**2)
   ENDDO    
ENDDO

DO i=1,natom-1
      z_iso_resraman(t_cor,i)=0.0_dp
ENDDO

DO i=1,t_cor-1
  DO j=1,natom-1
      z_iso_resraman(t_cor+i,j)=z_iso_resraman(t_cor-i,j)
  ENDDO
ENDDO

OPEN(UNIT=61,FILE='result_cvv_iso_resraman.txt',STATUS='unknown',IOSTAT=stat)

DO i=0,2*t_cor-1
  DO j=1,natom-1
        WRITE(61,*) z_iso_resraman(i,j)
  ENDDO
ENDDO
CLOSE(61)
DEALLOCATE(norm_iso)

!!!ANISOTROPIC!!!

norm_aniso=0.0_dp
DO t0=2,framecount
 t1=MIN(framecount,t0+t_cor)
  DO k=1,natom-1
  !!RE*RE
    z_aniso_resraman(0:t1-t0,k)=z_aniso_resraman(0:t1-t0,k)+(alpha_resraman_x_diff_re(t0,k,1)-alpha_resraman_y_diff_re(t0,k,2))&
            *(alpha_resraman_x_diff_re(t0:t1,k,1)-alpha_resraman_y_diff_re(t0:t1,k,2))/2.0_dp
    z_aniso_resraman(0:t1-t0,k)=z_aniso_resraman(0:t1-t0,k)+(alpha_resraman_y_diff_re(t0,k,2)-alpha_resraman_z_diff_re(t0,k,3))&
            *(alpha_resraman_y_diff_re(t0:t1,k,2)-alpha_resraman_z_diff_re(t0:t1,k,3))/2.0_dp
    z_aniso_resraman(0:t1-t0,k)=z_aniso_resraman(0:t1-t0,k)+(alpha_resraman_z_diff_re(t0,k,3)-alpha_resraman_x_diff_re(t0,k,1))&
            *(alpha_resraman_z_diff_re(t0:t1,k,3)-alpha_resraman_x_diff_re(t0:t1,k,1))/2.0_dp
    
   	z_aniso_resraman(0:t1-t0,k)=z_aniso_resraman(0:t1-t0,k)+(alpha_resraman_x_diff_re(t0,k,2)*0.50_dp+&
     	alpha_resraman_y_diff_re(t0,k,1)*0.50_dp)*(alpha_resraman_x_diff_re(t0:t1,k,2)*0.50_dp+&
    	alpha_resraman_y_diff_re(t0:t1,k,1)*0.50_dp)*3.0_dp
    
	z_aniso_resraman(0:t1-t0,k)=z_aniso_resraman(0:t1-t0,k)+(alpha_resraman_y_diff_re(t0,k,3)*0.50_dp+&
    	alpha_resraman_z_diff_re(t0,k,2)*0.50_dp)*(alpha_resraman_y_diff_re(t0:t1,k,3)*0.50_dp+&
    	alpha_resraman_z_diff_re(t0:t1,k,2)*0.50_dp)*3.0_dp
    
	z_aniso_resraman(0:t1-t0,k)=z_aniso_resraman(0:t1-t0,k)+(alpha_resraman_z_diff_re(t0,k,1)*0.50_dp+&
     	alpha_resraman_x_diff_re(t0,k,3)*0.50_dp)*(alpha_resraman_z_diff_re(t0:t1,k,1)*0.50_dp+&
     	alpha_resraman_x_diff_re(t0:t1,k,3)*0.50_dp)*3.0_dp
    
  !!IM*IM
    z_aniso_resraman(0:t1-t0,k)=z_aniso_resraman(0:t1-t0,k)+(alpha_resraman_x_diff_im(t0,k,1)-alpha_resraman_y_diff_im(t0,k,2))&
            *(alpha_resraman_x_diff_im(t0:t1,k,1)-alpha_resraman_y_diff_im(t0:t1,k,2))/2.0_dp
    z_aniso_resraman(0:t1-t0,k)=z_aniso_resraman(0:t1-t0,k)+(alpha_resraman_y_diff_im(t0,k,2)-alpha_resraman_z_diff_im(t0,k,3))&
            *(alpha_resraman_y_diff_im(t0:t1,k,2)-alpha_resraman_z_diff_im(t0:t1,k,3))/2.0_dp
    z_aniso_resraman(0:t1-t0,k)=z_aniso_resraman(0:t1-t0,k)+(alpha_resraman_z_diff_im(t0,k,3)-alpha_resraman_x_diff_im(t0,k,1))&
            *(alpha_resraman_z_diff_im(t0:t1,k,3)-alpha_resraman_x_diff_im(t0:t1,k,1))/2.0_dp
  
    z_aniso_resraman(0:t1-t0,k)=z_aniso_resraman(0:t1-t0,k)+(alpha_resraman_x_diff_im(t0,k,2)*0.50_dp+&
      alpha_resraman_y_diff_im(t0,k,1)*0.50_dp)*(alpha_resraman_x_diff_im(t0:t1,k,2)*0.50_dp+&
      alpha_resraman_y_diff_im(t0:t1,k,1)*0.50_dp)*3.0_dp
  
    z_aniso_resraman(0:t1-t0,k)=z_aniso_resraman(0:t1-t0,k)+(alpha_resraman_y_diff_im(t0,k,3)*0.50_dp+&
    	alpha_resraman_z_diff_im(t0,k,2)*0.50_dp)*(alpha_resraman_y_diff_im(t0:t1,k,3)*0.50_dp+&
		alpha_resraman_z_diff_im(t0:t1,k,2)*0.50_dp)*3.0_dp
  
    z_aniso_resraman(0:t1-t0,k)=z_aniso_resraman(0:t1-t0,k)+(alpha_resraman_z_diff_im(t0,k,1)*0.50_dp+&
	  alpha_resraman_x_diff_im(t0,k,3)*0.50_dp)*(alpha_resraman_z_diff_im(t0:t1,k,1)*0.50_dp+&
	  alpha_resraman_x_diff_im(t0:t1,k,3)*0.50_dp)*3.0_dp
  
  
  !!RE*IM
    z_aniso_resraman(0:t1-t0,k)=z_aniso_resraman(0:t1-t0,k)+((alpha_resraman_x_diff_re(t0,k,1)-alpha_resraman_y_diff_re(t0,k,2))&
            *(alpha_resraman_x_diff_im(t0:t1,k,1)-alpha_resraman_y_diff_im(t0:t1,k,2)))/2.0_dp*im_unit
    z_aniso_resraman(0:t1-t0,k)=z_aniso_resraman(0:t1-t0,k)+((alpha_resraman_y_diff_re(t0,k,2)-alpha_resraman_z_diff_re(t0,k,3))&
            *(alpha_resraman_y_diff_im(t0:t1,k,2)-alpha_resraman_z_diff_im(t0:t1,k,3)))/2.0_dp*im_unit
    z_aniso_resraman(0:t1-t0,k)=z_aniso_resraman(0:t1-t0,k)+((alpha_resraman_z_diff_re(t0,k,3)-alpha_resraman_x_diff_re(t0,k,1))&
            *(alpha_resraman_z_diff_im(t0:t1,k,3)-alpha_resraman_x_diff_im(t0:t1,k,1)))/2.0_dp*im_unit
   
   z_aniso_resraman(0:t1-t0,k)=z_aniso_resraman(0:t1-t0,k)+((alpha_resraman_x_diff_re(t0,k,2)*&
	0.50_dp+alpha_resraman_y_diff_re(t0,k,1)&
            *0.50_dp)*(alpha_resraman_x_diff_im(t0:t1,k,2)*0.50_dp+alpha_resraman_y_diff_im(t0:t1,k,1)*0.50_dp))*3.0_dp*im_unit
    z_aniso_resraman(0:t1-t0,k)=z_aniso_resraman(0:t1-t0,k)+((alpha_resraman_y_diff_re(t0,k,3)*&
	0.50_dp+alpha_resraman_z_diff_re(t0,k,2)&
            *0.50_dp)*(alpha_resraman_y_diff_im(t0:t1,k,3)*0.50_dp+alpha_resraman_z_diff_im(t0:t1,k,2)*0.50_dp))*3.0_dp*im_unit
    z_aniso_resraman(0:t1-t0,k)=z_aniso_resraman(0:t1-t0,k)+((alpha_resraman_z_diff_re(t0,k,1)*&
	0.50_dp+alpha_resraman_x_diff_re(t0,k,3)&
            *0.50_dp)*(alpha_resraman_z_diff_im(t0:t1,k,1)*0.50_dp+alpha_resraman_x_diff_im(t0:t1,k,3)*0.50_dp))*3.0_dp*im_unit
    
  !!IM*RE (SUBTRACT)
    z_aniso_resraman(0:t1-t0,k)=z_aniso_resraman(0:t1-t0,k)-((alpha_resraman_x_diff_im(t0,k,1)-alpha_resraman_y_diff_im(t0,k,2))&
            *(alpha_resraman_x_diff_re(t0:t1,k,1)-alpha_resraman_y_diff_re(t0:t1,k,2)))/2.0_dp*im_unit
    z_aniso_resraman(0:t1-t0,k)=z_aniso_resraman(0:t1-t0,k)-((alpha_resraman_y_diff_im(t0,k,2)-alpha_resraman_z_diff_im(t0,k,3))&
            *(alpha_resraman_y_diff_re(t0:t1,k,2)-alpha_resraman_z_diff_re(t0:t1,k,3)))/2.0_dp*im_unit
    z_aniso_resraman(0:t1-t0,k)=z_aniso_resraman(0:t1-t0,k)-((alpha_resraman_z_diff_im(t0,k,3)-alpha_resraman_x_diff_im(t0,k,1))&
            *(alpha_resraman_z_diff_re(t0:t1,k,3)-alpha_resraman_x_diff_re(t0:t1,k,1)))/2.0_dp*im_unit
    z_aniso_resraman(0:t1-t0,k)=z_aniso_resraman(0:t1-t0,k)-((alpha_resraman_x_diff_im(t0,k,2)*0.50_dp+&
	alpha_resraman_y_diff_im(t0,k,1)&
            *0.50_dp)*(alpha_resraman_x_diff_re(t0:t1,k,2)*0.50_dp+alpha_resraman_y_diff_re(t0:t1,k,1)*0.50_dp))*3.0_dp*im_unit
    z_aniso_resraman(0:t1-t0,k)=z_aniso_resraman(0:t1-t0,k)-((alpha_resraman_y_diff_im(t0,k,3)*0.50_dp+&
	alpha_resraman_z_diff_im(t0,k,2)&
            *0.50_dp)*(alpha_resraman_y_diff_re(t0:t1,k,3)*0.50_dp+alpha_resraman_z_diff_re(t0:t1,k,2)*0.50_dp))*3.0_dp*im_unit
    z_aniso_resraman(0:t1-t0,k)=z_aniso_resraman(0:t1-t0,k)-((alpha_resraman_z_diff_im(t0,k,1)*0.50_dp+&
	alpha_resraman_x_diff_im(t0,k,3)&
            *0.50_dp)*(alpha_resraman_z_diff_re(t0:t1,k,1)*0.50_dp+alpha_resraman_x_diff_re(t0:t1,k,3)*0.50_dp))*3.0_dp*im_unit
    
  norm_aniso(0:t1-t0,k)=norm_aniso(0:t1-t0,k)+1.0_dp
    ENDDO
ENDDO

z_aniso_resraman(:,:)=REAL(z_aniso_resraman(:,:)/norm_aniso(:,:),kind=dp)
z_aniso_resraman(:,:)=REAL(z_aniso_resraman(:,:)/(2.0_dp*pi),kind=dp)

DO i=0,t_cor-1
   DO j=1,natom-1
       z_aniso_resraman(i,j)=z_aniso_resraman(i,j)*0.5_dp*(1+COS(2.0_dp*3.14_dp*i/(2.0_dp*(t_cor-1))))
   ENDDO
ENDDO

DO i=1,natom-1
      z_aniso_resraman(t_cor,i)=0.0_dp
ENDDO

DO i=1,t_cor-1
  DO j=1,natom-1
      z_aniso_resraman(t_cor+i,j)=z_aniso_resraman(t_cor-i,j)
  ENDDO
ENDDO

OPEN(UNIT=62,FILE='result_cvv_aniso_resraman.txt',STATUS='unknown',IOSTAT=stat)

DO i=0,2*t_cor-1
  DO j=1,natom-1
        WRITE(62,*) z_aniso_resraman(i,j)
  ENDDO
ENDDO
CLOSE(62)
DEALLOCATE(norm_aniso)

END SUBROUTINE cvv_resraman

!***********************************************************************************
!***********************************************************************************

SUBROUTINE cvv_only_x(mol_num,natom,framecount,t_cor,z_para,z_ortho,alpha_diff_x,&
                alpha_diff_y,alpha_diff_z,dt,pi,direction)

CHARACTER(LEN=40),INTENT(IN)                             :: direction
INTEGER,INTENT(INOUT)                                    :: natom,framecount,mol_num,t_cor
REAL(kind=dp),INTENT(IN)                                  :: dt,pi
REAL(kind=dp),DIMENSION(:),ALLOCATABLE,INTENT(OUT)        :: z_para,z_ortho
REAL(kind=dp),DIMENSION(:,:,:),ALLOCATABLE,INTENT(INOUT)  :: alpha_diff_x,alpha_diff_y,alpha_diff_z

CHARACTER(LEN=40)                                        :: chara
CHARACTER(LEN=2),DIMENSION(:),ALLOCATABLE                :: element
INTEGER,DIMENSION(:),ALLOCATABLE                         :: norm
INTEGER                                                  :: stat,i,j,k,m,t0,t1

ALLOCATE(z_para(0:t_cor*2),norm(0:t_cor*2))
ALLOCATE(z_ortho(0:t_cor*2))

framecount=framecount-2

 norm=0.0_dp
 z_para=0.0_dp
 z_ortho=0.0_dp
 DO t0=1,framecount
  t1=MIN(framecount,t0+t_cor)
   DO k=1,mol_num

    IF (direction=='1') THEN
       z_para(0:t1-t0)=z_para(0:t1-t0)+alpha_diff_x(t0,k,1)*alpha_diff_x(t0:t1,k,1)!*18.01468_dp
       z_ortho(0:t1-t0)=z_ortho(0:t1-t0)+alpha_diff_x(t0,k,2)*alpha_diff_x(t0:t1,k,2)!*18.01468_dp
    ELSE IF (direction=='2') THEN 
       z_para(0:t1-t0)=z_para(0:t1-t0)+alpha_diff_y(t0,k,2)*alpha_diff_y(t0:t1,k,2)
       z_ortho(0:t1-t0)=z_ortho(0:t1-t0)+alpha_diff_y(t0,k,3)*alpha_diff_y(t0:t1,k,3)
    ELSEIF(direction=='3') THEN
       z_para(0:t1-t0)=z_para(0:t1-t0)+alpha_diff_z(t0,k,3)*alpha_diff_z(t0:t1,k,3)
       z_ortho(0:t1-t0)=z_ortho(0:t1-t0)+alpha_diff_z(t0,k,1)*alpha_diff_z(t0:t1,k,1)
    ENDIF

  ENDDO
    norm(0:t1-t0)=norm(0:t1-t0)+1.0_dp
 ENDDO

z_para(:)=z_para(:)/norm(:)
z_para(:)=z_para(:)/mol_num
z_para(:)=z_para(:)/(2.0_dp*pi)
z_ortho(:)=z_ortho(:)/norm(:)
z_ortho(:)=z_ortho(:)/mol_num
z_ortho(:)=z_ortho(:)/(2.0_dp*pi)

DO i=0,t_cor-1
       z_para(i)=z_para(i)*0.5_dp*(1+COS(2.0_dp*3.14_dp*i/(2.0_dp*(t_cor-1))))
       z_ortho(i)=z_ortho(i)*0.5_dp*(1+COS(2.0_dp*3.14_dp*i/(2.0_dp*(t_cor-1))))
ENDDO

z_para(t_cor)=0.0_dp
z_ortho(t_cor)=0.0_dp

DO i=1,t_cor-1                                                                                                                    
        z_para(t_cor+i)=z_para(t_cor-i)
        z_ortho(t_cor+i)=z_ortho(t_cor-i)
ENDDO

OPEN(UNIT=60,FILE='result_cvv_para.txt',STATUS='unknown',IOSTAT=stat)

DO i=0,2*t_cor-1
        WRITE(60,*) z_para(i)
ENDDO
CLOSE(60)


OPEN(UNIT=61,FILE='result_cvv_ortho.txt',STATUS='unknown',IOSTAT=stat)

DO i=0,2*t_cor-1
        WRITE(61,*) z_ortho(i)
ENDDO
CLOSE(61)
DEALLOCATE(norm)
END SUBROUTINE cvv_only_x

END MODULE vel_cor

