MODULE calc_spectra

USE setup,                             ONLY: constants, read_input, conversion
USE read_traj,                         ONLY: read_coord_frame
USE fin_diff,                          ONLY: central_diff, forward_diff
USE vel_cor,                           ONLY: cvv,cvv_iso,cvv_aniso,cvv_only_x,cvv_resraman
USE dipole_calc,                       ONLY: center_mass,solv_frag_index,wannier_frag,wannier
USE pade,                              ONLY: interpolate

USE,INTRINSIC                              :: iso_c_binding
USE OMP_LIB

IMPLICIT NONE

INCLUDE 'fftw3.f03'

PUBLIC :: spec_power,spec_ir,spec_raman,normal_mode_analysis,spec_static_raman,spec_abs,spec_static_resraman,spec_resraman

CONTAINS
SUBROUTINE spec_power(z,zhat,type_input,freq_range,natom,framecount,t_cor,dt,element,filename,coord_v,v,input_mass,&
           dom,mass_atom,read_function,mol_num,mass_tot,pi,coord,system,frag_type)

        
CHARACTER(LEN=40),INTENT(INOUT)                          :: read_function,system,frag_type
CHARACTER(LEN=40),INTENT(INOUT)                          :: type_input,filename,input_mass
CHARACTER(LEN=2),DIMENSION(:),ALLOCATABLE,INTENT(INOUT)  :: element
INTEGER,INTENT(INOUT)                                    :: natom,framecount,t_cor,mol_num
REAL(KIND=8),INTENT(INOUT)                               :: dt,dom
REAL(KIND=8),INTENT(IN)                                  :: pi,freq_range
REAL(KIND=8),INTENT(OUT)                                 :: mass_tot
REAL(KIND=8),DIMENSION(:),ALLOCATABLE,INTENT(INOUT)      :: z,mass_atom
REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT)    :: coord
REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE,INTENT(INOUT)  :: coord_v,v
COMPLEX(KIND=8),DIMENSION(:),ALLOCATABLE,INTENT(OUT)     :: zhat

CHARACTER(LEN=40)                                        :: chara
INTEGER                                                  :: stat,i,j,k,m,t0,t1
INTEGER(KIND=8)                                          :: plan

ALLOCATE(zhat(0:t_cor*2-1))

zhat=COMPLEX(0.d0,0.0d0)

CALL read_coord_frame(natom,framecount,element,filename,coord_v)

IF (type_input=='1') THEN   !!If it is from positions, do finite differences first
 CALL central_diff(dt,natom,framecount,coord_v,v,read_function,mol_num,system)
 CALL cvv(natom,framecount,t_cor,v,z,type_input,dt,input_mass,mass_atom,mass_tot,pi,&
      mol_num,read_function,system,frag_type)

ELSEIF (type_input=='2') THEN   !!If it is from velocities, compute autcorrelation directly
 CALL cvv(natom,framecount,t_cor,coord_v,z,type_input,dt,input_mass,mass_atom,mass_tot,pi,&
      mol_num,read_function,system,frag_type)
ENDIF
  
CALL dfftw_plan_dft_r2c_1d(plan,2*t_cor,z,zhat,FFTW_ESTIMATE) !!!FFT
CALL dfftw_execute_dft_r2c(plan,z,zhat)
CALL dfftw_destroy_plan(plan)
       
zhat=REAL(zhat,KIND=8) 

OPEN(UNIT=63,FILE='power_spec.txt',STATUS='unknown',IOSTAT=stat) !!write the output
DO i=0,2*t_cor-1
 zhat(i)=(zhat(i)*dt*7.211349d-9)/(natom*3.0d0) !!!unit conversion
 IF ((i*freq_range).GE.5000d0) CYCLE
 WRITE(63,*) i*freq_range,REAL(zhat(i),KIND=8)
ENDDO

CLOSE(63)

END SUBROUTINE spec_power 

!***********************************************************************************************!
!***********************************************************************************************!
SUBROUTINE spec_ir(z,zhat,freq_range,natom,framecount,t_cor,dt,element,filename,coord_v,v,input_mass,&
           dom,pi,mol_num,box_all,box_x,box_y,box_z,vec,vec_pbc,debye,periodic,mass_tot,mass_atom,&
           type_input,dip,read_function,coord,type_dipole,dipole,system,mass_tot_frag,sinc_const,&
           nfrag,frag_type)
        
CHARACTER(LEN=40),INTENT(INOUT)                          :: filename,input_mass,periodic,system,frag_type
CHARACTER(LEN=40),INTENT(INOUT)                          :: type_dipole,read_function,type_input
CHARACTER(LEN=2),DIMENSION(:),ALLOCATABLE,INTENT(INOUT)  :: element
INTEGER,INTENT(INOUT)                                    :: natom,framecount,t_cor,mol_num,nfrag
REAL(KIND=8),INTENT(INOUT)                               :: dt,dom,vec(3),vec_pbc(3),mass_tot
REAL(KIND=8),INTENT(IN)                                  :: pi,freq_range,sinc_const
REAL(KIND=8),INTENT(INOUT)                               :: debye,box_x,box_y,box_z,box_all
REAL(KIND=8),DIMENSION(:),ALLOCATABLE,INTENT(INOUT)      :: z,mass_atom
REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT)    :: coord
REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE,INTENT(INOUT)  :: coord_v,v,dip,dipole
COMPLEX(KIND=8),DIMENSION(:),ALLOCATABLE,INTENT(OUT)     :: zhat
REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT)    :: mass_tot_frag

CHARACTER(LEN=40)                                        :: chara
INTEGER                                                  :: stat,i,j,k,m,t0,t1
INTEGER(KIND=8)                                          :: plan
REAL(KIND=8)                                             :: f

ALLOCATE(zhat(0:2*t_cor-1))

zhat = COMPLEX(0.d0,0.0d0)

IF (system=='2' .AND. type_dipole=='1') THEN
    mol_num=1
ENDIF

IF (system=='1' .OR. (system=='2' .AND. type_dipole=='1')) THEN  !!fragment approach or the whole cell
 CALL central_diff(dt,mol_num,framecount,dipole,v,read_function,mol_num,system)
 CALL cvv(nfrag,framecount,t_cor,v,z,type_input,dt,input_mass,mass_atom,mass_tot,pi,&
     mol_num,read_function,system,frag_type)

 !IF (type_dipole=='1') THEN !!Wannier centers, but currently only works for water molecule!
 ! CALL wannier(element,filename,natom,box_all,box_x,box_y,box_z,vec,vec_pbc,debye,mol_num,&
 !      periodic,mass_tot,framecount,mass_atom,coord_v,dip)
 ! CALL central_diff(dt,natom,framecount,dip,v,read_function,mol_num,system)

ELSEIF (system=='2' .AND. type_dipole=='2') THEN !!molecular approach 
  CALL read_coord_frame(mol_num,framecount,element,filename,coord_v)
  CALL central_diff(dt,mol_num,framecount,coord_v,v,read_function,mol_num,system)
  CALL cvv(mol_num,framecount,t_cor,v,z,type_input,dt,input_mass,mass_atom,mass_tot,pi,&
     mol_num,read_function,system,frag_type)
ENDIF

CALL dfftw_plan_dft_r2c_1d(plan,2*t_cor,z(0:2*t_cor-1),zhat(0:2*t_cor-1),FFTW_ESTIMATE) !!FFT!!
CALL dfftw_execute_dft_r2c(plan,z,zhat)
CALL dfftw_destroy_plan(plan)

zhat=REAL(zhat,KIND=8)
OPEN(UNIT=61,FILE='IR_spectrum.txt',STATUS='unknown',IOSTAT=stat) !!write output
DO i=0,2*t_cor-1
 zhat(i)=zhat(i)*3047.2310d0*dt*(sinc_const*(i)/SIN(sinc_const*(i)))**2.d0 !!unit conv. & sinc func.
 IF ((i*freq_range).GE.5000d0) CYCLE
 zhat(0)=0.00d0
 WRITE(61,*) i*freq_range,-1.0d0*REAL(zhat(i),KIND=8) 
ENDDO 
CLOSE(61)

!IF (type_dipole=='1') THEN
      !  DEALLOCATE(dip)
!ENDIF

END SUBROUTINE spec_ir
!****************************************************************************************!
!****************************************************************************************!
SUBROUTINE spec_raman(natom,framecount,element,coord,wannier_free,wannier_x,wannier_y,wannier_z,&
                mass_atom,mass_tot,periodic,mol_num,dt,dom,t_cor,speed_light,coord_v,v,type_input,&
                box_all,box_x,box_y,box_z,vec,vec_pbc,debye,read_function,z_iso,z_aniso,&
                z_ortho,z_para,const_planck,const_boltz,const_permit,temp,laser_in,pi,filename,&
                averaging,direction,type_dipole,system,natom_frag,fragment,refpoint,dipole,cell_type,&
                mass_tot_frag,frag_type,nfrag,charge,mass_tot_cell)
                
CHARACTER(LEN=40),INTENT(INOUT)                          :: read_function,type_input,periodic,filename,frag_type
CHARACTER(LEN=40),INTENT(INOUT)                          :: averaging,direction,type_dipole,system,cell_type
CHARACTER(LEN=40),INTENT(INOUT)                          :: wannier_free,wannier_x,wannier_y,wannier_z
CHARACTER(LEN=2),DIMENSION(:),ALLOCATABLE,INTENT(INOUT)  :: element
INTEGER,INTENT(INOUT)                                    :: nfrag,natom,framecount,t_cor,mol_num
INTEGER,DIMENSION(:),ALLOCATABLE,INTENT(INOUT)           :: natom_frag
INTEGER,DIMENSION(:,:,:),ALLOCATABLE,INTENT(INOUT)       :: fragment
REAL(KIND=8),INTENT(INOUT)                               :: dt,dom,vec(3),vec_pbc(3),debye,mass_tot,mass_tot_cell
REAL(KIND=8),INTENT(INOUT)                               :: const_planck,const_boltz,const_permit,temp,laser_in,pi
REAL(KIND=8),INTENT(INOUT)                               :: speed_light
REAL(KIND=8),INTENT(INOUT)                               :: box_x,box_y,box_z,box_all
REAL(KIND=8),DIMENSION(:),ALLOCATABLE,INTENT(INOUT)      :: z_iso,z_aniso,z_ortho,z_para,mass_atom,charge
REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT)    :: coord
REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE,INTENT(INOUT)  :: coord_v,v,dipole,refpoint
REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT)    :: mass_tot_frag

CHARACTER(LEN=40)                                        :: chara
INTEGER                                                  :: stat,i,j,k,m,t0,t1
INTEGER(KIND=8)                                          :: plan
INTEGER,DIMENSION(:),ALLOCATABLE                         :: natom_frag_x,natom_frag_free
INTEGER,DIMENSION(:),ALLOCATABLE                         :: natom_frag_y,natom_frag_z
INTEGER,DIMENSION(:,:,:),ALLOCATABLE                     :: fragment_x,fragment_free
INTEGER,DIMENSION(:,:,:),ALLOCATABLE                     :: fragment_y,fragment_z
COMPLEX(KIND=8),DIMENSION(:),ALLOCATABLE                 :: zhat_iso,zhat_aniso,zhat_para
COMPLEX(KIND=8),DIMENSION(:),ALLOCATABLE                 :: zhat_ortho,zhat_unpol
REAL(KIND=8),DIMENSION(:),ALLOCATABLE                    :: zhat_unpol_x,zhat_depol_x,zhat_para_all,zhat_depol
REAL(KIND=8)                                             :: f,freq_range
REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE                :: refpoint_free,refpoint_x,refpoint_y,refpoint_z
REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE                :: dip_free,dip_x,dip_y,dip_z
REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE                :: alpha_x,alpha_y,alpha_z
REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE                :: alpha_diff_x,alpha_diff_y,alpha_diff_z

!!!!ALLOCATION!!!
ALLOCATE(alpha_x(framecount,mol_num,3),alpha_y(framecount,mol_num,3),alpha_z(framecount,mol_num,3))

IF (averaging=='1') THEN

!!FIELD_FREE!!!
IF (system=='1' .OR. type_dipole=='1') THEN
 IF (cell_type.NE.'3') THEN
  CALL read_coord_frame(natom,framecount,element,wannier_free,coord_v)
  CALL center_mass(natom_frag,natom,refpoint,coord_v,wannier_free,element,box_all,box_x,box_y,&
        box_z,vec,vec_pbc,fragment_free,mass_atom,framecount,cell_type,mass_tot_frag,frag_type,mol_num,nfrag,&
        type_dipole,system,mass_tot_cell)
  CALL wannier_frag(natom_frag,wannier_free,natom,element,coord_v,box_all,box_x,box_y,box_z,vec,vec_pbc,dip_free,&
          refpoint,fragment_free,framecount,mass_tot_frag,mol_num,system,type_dipole,charge,mass_tot_cell)
 ELSEIF (cell_type=='3') THEN
  CALL read_coord_frame(natom,framecount,element,wannier_free,coord_v)
  CALL solv_frag_index(natom,coord_v,wannier_free,element,box_all,vec,vec_pbc,&
                     box_x,box_y,box_z,mass_atom,framecount,cell_type,refpoint,natom_frag_free,fragment_free,mass_tot_frag)
  CALL wannier_frag(natom_frag_free,wannier_free,natom,element,coord_v,box_all,box_x,box_y,box_z,vec,vec_pbc,dip_free,&
        refpoint,fragment_free,framecount,mass_tot_frag,mol_num,system,type_dipole,charge,mass_tot_cell)
 ENDIF
ELSEIF (system=='2') THEN
!IF (type_dipole=='1') THEN
!CALL read_coord_frame(natom,framecount,element,wannier_free,coord_v)
!CALL wannier(element,filename,natom,box_all,box_x,box_y,box_z,vec,vec_pbc,debye,mol_num,&
 !      periodic,mass_tot,framecount,mass_atom,coord_v,dip_free)
IF (type_dipole=='2') THEN
CALL read_coord_frame(mol_num,framecount,element,wannier_free,dip_free)
ENDIF
ENDIF

!!!X-FIELD!!!
CALL read_coord_frame(natom,framecount,element,wannier_x,coord_v)
IF (system=='1' .OR. type_dipole=='1') THEN
 IF (cell_type.NE.'3') THEN
  CALL center_mass(natom_frag,natom,refpoint,coord_v,wannier_x,element,box_all,box_x,box_y,&
        box_z,vec,vec_pbc,fragment_x,mass_atom,framecount,cell_type,mass_tot_frag,frag_type,mol_num,nfrag,&
        type_dipole,system,mass_tot_cell)
  CALL wannier_frag(natom_frag,wannier_x,natom,element,coord_v,box_all,box_x,box_y,box_z,vec,vec_pbc,dip_x,&
        refpoint,fragment_x,framecount,mass_tot_frag,mol_num,system,type_dipole,charge,mass_tot_cell)
 ELSEIF (cell_type=='3') THEN
  CALL solv_frag_index(natom,coord_v,wannier_x,element,box_all,vec,vec_pbc,&
                     box_x,box_y,box_z,mass_atom,framecount,cell_type,refpoint,natom_frag_x,fragment_x,mass_tot_frag)
  CALL wannier_frag(natom_frag_x,wannier_x,natom,element,coord_v,box_all,box_x,box_y,box_z,vec,vec_pbc,dip_x,&
        refpoint,fragment_x,framecount,mass_tot_frag,mol_num,system,type_dipole,charge,mass_tot_cell)
 ENDIF
 IF (system=='1') THEN
    CALL forward_diff(mol_num,framecount,alpha_x,dip_free,dip_x,system,read_function)
 ELSEIF (system=='2' .AND. type_dipole=='1') THEN
    CALL forward_diff(nfrag,framecount,alpha_x,dip_free,dip_x,system,read_function)
 ENDIF
 ELSEIF (system=='2') THEN
!IF (type_dipole=='1') THEN
!CALL wannier(element,wannier_x,natom,box_all,box_x,box_y,box_z,vec,vec_pbc,debye,mol_num,&
 !       periodic,mass_tot,framecount,mass_atom,coord_v,dip_x)
!CALL forward_diff(mol_num,framecount,alpha_x,dip_free,dip_x,system,read_function)

IF (type_dipole=='2') THEN
CALL forward_diff(mol_num,framecount,alpha_x,dip_free,coord_v,system,read_function)

ELSEIF (type_dipole=='3') THEN
DO i=1,framecount
 DO j=1,1
 alpha_x(i,j,:)=coord_v(i,j,:)
 ENDDO
ENDDO 
alpha_x=REAL(alpha_x*((8.988D+15)/(5.142D+11*3.33564D-30)),KIND=8) !conversion to debye/E
ENDIF
ENDIF

IF (system=='2' .AND. type_dipole=='1') THEN
  CALL central_diff(dt,nfrag,framecount,alpha_x,alpha_diff_x,read_function,mol_num,system)
ELSE  
  CALL central_diff(dt,mol_num,framecount,alpha_x,alpha_diff_x,read_function,mol_num,system)
 ENDIF

!!!Y-FIELD!!!
CALL read_coord_frame(natom,framecount,element,wannier_y,coord_v)
IF (system=='1' .OR. type_dipole=='1') THEN
 IF (cell_type.NE.'3') THEN
 CALL center_mass(natom_frag,natom,refpoint,coord_v,wannier_y,element,box_all,box_x,box_y,&
        box_z,vec,vec_pbc,fragment_y,mass_atom,framecount,cell_type,mass_tot_frag,frag_type,mol_num,nfrag,&
        type_dipole,system,mass_tot_cell)
 CALL wannier_frag(natom_frag,wannier_y,natom,element,coord_v,box_all,box_x,box_y,box_z,vec,vec_pbc,dip_y,&
        refpoint,fragment_y,framecount,mass_tot_frag,mol_num,system,type_dipole,charge,mass_tot_cell)
 ELSEIF (cell_type=='3') THEN
  CALL solv_frag_index(natom,coord_v,wannier_y,element,box_all,vec,vec_pbc,&
                     box_x,box_y,box_z,mass_atom,framecount,cell_type,refpoint,natom_frag_y,fragment_y,mass_tot_frag)
  CALL wannier_frag(natom_frag_y,wannier_y,natom,element,coord_v,box_all,box_x,box_y,box_z,vec,vec_pbc,dip_y,&
        refpoint,fragment_y,framecount,mass_tot_frag,mol_num,system,type_dipole,charge,mass_tot_cell)
  ENDIF
 IF (system=='1') THEN
    CALL forward_diff(mol_num,framecount,alpha_y,dip_free,dip_y,system,read_function)
 ELSEIF (system=='2' .AND. type_dipole=='1') THEN
    CALL forward_diff(nfrag,framecount,alpha_y,dip_free,dip_y,system,read_function)
 ENDIF
ELSEIF (system=='2') THEN
!IF (type_dipole=='1') THEN
!CALL wannier(element,wannier_y,natom,box_all,box_x,box_y,box_z,vec,vec_pbc,debye,mol_num,&
 !        periodic,mass_tot,framecount,mass_atom,coord_v,dip_y)
!CALL forward_diff(mol_num,framecount,alpha_y,dip_free,dip_y,system,read_function)

IF (type_dipole=='2') THEN
CALL forward_diff(mol_num,framecount,alpha_y,dip_free,coord_v,system,read_function)

ELSEIF (type_dipole=='3') THEN
DO i=1,framecount
 DO j=1,1
 alpha_y(i,j,:)=coord_v(i,j,:)
 ENDDO
ENDDO 
alpha_y=REAL(alpha_y*((8.988D+15)/(5.142D+11*3.33564D-30)),KIND=8) !conversion to debye/E
ENDIF
ENDIF

IF (system=='2' .AND. type_dipole=='1') THEN
  CALL central_diff(dt,nfrag,framecount,alpha_y,alpha_diff_y,read_function,mol_num,system)
ELSE  
  CALL central_diff(dt,mol_num,framecount,alpha_y,alpha_diff_y,read_function,mol_num,system)
 ENDIF

!!!Z-FIELD!!! 
CALL read_coord_frame(natom,framecount,element,wannier_z,coord_v)
IF (system=='1' .OR. type_dipole=='1') THEN
 IF (cell_type.NE.'3') THEN
 CALL center_mass(natom_frag,natom,refpoint,coord_v,wannier_z,element,box_all,box_x,box_y,&
        box_z,vec,vec_pbc,fragment_z,mass_atom,framecount,cell_type,mass_tot_frag,frag_type,mol_num,nfrag,&
        type_dipole,system,mass_tot_cell)
 CALL wannier_frag(natom_frag,wannier_z,natom,element,coord_v,box_all,box_x,box_y,box_z,vec,vec_pbc,dip_z,&
        refpoint,fragment_z,framecount,mass_tot_frag,mol_num,system,type_dipole,charge,mass_tot_cell)
 ELSEIF (cell_type=='3') THEN
  CALL solv_frag_index(natom,coord_v,wannier_z,element,box_all,vec,vec_pbc,&
                     box_x,box_y,box_z,mass_atom,framecount,cell_type,refpoint,natom_frag_z,fragment_z,mass_tot_frag)
  CALL wannier_frag(natom_frag_z,wannier_z,natom,element,coord_v,box_all,box_x,box_y,box_z,vec,vec_pbc,dip_z,&
        refpoint,fragment_z,framecount,mass_tot_frag,mol_num,system,type_dipole,charge,mass_tot_cell)
 ENDIF
 IF (system=='1') THEN
    CALL forward_diff(mol_num,framecount,alpha_z,dip_free,dip_z,system,read_function)
 ELSEIF (system=='2' .AND. type_dipole=='1') THEN
    CALL forward_diff(nfrag,framecount,alpha_z,dip_free,dip_z,system,read_function)
 ENDIF

ELSEIF (system=='2') THEN
!IF (type_dipole=='1') THEN
!CALL wannier(element,wannier_z,natom,box_all,box_x,box_y,box_z,vec,vec_pbc,debye,mol_num,&
 !        periodic,mass_tot,framecount,mass_atom,coord_v,dip_z)
!CALL forward_diff(mol_num,framecount,alpha_z,dip_free,dip_z,system,read_function)

IF (type_dipole=='2') THEN
 CALL forward_diff(mol_num,framecount,alpha_z,dip_free,coord_v,system,read_function)

ELSEIF (type_dipole=='3') THEN
DO i=1,framecount
 DO j=1,1
 alpha_z(i,j,:)=coord_v(i,j,:)
 ENDDO
ENDDO 
alpha_z=REAL(alpha_z*((8.988D+15)/(5.142D+11*3.33564D-30)),KIND=8) !conversion to debye/E
ENDIF
ENDIF

IF (system=='2' .AND. type_dipole=='1') THEN
  CALL central_diff(dt,nfrag,framecount,alpha_z,alpha_diff_z,read_function,mol_num,system)
ELSE  
  CALL central_diff(dt,mol_num,framecount,alpha_z,alpha_diff_z,read_function,mol_num,system)
 ENDIF

!!!ACF AND FFT CALC!!!
print*,nfrag,'nfrag check'
        ALLOCATE(zhat_iso(0:t_cor*2),zhat_aniso(0:t_cor*2))
        ALLOCATE(zhat_unpol(0:t_cor*2),zhat_depol(0:t_cor*2),zhat_para_all(0:t_cor*2))
        
        zhat_iso = COMPLEX(0.d0,0.0d0)
        zhat_aniso = COMPLEX(0.d0,0.0d0)
        zhat_unpol = COMPLEX(0.d0,0.0d0)
        zhat_depol = COMPLEX(0.d0,0.0d0)
       
        IF (system=='1' .OR. (system=='2' .AND. type_dipole=='1')) THEN 
            CALL cvv_iso(nfrag,framecount,t_cor,z_iso,alpha_diff_x,alpha_diff_y,alpha_diff_z,dt,pi,frag_type)
        ELSE
            CALL cvv_iso(mol_num,framecount,t_cor,z_iso,alpha_diff_x,alpha_diff_y,alpha_diff_z,dt,pi,frag_type)
        ENDIF

        CALL dfftw_plan_dft_r2c_1d(plan,2*t_cor,z_iso,zhat_iso,FFTW_ESTIMATE)
        CALL dfftw_execute_dft_r2c(plan,z_iso,zhat_iso)
   
        IF (system=='1' .OR. (system=='2' .AND. type_dipole=='1')) THEN 
            CALL cvv_aniso(nfrag,natom,framecount,t_cor,z_aniso,alpha_diff_x,alpha_diff_y,&
                  alpha_diff_z,dt,pi,frag_type)
        ELSE
            CALL cvv_aniso(mol_num,natom,framecount,t_cor,z_aniso,alpha_diff_x,alpha_diff_y,&
                  alpha_diff_z,dt,pi,frag_type)
        ENDIF 
        
        CALL dfftw_plan_dft_r2c_1d(plan,2*t_cor,z_aniso,zhat_aniso,FFTW_ESTIMATE)
        CALL dfftw_execute_dft_r2c(plan,z_aniso,zhat_aniso)
        
        !!!ORTHOGONAL!!!
        OPEN(UNIT=63,FILE='result_fft_water_lib_ortho.txt',STATUS='unknown',IOSTAT=stat)
        zhat_aniso=REAL(zhat_aniso,KIND=8) 
        freq_range=REAL(dom/(2.0d0*t_cor),KIND=8)
        f=freq_range*dt*1.883652d-4

        DO i=0,2*t_cor-2
        zhat_aniso(i+1)=REAL(zhat_aniso(i+1),KIND=8)*(f*(i+1)/SIN(f*(i+1)))**2.d0
        zhat_aniso(i)=zhat_aniso(i)*((const_planck)/(8.0d0*const_boltz*const_permit*const_permit)*1d-29*0.421d0*dt&
        *(((laser_in-((i)*freq_range))**4)/((i)*freq_range))*(1.0d0/(1.0d0-EXP((-1.438777d0*&
           ((i)*freq_range))/temp))))*2.0d0*2.0d0*pi!*((-1.438777d0*i*freq_range)/temp)*(1.0d0/(1.0d0-EXP((-1.438777d0*&
                      !((i)*freq_range))/temp)))
                
        zhat_aniso(0)=0.0d0
        IF ((i*freq_range).GE.5000.0d0) CYCLE
        WRITE(63,*) i*freq_range,((REAL(zhat_aniso(i),KIND=8))/15.0d0)
        ENDDO
        CLOSE(63)

        !!!PARALLEL!!!
        OPEN(UNIT=64,FILE='result_fft_water_lib_para.txt',STATUS='unknown',IOSTAT=stat)
        zhat_iso=REAL(zhat_iso,KIND=8) 
        zhat_para_all=REAL(zhat_para_all,KIND=8) 
        freq_range=REAL(dom/(2.0d0*t_cor),KIND=8)
        f=freq_range*dt*1.883652d-4
           
        DO i=0,2*t_cor-2
 
         zhat_iso(i+1)=REAL(zhat_iso(i+1),KIND=8)*(f*(i+1)/SIN(f*(i+1)))**2.d0 
         zhat_iso(i)=zhat_iso(i)*((const_planck)/(8.0d0*const_boltz*const_permit*const_permit)*1d-29*0.421d0*dt&
         *(((laser_in-((i)*freq_range))**4)/((i)*freq_range))*(1.0d0/(1.0d0-EXP((-1.438777d0*&
          ((i)*freq_range))/temp))))*2.0d0*2.0d0*pi
  
        zhat_para_all(i)=zhat_iso(i)+(zhat_aniso(i)*4.0d0/45.0d0)
        zhat_para_all(0)=0.0d0
        IF ((i*freq_range).GE.5000.0d0) CYCLE
        WRITE(64,*) i*freq_range,REAL(zhat_para_all(i),KIND=8)
        ENDDO
        CLOSE(64)

        !!!UNPOL!!!
        OPEN(UNIT=65,FILE='result_fft_water_lib_unpol.txt',STATUS='unknown',IOSTAT=stat)
        zhat_unpol=REAL(zhat_unpol,KIND=8) 
        freq_range=REAL(dom/(2.0d0*t_cor),KIND=8)

        DO i=0,2*t_cor-2
         zhat_unpol(i)=zhat_iso(i)+(zhat_aniso(i)*7.0d0/45.0d0)
         zhat_unpol(0)=0.00d0
         IF ((i*freq_range).GE.5000.0d0) CYCLE
         WRITE(65,*) i*freq_range,REAL(zhat_unpol(i),KIND=8)
        ENDDO
        CLOSE(65)

        !!!DEPOL RATIO!!!
        OPEN(UNIT=66,FILE='result_fft_water_lib_depol.txt',STATUS='unknown',IOSTAT=stat)
        zhat_depol=REAL(zhat_depol,KIND=8) 
        freq_range=REAL(dom/(2.0d0*t_cor),KIND=8)

        DO i=0,2*t_cor-2        
         zhat_depol(i)=(REAL(zhat_aniso(i),KIND=8)/15.0d0)/REAL(zhat_para_all(i),KIND=8)
         IF ((i*freq_range).GE.5000.0d0) CYCLE
         WRITE(66,*) i*freq_range,REAL(zhat_depol(i),KIND=8)
        ENDDO

        CLOSE(66)
        DEALLOCATE(z_iso,z_aniso)
        DEALLOCATE(zhat_depol,zhat_para_all,zhat_unpol)

ELSEIF (averaging=='2') THEN

IF (type_dipole=='1') THEN
CALL read_coord_frame(natom,framecount,element,wannier_free,coord_v)
CALL wannier(element,filename,natom,box_all,box_x,box_y,box_z,vec,vec_pbc,debye,mol_num,&
       periodic,mass_tot,framecount,mass_atom,coord_v,dip_free)
ELSEIF (type_dipole=='2') THEN
CALL read_coord_frame(natom,framecount,element,wannier_free,dip_free)
ENDIF

!!!X-FIELD!!!
        IF (direction=='1') THEN
        !!!X-FIELD!!!
        CALL read_coord_frame(natom,framecount,element,wannier_x,coord_v)
        IF (type_dipole=='1') THEN
        CALL wannier(element,wannier_x,natom,box_all,box_x,box_y,box_z,vec,vec_pbc,debye,mol_num,&
          periodic,mass_tot,framecount,mass_atom,coord_v,dip_x)
        CALL forward_diff(mol_num,framecount,alpha_x,dip_free,dip_x,system,read_function)
        ELSEIF (type_dipole=='2') THEN
        CALL forward_diff(mol_num,framecount,alpha_x,dip_free,coord_v,system,read_function)
        ENDIF
        CALL central_diff(dt,natom,framecount,alpha_x,alpha_diff_x,read_function,mol_num,system)
        !!!Y-FIELD!!!
       ELSEIF (direction=='2') THEN
        CALL read_coord_frame(natom,framecount,element,wannier_y,coord_v)
        IF (type_dipole=='1') THEN
        CALL wannier(element,wannier_y,natom,box_all,box_x,box_y,box_z,vec,vec_pbc,debye,mol_num,&
         periodic,mass_tot,framecount,mass_atom,coord_v,dip_y)
        CALL forward_diff(mol_num,framecount,alpha_y,dip_free,dip_y,system,read_function)
        ELSEIF (type_dipole=='2') THEN
        CALL forward_diff(mol_num,framecount,alpha_y,dip_free,coord_v,system,read_function)
        ENDIF
        CALL central_diff(dt,natom,framecount,alpha_y,alpha_diff_y,read_function,mol_num,system)
       ELSEIF (direction=='3') THEN
        !!!Z-FIELD!!! 
        CALL read_coord_frame(natom,framecount,element,wannier_z,coord_v)
        IF (type_dipole=='1') THEN
        CALL wannier(element,wannier_z,natom,box_all,box_x,box_y,box_z,vec,vec_pbc,debye,mol_num,&
         periodic,mass_tot,framecount,mass_atom,coord_v,dip_z)
        CALL forward_diff(mol_num,framecount,alpha_z,dip_free,dip_z,system,read_function)
        ELSEIF (type_dipole=='2') THEN
        CALL forward_diff(mol_num,framecount,alpha_z,dip_free,coord_v,system,read_function)
        ENDIF
        CALL central_diff(dt,natom,framecount,alpha_z,alpha_diff_z,read_function,mol_num,system)
     
       ENDIF
        
       ALLOCATE(zhat_para(0:t_cor*2),zhat_unpol_x(0:t_cor*2),zhat_ortho(0:t_cor*2),zhat_depol_x(0:t_cor*2))

        zhat_para = COMPLEX(0.d0,0.0d0)
        zhat_ortho = COMPLEX(0.d0,0.0d0)
        zhat_unpol_x = COMPLEX(0.d0,0.0d0)
        zhat_depol_x = COMPLEX(0.d0,0.0d0)

        !!IF ONLY ISOTROPIC AVERAGING IS CONSIDERED!!
        CALL cvv_only_x(mol_num,natom,framecount,t_cor,z_para,z_ortho,alpha_diff_x,&
                alpha_diff_y,alpha_diff_z,dt,pi,direction)

        CALL dfftw_plan_dft_r2c_1d(plan,2*t_cor,z_para,zhat_para,FFTW_ESTIMATE)
        CALL dfftw_execute_dft_r2c(plan,z_para,zhat_para)

        CALL dfftw_plan_dft_r2c_1d(plan,2*t_cor,z_ortho,zhat_ortho,FFTW_ESTIMATE)
        CALL dfftw_execute_dft_r2c(plan,z_ortho,zhat_ortho)


        !!ORTHOGONAL!!
        OPEN(UNIT=68,FILE='result_fft_water_lib_ortho_iso.txt',STATUS='unknown',IOSTAT=stat)
        zhat_ortho=REAL(zhat_ortho,KIND=8) 
        freq_range=REAL(dom/t_cor,KIND=8)
        f=freq_range*dt*1.883652d-4
         
        DO i=0,2*t_cor-2
         zhat_ortho(i+1)=REAL(zhat_ortho(i+1),KIND=8)*(f*(i+1)/SIN(f*(i+1)))**2.d0

         zhat_ortho(i)=REAL(zhat_ortho(i),KIND=8)*((const_planck)/(8.0d0*const_boltz*const_permit*const_permit)*1d-29*0.421d0*dt&
                  *(((laser_in-((i)*freq_range))**4)/((i)*freq_range))*(1.0d0/(1.0d0-EXP((-1.438777d0*&
                ((i)*freq_range))/temp))))*2.0d0*pi*2.0d0
     
         zhat_ortho(0)=0.0d0
        IF ((i*freq_range).GE.5000d0) CYCLE
        WRITE(68,*) i*freq_range,(REAL(zhat_ortho(i),KIND=8))
        ENDDO
        CLOSE(68)

        !!PARALLEL!!
        OPEN(UNIT=67,FILE='result_fft_water_lib_para_iso.txt',STATUS='unknown',IOSTAT=stat)
        zhat_para=REAL(zhat_para,KIND=8) 
        freq_range=REAL(dom/t_cor,KIND=8)
        f=freq_range*dt*1.883652d-4

        DO i=0,2*t_cor-2
         zhat_para(i+1)=REAL(zhat_para(i+1),KIND=8)*(f*(i+1)/SIN(f*(i+1)))**2.d0

         zhat_para(i)=REAL(zhat_para(i),KIND=8)*((const_planck)/(8.0d0*const_boltz*const_permit*const_permit)*1d-29*0.421d0*dt&
                 *(((laser_in-((i)*freq_range))**4)/((i)*freq_range))*(1.0d0/(1.0d0-EXP((-1.438777d0*&
                 ((i)*freq_range))/temp))))*2.0d0*pi*2.0d0    
         zhat_para(0)=0.0d0
        IF ((i*freq_range).GE.5000d0) CYCLE
        WRITE(67,*) (i)*freq_range,(REAL(zhat_para(i),KIND=8))!,REAL(integral(i),KIND=8)
        ENDDO
        CLOSE(67)

        !!UNPOL!!
        OPEN(UNIT=69,FILE='result_fft_water_lib_unpol_iso.txt',STATUS='unknown',IOSTAT=stat)
        freq_range=REAL(dom/t_cor,KIND=8)

        DO i=0,2*t_cor-2       
          zhat_unpol_x(i)=zhat_para(i)+zhat_ortho(i) 
          IF ((i*freq_range).GE.5000d0) CYCLE
          WRITE(69,*) i*freq_range,REAL(zhat_unpol_x(i),KIND=8)
        ENDDO
        CLOSE(69)

       !!DEPOL RATIO!!
       OPEN(UNIT=70,FILE='result_fft_water_lib_depol_iso.txt',STATUS='unknown',IOSTAT=stat)

       DO i=0,2*t_cor-2
         zhat_depol_x(i)=REAL(zhat_ortho(i),KIND=8)/REAL(zhat_para(i),KIND=8) 
         IF ((i*freq_range).GE.5000d0) CYCLE
         WRITE(70,*) i*freq_range,REAL(zhat_depol_x(i),KIND=8)!REAL(zhat_ortho(i),KIND=8)/REAL(zhat_para(i),KIND=8)
       ENDDO
       CLOSE(70)

       DEALLOCATE(z_para,z_ortho,zhat_unpol_x,zhat_depol_x,zhat_para,zhat_ortho)
ENDIF

CALL dfftw_destroy_plan(plan)

!DEALLOCATE(dip_free,dip_x,dip_y,dip_z)
IF (system=='1') THEN
    DEALLOCATE(fragment_x,fragment_y,fragment_z,fragment_free)
  !  DEALLOCATE(natom_frag_x,natom_frag_y,natom_frag_z,natom_frag_free)
ENDIF

DEALLOCATE(alpha_x,alpha_y,alpha_z)
DEALLOCATE(alpha_diff_x,alpha_diff_y,alpha_diff_z)

END SUBROUTINE spec_raman

!....................................................................................................................!
!....................................................................................................................!

SUBROUTINE normal_mode_analysis(natom,force,dx,hartreebohr2evang,hessian_factor,mass_mat)

INTEGER,INTENT(IN)                                          :: natom
REAL(KIND=8),INTENT(IN)                                     :: dx,hartreebohr2evang,hessian_factor
REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE,INTENT(IN)          :: mass_mat
REAL(KIND=8),DIMENSION(:,:,:,:,:),ALLOCATABLE,INTENT(INOUT) :: force

INTEGER                                                     :: i,j,m,n,p,k
REAL(KIND=8)                                                :: factor
REAL(KIND=8),DIMENSION(:,:,:,:),ALLOCATABLE                 :: hessian
REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE                     :: hessian_new

factor=REAL(1.0d0/(2.0d0*dx),KIND=8)

ALLOCATE(hessian(0:natom-1,0:2,0:natom-1,0:2),hessian_new(0:natom*3-1,0:natom*3-1))

hessian=hartreebohr2evang*factor*hessian_factor*(force(2,:,:,:,:)-force(1,:,:,:,:))

p=0
DO i=0,natom-1
 DO m=0,2
  k=0
  DO j=0,natom-1
   DO n=0,2
    hessian_new(i+m+p,j+n+k)=mass_mat(i+1,j+1)*hessian(i,m,j,n) 
   ENDDO
   k=k+2
  ENDDO
 ENDDO
 p=p+2
ENDDO

!hessian_new(:,:)=RESHAPE(hessian(:,:,:,:), (/3*natom, 3*natom/))
hessian_new(:,:)=REAL((hessian_new(:,:)+TRANSPOSE(hessian_new(:,:)))/2.0d0,KIND=8)

print*,hessian_new(2,8)
print*,hessian_new(4,5)

END SUBROUTINE normal_mode_analysis

!....................................................................................................................!
!....................................................................................................................!

SUBROUTINE spec_static_raman(nmodes,pol_dq,laser_in,freq,temp,raman_int,pi,element,coord,disp,bohr2ang,natom)

INTEGER,INTENT(INOUT)                                    :: nmodes,natom
CHARACTER(LEN=2),DIMENSION(:),ALLOCATABLE,INTENT(INOUT)  :: element
REAL(KIND=8),INTENT(INOUT)                               :: laser_in,temp,pi,bohr2ang
REAL(KIND=8),DIMENSION(:),ALLOCATABLE,INTENT(INOUT)      :: freq
REAL(KIND=8),DIMENSION(:),ALLOCATABLE,INTENT(OUT)        :: raman_int
REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT)    :: coord
REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE,INTENT(INOUT)  :: pol_dq
REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE,INTENT(INOUT)  :: disp

INTEGER                                                  :: stat,i,j,k,m,x,freq_range
INTEGER                                                  :: start_freq,end_freq
REAL(KIND=8),DIMENSION(:),ALLOCATABLE                    :: iso_sq,aniso_sq,data1,data2!,broad
REAL(KIND=8)                                             :: omega,broad

ALLOCATE(iso_sq(nmodes),aniso_sq(nmodes))
ALLOCATE(raman_int(nmodes))

start_freq=1.0d0
end_freq=INT(MAXVAL(freq)+1000.0D0)
freq_range=INT(end_freq-start_freq)
omega=5.0D0
data1=0.0d0
data2=0.0d0

ALLOCATE(data1(freq_range*nmodes),data2(freq_range*nmodes))

!!!Isotropic and anisotropic contributions!!
DO k=1,nmodes
 iso_sq(k)=REAL((pol_dq(k,1,1)+pol_dq(k,2,2)+pol_dq(k,3,3))/3.0d0,KIND=8)**2.0d0
     
 aniso_sq(k)=(0.50d0*(((pol_dq(k,1,1)-pol_dq(k,2,2))**2.0d0)+((pol_dq(k,2,2)-pol_dq(k,3,3))**2.0d0)&
 +((pol_dq(k,3,3)-pol_dq(k,1,1))**2.0d0)))+(3.0d0*((pol_dq(k,1,2)**2.0d0)+(pol_dq(k,2,3)**2.0d0)&
 +(pol_dq(k,3,1)**2.0d0)))
ENDDO

!!!Calculation of the intensities!!
DO k=1,nmodes
 raman_int(k)=REAL(((7.0d0*aniso_sq(k))+(45.0d0*iso_sq(k)))/45.0d0,KIND=8)*REAL(((laser_in-freq(k))**4.0d0)/freq(k),KIND=8)&
 *REAL(1.0d0/(1.0d0-EXP(-1.438777d0*freq(k)/temp)),KIND=8)
ENDDO

!!!Broadening the spectrum!!
DO x=start_freq,end_freq
 broad = 0.0d0
 DO i=1,nmodes
  broad=broad+(raman_int(i)*(1.0d0/(omega*SQRT(2.0d0*pi)))*EXP(-0.50d0*((x-freq(i))/omega)**2.0d0))
 ENDDO
 data2(x)=data2(x)+broad
ENDDO
       
OPEN(UNIT=98,FILE='result_static_raman.txt',STATUS='unknown',IOSTAT=stat)
DO i=start_freq,end_freq
 WRITE(98,*) i,data2(i)
ENDDO
CLOSE(98)

!!Write Molden output
raman_int=REAL(raman_int/MINVAL(raman_int),KIND=8)
OPEN(UNIT=15,FILE='raman.mol',STATUS='unknown',IOSTAT=stat)
WRITE(15,*) "[Molden Format]"
WRITE(15,*) "[GEOMETRIES] XYZ"
WRITE(15,*) natom
WRITE(15,*) 
DO i=1,natom
 WRITE(15,*) element(i),coord(i,1),coord(i,2),coord(i,3)
ENDDO
WRITE(15,*) "[FREQ]"
DO i=1,nmodes
 WRITE(15,*) freq(i)
ENDDO
WRITE(15,*) "[INT]"
DO i=1,nmodes
 WRITE(15,*) raman_int(i)
ENDDO
WRITE(15,*) "[FR-COORD]"
WRITE(15,*) natom
WRITE(15,*) 
DO i=1,natom
 WRITE(15,*) element(i),coord(i,1)/bohr2ang,coord(i,2)/bohr2ang,coord(i,3)/bohr2ang
ENDDO
WRITE(15,*) "[FR-NORM-COORD]"
DO i=1,nmodes
 WRITE(15,*) "vibration", i 
 DO j=1,natom
  WRITE(15,*) disp(j,i,1)/bohr2ang,disp(j,i,2)/bohr2ang,disp(j,i,3)/bohr2ang
 ENDDO
ENDDO
CLOSE(15)

DEALLOCATE(iso_sq,aniso_sq,data1,data2,raman_int)

END SUBROUTINE spec_static_raman
!....................................................................................................................!
!....................................................................................................................!

SUBROUTINE spec_abs(nmodes,natom,pol_rtp,freq,pi,framecount_rtp,speed_light,framecount_rtp_pade,reccm2ev,check_pade&
           ,dom_rtp)

INTEGER,INTENT(INOUT)                                         :: natom,nmodes,framecount_rtp,framecount_rtp_pade
CHARACTER(LEN=40),INTENT(IN)                                  :: check_pade
REAL(KIND=8),INTENT(IN)                                       :: speed_light,reccm2ev,pi,dom_rtp
REAL(KIND=8),DIMENSION(:),ALLOCATABLE,INTENT(INOUT)           :: freq
REAL(KIND=8),DIMENSION(:,:,:,:,:,:),ALLOCATABLE,INTENT(INOUT) :: pol_rtp

INTEGER                                                       :: stat,i,j,k,m,x,freq_range,l,o,n,r
INTEGER(KIND=8)                                               :: plan
REAL(KIND=8)                                                  :: rtp_freq_range
REAL(KIND=8),DIMENSION(:,:,:,:),ALLOCATABLE                   :: trace,abs_intens
COMPLEX(KIND=8),DIMENSION(:,:,:,:,:,:),ALLOCATABLE            :: y_out,zhat_pol_rtp

ALLOCATE(zhat_pol_rtp(natom,3,2,3,3,framecount_rtp))
ALLOCATE(y_out(natom,3,2,3,3,framecount_rtp_pade)) 

zhat_pol_rtp = COMPLEX(0.d0,0.0d0)

!!!FFT of the RTP polarizabilities
DO j=1,natom
 DO i=1,3
  DO k=1,2
   DO m=1,3
    DO o=1,3
     CALL dfftw_plan_dft_r2c_1d(plan,framecount_rtp,pol_rtp(j,i,k,m,o,1:framecount_rtp),&
     zhat_pol_rtp(j,i,k,m,o,1:framecount_rtp),FFTW_ESTIMATE)
     CALL dfftw_execute_dft_r2c(plan,pol_rtp(j,i,k,m,o,1:framecount_rtp),zhat_pol_rtp(j,i,k,m,o,1:framecount_rtp))
     CALL dfftw_destroy_plan(plan)
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO

IF (check_pade=='y') THEN
!!Call Pade
 !$OMP PARALLEL DO COLLAPSE(5)
 DO j=1,natom !!natom
  DO i=1,3 !!dims
   DO k=1,2 !!disp
    DO m=1,3
     DO o=1,3
      CALL interpolate(framecount_rtp,zhat_pol_rtp(j,i,k,m,o,1:framecount_rtp),framecount_rtp_pade,y_out(j,i,k,m,o,:))
     ENDDO
    ENDDO
   ENDDO 
  ENDDO 
 ENDDO
 !$OMP END PARALLEL DO
framecount_rtp=framecount_rtp_pade
ENDIF
print*,framecount_rtp,"checkpoint"
!!!Finding laser frequency
rtp_freq_range=REAL(dom_rtp/framecount_rtp,KIND=8)

ALLOCATE(trace(natom,3,2,framecount_rtp))
ALLOCATE(abs_intens(natom,3,2,framecount_rtp))
trace=0.0d0

!!!Calculate absorption spectra
DO j=1,natom !!atom_num
 DO i=1,3 !!dims
  DO k=1,2 !!direction
   DO o=1,framecount_rtp
    IF (check_pade=='n') THEN
     trace(j,i,k,o)=DIMAG(zhat_pol_rtp(j,i,k,1,1,o))+DIMAG(zhat_pol_rtp(j,i,k,2,2,o))+DIMAG(zhat_pol_rtp(j,i,k,3,3,o))
    ELSEIF (check_pade=='y') THEN
     trace(j,i,k,o)=DIMAG(y_out(j,i,k,1,1,o))+DIMAG(y_out(j,i,k,2,2,o))+DIMAG(y_out(j,i,k,3,3,o))
    ENDIF
    abs_intens(j,i,k,o)=(4.0d0*pi*trace(j,i,k,o))/(3.0d0*speed_light)
   ENDDO
  ENDDO
 ENDDO
ENDDO

OPEN(UNIT=13,FILE='absorption_spectra_pade.txt',STATUS='unknown',IOSTAT=stat)
DO j=1,1 !!atom_num: 1st atom
 DO i=1,1 !!dims: x dimension
  DO k=1,1 !! + direction
   DO o=1,framecount_rtp
    WRITE(13,*) o*rtp_freq_range*reccm2ev,abs_intens(j,i,k,o)*o*rtp_freq_range*(-1.0d0)
   ENDDO
  ENDDO
 ENDDO
ENDDO
CLOSE(13)

DEALLOCATE(zhat_pol_rtp,trace,abs_intens,y_out)

END SUBROUTINE spec_abs
!....................................................................................................................!
!....................................................................................................................!
SUBROUTINE spec_static_resraman(nmodes,natom,pol_rtp,laser_in_resraman,freq,temp,pi,framecount_rtp,dom_rtp,dx,bohr2ang,&
           disp,speed_light,framecount_rtp_pade,check_pade,atom_mass_inv_sqrt)

CHARACTER(LEN=40),INTENT(IN)                                  :: check_pade
INTEGER,INTENT(INOUT)                                         :: natom,nmodes,framecount_rtp,framecount_rtp_pade
REAL(KIND=8),INTENT(IN)                                       :: dx,bohr2ang,speed_light
REAL(KIND=8),INTENT(INOUT)                                    :: laser_in_resraman,temp,pi,dom_rtp
REAL(KIND=8),DIMENSION(:),ALLOCATABLE,INTENT(INOUT)           :: freq,atom_mass_inv_sqrt
REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE,INTENT(INOUT)       :: disp
REAL(KIND=8),DIMENSION(:,:,:,:,:,:),ALLOCATABLE,INTENT(INOUT) :: pol_rtp

INTEGER                                                       :: stat,i,j,k,m,x,freq_range,l,o,n,r
INTEGER                                                       :: start_freq,end_freq,rtp_point
INTEGER(KIND=8)                                               :: plan
REAL(KIND=8)                                                  :: omega,broad,factor
REAL(KIND=8)                                                  :: rtp_freq_range,pade_freq_range
REAL(KIND=8),DIMENSION(:),ALLOCATABLE                         :: data1,data2
REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE                       :: iso_sq,aniso_sq,raman_int
REAL(KIND=8),DIMENSION(:,:,:,:),ALLOCATABLE                   :: zhat_pol_dq_rtp
REAL(KIND=8),DIMENSION(:,:,:,:,:),ALLOCATABLE                 :: zhat_pol_dxyz_rtp
COMPLEX(KIND=8),DIMENSION(:,:,:,:,:,:),ALLOCATABLE            :: y_out,zhat_pol_rtp

ALLOCATE(zhat_pol_rtp(natom,3,2,3,3,framecount_rtp))
ALLOCATE(y_out(natom,3,2,3,3,framecount_rtp_pade)) 

zhat_pol_rtp=COMPLEX(0.d0,0.0d0)
omega=5.0D0
data1=0.0d0
data2=0.0d0
broad=0.0d0
zhat_pol_dq_rtp=0.0d0
factor=1.d0/(2.0d0*dx)
start_freq=1.0d0
end_freq=INT(MAXVAL(freq)+1000.0D0)
freq_range=INT(end_freq-start_freq)
!mass_inv_sqrt(:)=REAL(1.0d0/SQRT(mass_atom(:)),KIND=8)

ALLOCATE(data1(freq_range*nmodes),data2(freq_range*nmodes))

!!!FFT of the RTP polarizabilities
DO j=1,natom
 DO i=1,3
  DO k=1,2
   DO m=1,3
    DO o=1,3
     CALL dfftw_plan_dft_r2c_1d(plan,framecount_rtp,pol_rtp(j,i,k,m,o,1:framecount_rtp),&
     zhat_pol_rtp(j,i,k,m,o,1:framecount_rtp),FFTW_ESTIMATE)
     CALL dfftw_execute_dft_r2c(plan,pol_rtp(j,i,k,m,o,1:framecount_rtp),zhat_pol_rtp(j,i,k,m,o,1:framecount_rtp))
     CALL dfftw_destroy_plan(plan)
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO

IF (check_pade=='y') THEN
!!Call Pade
 !$OMP PARALLEL DO COLLAPSE(5)
 DO j=1,natom !!natom
  DO i=1,3 !!dims
   DO k=1,2 !!disp
    DO m=1,3
     DO o=1,3
      CALL interpolate(framecount_rtp,zhat_pol_rtp(j,i,k,m,o,1:framecount_rtp),framecount_rtp_pade,y_out(j,i,k,m,o,:))
     ENDDO
    ENDDO
   ENDDO 
  ENDDO 
 ENDDO
 !$OMP END PARALLEL DO
framecount_rtp=framecount_rtp_pade
ENDIF

ALLOCATE(zhat_pol_dxyz_rtp(natom,3,3,3,framecount_rtp))
ALLOCATE(zhat_pol_dq_rtp(nmodes,3,3,framecount_rtp))
ALLOCATE(iso_sq(nmodes,framecount_rtp),aniso_sq(nmodes,framecount_rtp))
ALLOCATE(raman_int(nmodes,framecount_rtp))

!!!Finding laser frequency
rtp_freq_range=REAL(dom_rtp/framecount_rtp,KIND=8)
rtp_point=ANINT(laser_in_resraman/rtp_freq_range,KIND=8)
print*,laser_in_resraman,"laser_in_Resraman",rtp_freq_range,"rtp_freq_range",dom_rtp,"dom_rtp",rtp_point,'rtp_point'

!!!Finite differences
DO j=1,natom
 DO i=1,3
  DO m=1,3
   DO o=1,3
    DO n=1,framecount_rtp 
     IF (check_pade=='n') THEN
      zhat_pol_dxyz_rtp(j,i,m,o,n)=REAL((REAL(zhat_pol_rtp(j,i,2,m,o,n),KIND=8)-REAL(zhat_pol_rtp(j,i,1,m,o,n),&
              KIND=8))*factor,KIND=8)
     ELSEIF (check_pade=='y') THEN
      zhat_pol_dxyz_rtp(j,i,m,o,n)=REAL((REAL(y_out(j,i,2,m,o,n),KIND=8)-REAL(y_out(j,i,1,m,o,n),KIND=8))*factor,KIND=8)
     ENDIF
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO

!!!Derivatives w.r.t. mass weighted normal coordinates
DO i=1,nmodes
 DO j=1,natom
  DO o=1,framecount_rtp
   zhat_pol_dq_rtp(i,:,:,o)=zhat_pol_dq_rtp(i,:,:,o)+(zhat_pol_dxyz_rtp(j,1,:,:,o)*disp(j,i,1)*atom_mass_inv_sqrt(j))&
   +(zhat_pol_dxyz_rtp(j,2,:,:,o)*disp(j,i,2)*atom_mass_inv_sqrt(j))+(zhat_pol_dxyz_rtp(j,3,:,:,o)*&
   disp(j,i,3)*atom_mass_inv_sqrt(j))
  ENDDO
 ENDDO
ENDDO

!!!Isotropic and anisotropic contributions!!
DO k=1,nmodes
 DO o=1,framecount_rtp
  iso_sq(k,o)=REAL((zhat_pol_dq_rtp(k,1,1,o)+zhat_pol_dq_rtp(k,2,2,o)+zhat_pol_dq_rtp(k,3,3,o))/3.0d0,KIND=8)**2.0d0
        
  aniso_sq(k,o)=(0.50d0*(((zhat_pol_dq_rtp(k,1,1,o)-zhat_pol_dq_rtp(k,2,2,o))**2.0d0)+&
  ((zhat_pol_dq_rtp(k,2,2,o)-zhat_pol_dq_rtp(k,3,3,o))**2.0d0)+((zhat_pol_dq_rtp(k,3,3,o)-&
  zhat_pol_dq_rtp(k,1,1,o))**2.0d0)))+(3.0d0*((zhat_pol_dq_rtp(k,1,2,o)**2.0d0)+&
  (zhat_pol_dq_rtp(k,2,3,o)**2.0d0)+(zhat_pol_dq_rtp(k,3,1,o)**2.0d0)))
 ENDDO
ENDDO

!!!Calculation of the intensities!!
DO k=1,nmodes
 raman_int(k,rtp_point)=REAL(((7.0d0*aniso_sq(k,rtp_point))+(45.0d0*iso_sq(k,rtp_point)))/45.0d0,KIND=8)*&
 REAL(((laser_in_resraman-freq(k))**4.0d0)/freq(k),KIND=8)*REAL(1.0d0/(1.0d0-EXP(-1.438777d0*freq(k)/temp)),KIND=8)
ENDDO

!!!Broadening the spectrum!!
DO x=start_freq,end_freq
 DO i=1,nmodes
  broad=broad+(raman_int(i,rtp_point)*(1.0d0/(omega*SQRT(2.0d0*pi)))*EXP(-0.50d0*((x-freq(i))/omega)**2.0d0))
 ENDDO
 data2(x)=data2(x)+broad
ENDDO
       
OPEN(UNIT=98,FILE='result_static_raman.txt',STATUS='unknown',IOSTAT=stat)
DO i=start_freq,end_freq
 WRITE(98,*) i,data2(i)
ENDDO
CLOSE(98)

DEALLOCATE(iso_sq,aniso_sq,data1,data2,pol_rtp)
DEALLOCATE(zhat_pol_rtp,zhat_pol_dxyz_rtp,zhat_pol_dq_rtp)!,mass_inv_sqrt)
DEALLOCATE(raman_int,y_out)

END SUBROUTINE spec_static_resraman

!....................................................................................................................!
!....................................................................................................................!

SUBROUTINE spec_resraman(natom,framecount,element,rtp_dipole_x,rtp_dipole_y,rtp_dipole_z,type_input,mol_num,system,&
                read_function,dt,t_cor,pi,z_iso_resraman,z_aniso_resraman,dom,speed_light,const_planck,const_boltz,&
                const_permit,temp,dom_rtp,laser_in_resraman,y_out)

CHARACTER(LEN=40),INTENT(INOUT)                          :: read_function,system
CHARACTER(LEN=40),INTENT(INOUT)                          :: type_input
CHARACTER(LEN=40),INTENT(INOUT)                          :: rtp_dipole_x,rtp_dipole_y,rtp_dipole_z
CHARACTER(LEN=2),DIMENSION(:),ALLOCATABLE,INTENT(INOUT)  :: element
INTEGER,INTENT(INOUT)                                    :: natom,framecount,mol_num,t_cor
REAL(KIND=8),INTENT(INOUT)                               :: dt,pi,dom,dom_rtp
REAL(KIND=8),INTENT(INOUT)                               :: speed_light,laser_in_resraman
REAL(KIND=8),INTENT(INOUT)                               :: const_planck,const_boltz,const_permit,temp
COMPLEX(KIND=8),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT) :: z_iso_resraman,z_aniso_resraman
COMPLEX(KIND=8),DIMENSION(:),ALLOCATABLE,INTENT(INOUT)   :: y_out

CHARACTER(LEN=40)                                        :: chara
INTEGER                                                  :: stat,i,j,k,m,t0,t1
INTEGER(KIND=8)                                          :: plan
COMPLEX(KIND=8),DIMENSION(:,:,:),ALLOCATABLE             :: yx_out,yy_out,yz_out
COMPLEX(KIND=8),DIMENSION(:,:),ALLOCATABLE               :: zhat_iso_resraman,zhat_aniso_resraman
COMPLEX(KIND=8),DIMENSION(:,:,:),ALLOCATABLE             :: zhat_resraman_x,zhat_resraman_y,zhat_resraman_z
REAL(KIND=8)                                             :: f,freq_range,rtp_freq_range,pade_freq_range,laser_in
REAL(KIND=8),DIMENSION(:),ALLOCATABLE                    :: trace,abs_intens,trace_pade,abs_intens_pade
REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE                  :: zhat_unpol_resraman
REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE                :: alpha_resraman_x,alpha_resraman_y,alpha_resraman_z 
REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE                :: alpha_resraman_x_diff_re,alpha_resraman_y_diff_re
REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE                :: alpha_resraman_z_diff_re
REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE                :: alpha_resraman_x_diff_im,alpha_resraman_y_diff_im
REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE                :: alpha_resraman_z_diff_im
REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE                :: alpha_resraman_x_im,alpha_resraman_y_im
REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE                :: alpha_resraman_z_im
REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE                :: alpha_resraman_x_re,alpha_resraman_y_re
REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE                :: alpha_resraman_z_re
REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE                :: alpha_x,alpha_y,alpha_z
REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE                :: dip_x,dip_y,dip_z
      
ALLOCATE(zhat_resraman_x(framecount,natom,3),zhat_resraman_y(framecount,natom,3),zhat_resraman_z(framecount,natom,3))
ALLOCATE(alpha_resraman_x(framecount,natom,3),alpha_resraman_y(framecount,natom,3),alpha_resraman_z(framecount,natom,3))
ALLOCATE(alpha_resraman_x_diff_re(framecount-2,natom,3),alpha_resraman_y_diff_re(framecount-2,natom,3),&
                alpha_resraman_z_diff_re(framecount-2,natom,3))
ALLOCATE(alpha_resraman_x_diff_im(framecount-2,natom,3),alpha_resraman_y_diff_im(framecount-2,natom,3),&
                alpha_resraman_z_diff_im(framecount-2,natom,3))
ALLOCATE(alpha_x(framecount,natom,3),alpha_y(framecount,natom,3),alpha_z(framecount,natom,3))
ALLOCATE(yx_out(framecount,10000,3),yy_out(framecount,10000,3),yz_out(framecount,10000,3)) 

!!X-Field!!
CALL read_coord_frame(mol_num,framecount,element,rtp_dipole_x,dip_x)
CALL forward_diff(mol_num,framecount,alpha_x,dip_x,dip_x,system,read_function)
        
zhat_resraman_x = COMPLEX(0.d0,0.0d0)
zhat_resraman_y = COMPLEX(0.d0,0.0d0)
zhat_resraman_z = COMPLEX(0.d0,0.0d0)
        
DO i=1,framecount
  DO j=1,3
    CALL dfftw_plan_dft_r2c_1d(plan,natom,alpha_x(i,1:natom,j),zhat_resraman_x(i,1:natom,j),FFTW_ESTIMATE)
    CALL dfftw_execute_dft_r2c(plan,alpha_x(i,1:natom,j),zhat_resraman_x(i,1:natom,j)) !!!important to specify arrays!!
    CALL dfftw_destroy_plan(plan)
  ENDDO
ENDDO

alpha_resraman_x_re=REAL(zhat_resraman_x,KIND=8)
alpha_resraman_x_im=AIMAG(zhat_resraman_x)

!!Call Pade
DO i=1,1 !!framecount
 DO j=1,3
  CALL interpolate(natom-1,zhat_resraman_x(i,1:natom,j), 10000, yx_out(i,:,j))
 ENDDO
ENDDO

!OPEN(UNIT=40,FILE='y_out.txt',STATUS='unknown',IOSTAT=stat)
!DO i=1,10000
 ! WRITE(40,*) i/10000.d0,REAL(yx_out(i),KIND=8),AIMAG(yx_out(i))
!ENDDO
!CLOSE(40)

!OPEN(UNIT=40,FILE='zhat_out.txt',STATUS='unknown',IOSTAT=stat)
!DO i=1,natom-1
 ! WRITE(40,*) i/(natom-1.d0),REAL(zhat_resraman_x(1,i,1),KIND=8),AIMAG(zhat_resraman_x(1,i,1))
!ENDDO
!CLOSE(40)

CALL central_diff(dt,natom,framecount,alpha_resraman_x_re,alpha_resraman_x_diff_re,read_function,mol_num,system)
CALL central_diff(dt,natom,framecount,alpha_resraman_x_im,alpha_resraman_x_diff_im,read_function,mol_num,system)

!!Y-Field!!

CALL read_coord_frame(natom,framecount,element,rtp_dipole_y,dip_y)
CALL forward_diff(mol_num,framecount,alpha_y,dip_y,dip_y,system,read_function)

DO i=1,framecount
 DO j=1,3
  CALL dfftw_plan_dft_r2c_1d(plan,natom,alpha_y(i,1:natom,j),zhat_resraman_y(i,1:natom,j),FFTW_ESTIMATE)
  CALL dfftw_execute_dft_r2c(plan,alpha_y(i,1:natom,j),zhat_resraman_y(i,1:natom,j)) !!!important to specify arrays!!
  CALL dfftw_destroy_plan(plan)
 ENDDO
ENDDO

alpha_resraman_y_re=REAL(zhat_resraman_x,KIND=8)
alpha_resraman_y_im=AIMAG(zhat_resraman_x)

!!Call Pade
DO i=1,1 !!framecount
 DO j=1,3
  CALL interpolate(natom,zhat_resraman_y(i,1:natom,j), 10000, yy_out(i,:,j))
 ENDDO
ENDDO

!OPEN(UNIT=40,FILE='yy_out.txt',STATUS='unknown',IOSTAT=stat)
!DO i=1,10000
!  WRITE(40,*) i/10000.d0,REAL(yy_out(i),KIND=8),AIMAG(yy_out(i))
!ENDDO
!CLOSE(40)

!pade_interpolation.f90OPEN(UNIT=40,FILE='zhaty_out.txt',STATUS='unknown',IOSTAT=stat)
!DO i=1,natom-1
!  WRITE(40,*) i/(natom-1.d0),REAL(zhat_resraman_y(1,i,2),KIND=8),AIMAG(zhat_resraman_y(1,i,2))
!ENDDO
!CLOSE(40)

CALL central_diff(dt,natom,framecount,alpha_resraman_y_re,alpha_resraman_y_diff_re,read_function,mol_num,system)
CALL central_diff(dt,natom,framecount,alpha_resraman_y_im,alpha_resraman_y_diff_im,read_function,mol_num,system)

!!Z-Field!!

CALL read_coord_frame(natom,framecount,element,rtp_dipole_z,dip_z)
CALL forward_diff(mol_num,framecount,alpha_z,dip_z,dip_z,system,read_function)

DO i=1,framecount
 DO j=1,3
  CALL dfftw_plan_dft_r2c_1d(plan,natom,alpha_z(i,1:natom-1,j),zhat_resraman_z(i,1:natom,j),FFTW_ESTIMATE)
  CALL dfftw_execute_dft_r2c(plan,alpha_z(i,1:natom,j),zhat_resraman_z(i,1:natom,j)) !!!important to specify arrays!!
  CALL dfftw_destroy_plan(plan)
 ENDDO
ENDDO
       
alpha_resraman_z_re=REAL(zhat_resraman_x,KIND=8)
alpha_resraman_z_im=AIMAG(zhat_resraman_x)

!!Call Pade
DO i=1,1 !!framecount
 DO j=1,3
  CALL interpolate(natom,zhat_resraman_z(i,1:natom,j), 10000, yz_out(i,:,j))
 ENDDO
ENDDO

!OPEN(UNIT=40,FILE='yz_out.txt',STATUS='unknown',IOSTAT=stat)
!DO i=1,10000
!  WRITE(40,*) i/10000.d0,REAL(yz_out(i),KIND=8),AIMAG(yz_out(i))
!ENDDO
!CLOSE(40)

!OPEN(UNIT=40,FILE='zhatz_out.txt',STATUS='unknown',IOSTAT=stat)
!DO i=1,natom-1
!  WRITE(40,*) i/(natom-1.d0),REAL(zhat_resraman_z(1,i,3),KIND=8),AIMAG(zhat_resraman_z(1,i,3))
!ENDDO
!CLOSE(40)

CALL central_diff(dt,natom,framecount,alpha_resraman_z_re,alpha_resraman_z_diff_re,read_function,mol_num,system)
CALL central_diff(dt,natom,framecount,alpha_resraman_z_im,alpha_resraman_z_diff_im,read_function,mol_num,system)

!!!Calculate absorption spectra

rtp_freq_range=REAL(dom_rtp/(natom),KIND=8)
pade_freq_range=REAL(dom_rtp/(10000),KIND=8)

ALLOCATE(abs_intens(natom),trace(natom))
ALLOCATE(abs_intens_pade(10000),trace_pade(10000))

DO i=1,1
 DO j=1,natom
   trace(j)=DIMAG(zhat_resraman_x(i,j,1))+DIMAG(zhat_resraman_y(i,j,2))+DIMAG(zhat_resraman_z(i,j,3))
   abs_intens(j)=(4.0d0*pi*trace(j))/(3.0d0*speed_light)
 ENDDO
ENDDO

OPEN(UNIT=41,FILE='absorption_spectra.txt',STATUS='unknown',IOSTAT=stat)
DO i=1,natom
  WRITE(41,*) i*rtp_freq_range*1.23984198E-4,abs_intens(i)*i*rtp_freq_range
ENDDO
CLOSE(41)

DO i=1,1
 DO j=1,10000
   trace_pade(j)=DIMAG(yx_out(i,j,1))+DIMAG(yy_out(i,j,2))+DIMAG(yz_out(i,j,3))
   abs_intens_pade(j)=(4.0d0*pi*trace_pade(j))/(3.0d0*speed_light)
 ENDDO
ENDDO

OPEN(UNIT=42,FILE='absorption_spectra_pade.txt',STATUS='unknown',IOSTAT=stat)
DO i=1,10000
  WRITE(42,*) i*pade_freq_range*1.23984198E-4,abs_intens_pade(i)*i*pade_freq_range
ENDDO
CLOSE(42)

!!Generate the spectrum!!

CALL cvv_resraman(framecount,natom,t_cor,dt,pi,alpha_resraman_x_diff_re,alpha_resraman_y_diff_re,&
          alpha_resraman_z_diff_re,alpha_resraman_x_diff_im,alpha_resraman_y_diff_im,alpha_resraman_z_diff_im,&
          z_iso_resraman,z_aniso_resraman)

ALLOCATE(zhat_iso_resraman(0:t_cor*2,natom),zhat_aniso_resraman(0:t_cor*2,natom))
ALLOCATE(zhat_unpol_resraman(0:t_cor*2,natom))

zhat_iso_resraman = COMPLEX(0.d0,0.0d0)
zhat_aniso_resraman = COMPLEX(0.d0,0.0d0)

DO j=1,natom
  CALL dfftw_plan_dft_1d(plan,2*t_cor,z_iso_resraman(0:t_cor*2,j),zhat_iso_resraman(0:t_cor*2,j),FFTW_FORWARD,FFTW_ESTIMATE)
  CALL dfftw_execute_dft(plan,z_iso_resraman(0:t_cor*2,j),zhat_iso_resraman(0:t_cor*2,j)) !!!important to specify arrays!!
  CALL dfftw_destroy_plan(plan)

  CALL dfftw_plan_dft_1d(plan,2*t_cor,z_aniso_resraman(0:t_cor*2,j),zhat_aniso_resraman(0:t_cor*2,j),FFTW_FORWARD,FFTW_ESTIMATE)
  CALL dfftw_execute_dft(plan,z_aniso_resraman(0:t_cor*2,j),zhat_aniso_resraman(0:t_cor*2,j)) !!!important to specify arrays!!
  CALL dfftw_destroy_plan(plan)
ENDDO

freq_range=REAL(dom/(2*t_cor),KIND=8)
j=ANINT(laser_in_resraman/rtp_freq_range,KIND=8)

!!!!UNPOLARIZED!!!!

!zhat_iso_resraman=AIMAG(zhat_iso_resraman) 
!zhat_aniso_resraman=AIMAG(zhat_aniso_resraman) 

!OPEN(UNIT=30,FILE='zhat_aimag.txt',STATUS='unknown',IOSTAT=stat)
!DO i=0,2*t_cor-2
! WRITE(30,*) REAL(zhat_aniso_resraman(i,j),KIND=8),AIMAG(zhat_aniso_resraman(i,j))
!ENDDO
!CLOSE(30)

f=freq_range*dt*1.883652d-4
OPEN(UNIT=73,FILE='o-NP_resraman.txt',STATUS='unknown',IOSTAT=stat)
DO i=0,2*t_cor-2
  ! j=22
   !zhat_iso_resraman(i+1,j),AIMAG(zhat_iso_resraman(i+1,j),KIND=8)*(f*(i+1)/SIN(f*(i+1)))**2.d0
   zhat_iso_resraman(i+1,j)=(zhat_iso_resraman(i+1,j))*(f*(i+1)/SIN(f*(i+1)))**2.d0
   zhat_aniso_resraman(i+1,j)=(zhat_aniso_resraman(i+1,j))*(f*(i+1)/SIN(f*(i+1)))**2.d0
        
   zhat_unpol_resraman(i,j)=(zhat_iso_resraman(i,j))+(zhat_aniso_resraman(i,j)*7.0d0/45.0d0)*&
           ((const_planck)/(8.0d0*const_boltz*const_permit*const_permit)*1d-29*0.421d0*dt*((((laser_in*j)-((i)*freq_range))**4)/&
           ((i)*freq_range))*(1.0d0/(1.0d0-EXP((-1.438777d0*((i)*freq_range))/temp))))*2.0d0*2.0d0*pi
    zhat_unpol_resraman(0,j)=0.0d0
    IF ((i*freq_range).GE.5000.0d0) CYCLE
    !WRITE(73,*) i*freq_range,REAL(zhat_unpol_resraman(i,j),KIND=8),j
    WRITE(73,*) i*freq_range,zhat_unpol_resraman(i,j),j
! ENDDO
ENDDO                
CLOSE(73)

END SUBROUTINE spec_resraman
END MODULE calc_spectra

