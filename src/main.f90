PROGRAM vib2d

USE,INTRINSIC           :: iso_c_binding
USE setup,              ONLY: constants,read_input,masses_charges,conversion,pbc_orthorombic,pbc_hexagonal      
USE read_traj,          ONLY: read_coord,read_coord_frame,read_static,read_static_resraman
USE dipole_calc,        ONLY: center_mass,wannier,wannier_frag,solv_frag_index
USE vel_cor,            ONLY: cvv,cvv_iso,cvv_aniso,cvv_only_x,cvv_resraman
USE fin_diff,           ONLY: central_diff,forward_diff,finite_diff_static,finite_diff_static_resraman
USE calc_spectra,       ONLY: spec_power,spec_ir,spec_raman,normal_mode_analysis,spec_static_raman,spec_abs,&
                              spec_static_resraman,spec_resraman
USE config_info,        ONLY: output_config_info


IMPLICIT NONE

INCLUDE 'fftw3.f03'

INTEGER                                         :: b,i,j,k,natom,framecount,framecount_rtp_pade,t0,t1
INTEGER                                         :: frm,t_cor,nu,tau,stat,mol_num,nmodes,framecount_rtp,nfrag
INTEGER                                         :: count_0, count_1, count_rate, count_max
INTEGER(KIND=8)                                 :: plan
INTEGER,DIMENSION(:),ALLOCATABLE                :: natom_frag,natom_frag_x,natom_frag_free,nfrag_BO,nfrag_BC,nfrag_Ph
INTEGER,DIMENSION(:,:,:),ALLOCATABLE            :: fragment
CHARACTER(LEN=40)                               :: system,filename,static_pol,read_function,length,type_input,output_dip
CHARACTER(LEN=40)                               :: filename_dip,type_dipole,direction,averaging,cell_type,coord_file
CHARACTER(LEN=40)                               :: normal_freq_file,normal_displ_file,frag_type,type_static,force_file
CHARACTER(LEN=40)                               :: static_dip_free_file,static_dip_x_file,static_dip_y_file,static_dip_z_file
CHARACTER(LEN=40)                               :: rtp_dipole_x,rtp_dipole_y,rtp_dipole_z,check_pade,charac 
CHARACTER(LEN=40)                               :: output_dip_free,output_dip_x,output_dip_y,output_dip_z
CHARACTER(LEN=40)                               :: dipole_free,dipole_x,dipole_y,dipole_z,output_findif_dip
CHARACTER(LEN=40)                               :: output_findif_x,output_findif_y,output_findif_z,input_mass
CHARACTER(LEN=40)                               :: wannier_free,wannier_x,wannier_y,wannier_z,periodic
CHARACTER(LEN=2),DIMENSION(:),ALLOCATABLE       :: element
REAL(KIND=8)                                    :: dist,box_all,box_x,box_y,box_z,debye,mass_tot,fs2s,reccm2ev
REAL(KIND=8)                                    :: freq_range,dom,ce,co,h_kbT,raman_eq,a,dt_rtp,dom_rtp
REAL(KIND=8)                                    :: const_planck,const_permit,speed_light,const_boltz,temp,laser_in
REAL(KIND=8)                                    :: f,tmax,omega,theta,sinth,costh,sinsq,const_charge,laser_in_resraman
REAL(KIND=8)                                    :: cossq,thsq,thcub,alpha,beta,gamma0,dt,pi,multiplier,dx,bohr2ang
REAL(KIND=8)                                    :: time_init, time_final, elapsed_time
REAL(KIND=8)                                    :: hessian_factor,ang,at_u,sinc_const,mass_tot_cell
REAL(KIND=8)                                    :: damping_constant,joule_unit,ev_unit,action_unit,hartreebohr2evang
REAL(KIND=8),DIMENSION(3)                       :: vec,vec_pbc,coord2,coord1
REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE       :: refpoint,refpoint_free,refpoint_x,refpoint_y,refpoint_z
REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE       :: alpha_resraman_x,alpha_resraman_y,alpha_resraman_z
REAL(KIND=8),DIMENSION(:),ALLOCATABLE           :: kissfft,z,norm,mass_atom,z_aniso,z_iso,z_ortho,z_para,zhat_depol
REAL(KIND=8),DIMENSION(:),ALLOCATABLE           :: atom_mass_inv_sqrt,charge
REAL(KIND=8),DIMENSION(:),ALLOCATABLE           :: zhat_para_all,zhat_depol_x,zhat_unpol_x,freq,raman_int,test_x
REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE         :: test_x2,mass_mat
COMPLEX(KIND=8),DIMENSION(:),ALLOCATABLE        :: zhat,zhat_aniso,zhat_iso,zhat_para,zhat_ortho,zhat_unpol,integral,zhat_test
COMPLEX(KIND=8),DIMENSION(:,:,:),ALLOCATABLE    :: zhat_resraman_x
!COMPLEX(KIND=8),DIMENSION(:),ALLOCATABLE        :: zhat_resraman_x
COMPLEX(KIND=8),DIMENSION(:,:),ALLOCATABLE      :: zhat_test2
COMPLEX(KIND=8),DIMENSION(:,:),ALLOCATABLE      :: z_iso_resraman,z_aniso_resraman
COMPLEX(KIND=8),DIMENSION(:),ALLOCATABLE        :: y_out
REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE         :: coord,mass,dipole2,refpoint2,mass_tot_frag
REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE       :: dip,dip_free,dip_x,dip_y,dip_z
REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE       :: disp,pol_dq
REAL(KIND=8),DIMENSION(:,:,:,:),ALLOCATABLE     :: com,pol_dq_rtp
REAL(KIND=8),DIMENSION(:,:,:,:),ALLOCATABLE     :: static_dipole_x,static_dipole_y,static_dipole_z,static_dipole_free
REAL(KIND=8),DIMENSION(:,:,:,:,:),ALLOCATABLE   :: pol,force
REAL(KIND=8),DIMENSION(:,:,:,:,:),ALLOCATABLE       :: static_dipole_x_rtp,static_dipole_y_rtp,static_dipole_z_rtp
REAL(KIND=8),DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: pol_rtp
REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE       :: coord_v, coord_dip,coord_f,coord_x,coord_y,coord_z,dipole,coord_v_x
REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE       :: coord_v_free,alpha_x,alpha_y,alpha_z,v
REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE       :: alpha_diff_x,alpha_diff_y,alpha_diff_z

CALL SYSTEM_CLOCK(count_0, count_rate, count_max) !Starting time
time_init = count_0*1.0d0/count_rate

CALL output_config_info()

CALL constants(const_charge,debye,t_cor,const_planck,const_permit,speed_light,const_boltz,temp,pi,dx,bohr2ang,fs2s,&
     damping_constant,joule_unit,ev_unit,action_unit,hartreebohr2evang,hessian_factor,at_u,ang,framecount_rtp_pade,reccm2ev)

CALL read_input(filename,static_pol,static_dip_free_file,static_dip_x_file,static_dip_y_file,static_dip_z_file,normal_freq_file,&
     normal_displ_file,read_function,system,length,box_all,box_x,box_y,box_z,dt,type_input,wannier_free,wannier_x,wannier_y,&
     wannier_z,input_mass,periodic,direction,averaging,type_dipole,cell_type,rtp_dipole_x,rtp_dipole_y,rtp_dipole_z,&
     framecount_rtp,dt_rtp,laser_in_resraman,frag_type,type_static,force_file,laser_in,check_pade)

CALL conversion(dt,dom,dt_rtp,dom_rtp,speed_light,freq_range,t_cor,sinc_const)

IF (read_function=='P') THEN
    CALL read_coord(natom,framecount,element,coord,filename,periodic,mol_num,system,read_function,&
            framecount_rtp,type_dipole)
    CALL masses_charges(natom,mass_atom,atom_mass_inv_sqrt,mass_mat,element,mass_tot,charge)
    CALL spec_power(z,zhat,type_input,freq_range,natom,framecount,t_cor,dt,element,filename,coord_v,v,input_mass,dom,&
         mass_atom,read_function,mol_num,mass_tot,pi,coord,system,frag_type)
    DEALLOCATE(element,z,zhat,coord_v,mass_atom,coord)

ELSEIF (read_function=='MD-IR') THEN             
        IF (system=='1' .OR. system=='2' .AND. type_dipole=='1') THEN !!fragment approach or whole supercell
        IF (cell_type=='1' .OR. cell_type=='2') THEN !!KP or SC
            CALL read_coord(natom,framecount,element,coord,filename,periodic,mol_num,system,read_function,&
                    framecount_rtp,type_dipole)
            CALL masses_charges(natom,mass_atom,atom_mass_inv_sqrt,mass_mat,element,mass_tot,charge)
            CALL read_coord_frame(natom,framecount,element,filename,coord_v)
              !  CALL wannier(element,filename,natom,box_all,box_x,box_y,box_z,vec,vec_pbc,debye,mol_num,&
               !        periodic,mass_tot,framecount,mass_atom,coord_v,dip)
            !    CALL frag_index(natom,filename,element,coord_v,fragment,natom_frag,framecount)
            CALL center_mass(natom_frag,natom,refpoint,coord_v,filename,element,box_all,box_x,box_y,&
                 box_z,vec,vec_pbc,fragment,mass_atom,framecount,cell_type,mass_tot_frag,frag_type,mol_num,&
                 nfrag,type_dipole,system,mass_tot_cell)
            CALL wannier_frag(natom_frag,filename,natom,element,coord_v,box_all,box_x,box_y,box_z,vec,vec_pbc,dipole,&
                    refpoint,fragment,framecount,mass_tot_frag,mol_num,system,type_dipole,charge,mass_tot_cell)
            CALL spec_ir(z,zhat,freq_range,natom,framecount,t_cor,dt,element,filename,coord_v,v,input_mass,dom,pi,&
                 mol_num,box_all,box_x,box_y,box_z,vec,vec_pbc,debye,periodic,mass_tot,mass_atom,type_input,dip,&
                 read_function,coord,type_dipole,dipole,system,mass_tot_frag,sinc_const,nfrag,frag_type)
        ELSEIF (cell_type=='3') THEN !!SC with solvent
            CALL read_coord(natom,framecount,element,coord,filename,periodic,mol_num,system,&
                 read_function,framecount_rtp,type_dipole)
            CALL read_coord_frame(natom,framecount,element,filename,coord_v)
            CALL solv_frag_index(natom,coord_v,filename,element,box_all,vec,vec_pbc,&
                 box_x,box_y,box_z,mass_atom,framecount,cell_type,refpoint,natom_frag,fragment,mass_tot_frag)
            CALL wannier_frag(natom_frag,filename,natom,element,coord_v,box_all,box_x,box_y,box_z,vec,vec_pbc,dipole,&
                    refpoint,fragment,framecount,mass_tot_frag,mol_num,system,type_dipole,charge,mass_tot_cell)
            CALL spec_ir(z,zhat,freq_range,natom,framecount,t_cor,dt,element,filename,coord_v,v,input_mass,dom,pi,&
                 mol_num,box_all,box_x,box_y,box_z,vec,vec_pbc,debye,periodic,mass_tot,mass_atom,type_input,dip,&
                 read_function,coord,type_dipole,dipole,system,mass_tot_frag,sinc_const,nfrag,frag_type)
        ENDIF
    ELSEIF (system=='2') THEN !!molecular approach (Berry phase)
            CALL read_coord(natom,framecount,element,coord,filename,periodic,mol_num,system,read_function,&
                    framecount_rtp,type_dipole)
            CALL spec_ir(z,zhat,freq_range,natom,framecount,t_cor,dt,element,filename,coord_v,v,input_mass,dom,pi,&
                 mol_num,box_all,box_x,box_y,box_z,vec,vec_pbc,debye,periodic,mass_tot,mass_atom,type_input,dip,&
                 read_function,coord,type_dipole,dipole,system,mass_tot_frag,sinc_const,nfrag,frag_type)
    ENDIF

ELSEIF (read_function=='MD-R') THEN
    CALL read_coord(natom,framecount,element,coord,wannier_free,periodic,mol_num,system,read_function,framecount_rtp,type_dipole)
    CALL masses_charges(natom,mass_atom,atom_mass_inv_sqrt,mass_mat,element,mass_tot,charge)
    CALL spec_raman(natom,framecount,element,coord,wannier_free,wannier_x,wannier_y,wannier_z,mass_atom,mass_tot,periodic,&
         mol_num,dt,dom,t_cor,speed_light,coord_v,v,type_input,box_all,box_x,box_y,box_z,vec,vec_pbc,debye,read_function,&
         z_iso,z_aniso,z_ortho,z_para,const_planck,const_boltz,const_permit,temp,laser_in,pi,filename,averaging,direction,&
         type_dipole,system,natom_frag,fragment,refpoint,dipole,cell_type,mass_tot_frag,frag_type,nfrag,charge,mass_tot_cell)
    DEALLOCATE(coord,coord_v,mass_atom,element)

ELSEIF (read_function=='NMA') THEN
!    CALL read_coord(natom,framecount,element,coord,filename,periodic,mol_num,system,read_function,&
!         framecount_rtp,type_dipole)
!    CALL masses_charges(natom,mass_atom,atom_mass_inv_sqrt,mass_mat,element,mass_tot,charge)
!    CALL read_static(natom,element,normal_freq_file,normal_displ_file,static_pol,pol,freq,disp,&
 !        nmodes,static_dip_free_file,static_dip_x_file,static_dip_y_file,static_dip_z_file,type_dipole,&
 !        static_dipole_x,static_dipole_y,static_dipole_z,static_dipole_free,read_function,type_static,force_file,force)
 !   CALL normal_mode_analysis(natom,force,dx,hartreebohr2evang,hessian_factor,mass_mat)
        
ELSEIF (read_function=='R') THEN
    CALL read_coord(natom,framecount,element,coord,filename,periodic,mol_num,system,read_function,framecount_rtp,type_dipole)
    CALL read_static(natom,element,normal_freq_file,normal_displ_file,static_pol,pol,freq,disp,nmodes,static_dip_free_file,&
         static_dip_x_file,static_dip_y_file,static_dip_z_file,type_dipole,static_dipole_x,static_dipole_y,static_dipole_z,&
         static_dipole_free,read_function,type_static,force_file,force)
    CALL finite_diff_static(natom,nmodes,pol,pol_dq,disp,mass_atom,dx,bohr2ang,static_dipole_free,static_dipole_x,&
         static_dipole_y,static_dipole_z,type_dipole)
    CALL spec_static_raman(nmodes,pol_dq,laser_in,freq,temp,raman_int,pi,element,coord,disp,bohr2ang,natom)
    DEALLOCATE(freq,disp)
    DEALLOCATE(pol,pol_dq,raman_int)
    DEALLOCATE(element,coord,mass_atom)

ELSEIF (read_function=='MD-RR') THEN
    !CALL read_coord(natom,framecount,element,coord,rtp_dipole_x,periodic,mol_num,system,&
     !    read_function,framecount_rtp,type_dipole)
   ! CALL spec_resraman(natom,framecount,element,rtp_dipole_x,rtp_dipole_y,rtp_dipole_z,type_input,mol_num,system,&
    !     read_function,dt,t_cor,pi,z_iso_resraman,z_aniso_resraman,dom,speed_light,const_planck,const_boltz,&
     !    const_permit,temp,dom_rtp,laser_in_resraman,y_out)

ELSEIF (read_function=='ABS') THEN
        CALL read_coord(natom,framecount,element,coord,filename,periodic,mol_num,system,read_function,framecount_rtp,type_dipole)
        CALL read_static_resraman(natom,element,normal_freq_file,normal_displ_file,freq,disp,&
             nmodes,static_dip_x_file,static_dip_y_file,static_dip_z_file,framecount_rtp,&
             static_dipole_x_rtp,static_dipole_y_rtp,static_dipole_z_rtp)
        !CALL read_coord_frame(framecount_rtp+1,natom*6,charac,static_dip_x_file,static_dipole_x_rtp)
        !print*,static_dipole_x_rtp(2,3,:),"checkpoint2"
        CALL finite_diff_static_resraman(natom,pol_rtp,static_dipole_x_rtp,static_dipole_y_rtp,&
             static_dipole_z_rtp,framecount_rtp,speed_light,fs2s,damping_constant,joule_unit,ev_unit,action_unit,dt_rtp)
        CALL spec_abs(nmodes,natom,pol_rtp,freq,pi,framecount_rtp,speed_light,framecount_rtp_pade,reccm2ev,check_pade,&
             dom_rtp)

ELSEIF (read_function=='RR') THEN
        CALL read_coord(natom,framecount,element,coord,filename,periodic,mol_num,system,read_function,framecount_rtp,type_dipole)
        CALL masses_charges(natom,mass_atom,atom_mass_inv_sqrt,mass_mat,element,mass_tot,charge)
        CALL read_static_resraman(natom,element,normal_freq_file,normal_displ_file,freq,disp,&
             nmodes,static_dip_x_file,static_dip_y_file,static_dip_z_file,framecount_rtp,&
             static_dipole_x_rtp,static_dipole_y_rtp,static_dipole_z_rtp)
        CALL finite_diff_static_resraman(natom,pol_rtp,static_dipole_x_rtp,static_dipole_y_rtp,&
             static_dipole_z_rtp,framecount_rtp,speed_light,fs2s,damping_constant,joule_unit,ev_unit,action_unit,dt_rtp)
        CALL spec_static_resraman(nmodes,natom,pol_rtp,laser_in_resraman,freq,temp,pi,framecount_rtp,dom_rtp,dx,bohr2ang,&
           disp,speed_light,framecount_rtp_pade,check_pade,atom_mass_inv_sqrt)

!ELSEIF (read_function=='7') THEN
     !   CALL read_coord(natom,framecount,element,coord,filename,mass_atom,mass_tot,periodic,mol_num,system,&
     !           read_function,framecount_rtp)
    !    CALL conversion(framecount,dt,dom,t_cor,speed_light)
   !     CALL central_diff(dt,natom,framecount,coord_v,v,read_function,mol_num,system)
  !      DEALLOCATE(coord_v,v)
 !       DEALLOCATE(element,mass_atom)
ENDIF

CALL SYSTEM_CLOCK(count_1, count_rate, count_max) !Ending time
time_final = count_1*1.0d0/count_rate
elapsed_time = time_final - time_init !Elapsed time

WRITE(*,1003) INT(elapsed_time),elapsed_time-INT(elapsed_time) !Write elapsed time
1003 FORMAT('  Wall Clock = ',i0,F0.9)

END PROGRAM vib2d 
