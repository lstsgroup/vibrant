MODULE setup

IMPLICIT NONE

PRIVATE

PUBLIC :: constants,read_input,masses_charges,conversion,pbc_orthorombic,pbc_hexagonal

CONTAINS
        
SUBROUTINE constants(const_charge,debye,t_cor,const_planck,const_permit,speed_light,const_boltz,&
           temp,pi,dx,bohr2ang,fs2s,damping_constant,joule_unit,ev_unit,action_unit,hartreebohr2evang,&
           hessian_factor,at_u,ang,framecount_rtp_pade,reccm2ev)

INTEGER,INTENT(OUT)                                 :: t_cor,framecount_rtp_pade
REAL(KIND=8),INTENT(OUT)                            :: debye,const_planck,const_permit,pi,const_charge,fs2s
REAL(KIND=8),INTENT(OUT)                            :: speed_light,const_boltz,temp,dx,bohr2ang,reccm2ev
REAL(KIND=8),INTENT(OUT)                            :: damping_constant,joule_unit,ev_unit,action_unit
REAL(KIND=8),INTENT(OUT)                            :: hartreebohr2evang,hessian_factor,at_u,ang

dx=0.00265d0  !!for r-met
!dx=0.001d0  !!for COF-1... etc
const_charge=1.602176565E-19  
bohr2ang=0.5291772109d0 !!bohr two angstrom
hartreebohr2evang=51.42208619083232d0 !!hartree/bohr to eV/angstrom
damping_constant=0.10d0 !! ev
joule_unit=4.359744722E-18 !! J
ev_unit=27.211386d0 !! ev
action_unit=1.054571817E-34 !J.s
reccm2ev=0.000124d0 !cm^-1 to eV
pi=3.14159d0        
temp=300.0d0        !K
debye=0.393456d0
t_cor=1024
const_planck=6.62607015E-34 !m^2*kg/s or J.s
const_permit=8.8541878128E-12 !F*m^−1        
speed_light=2.9979246E+10  !cm/s
const_boltz=1.380649E-23 !m^2*kg*s^-2*K-1 or J/K
fs2s=1.0E-15
ang=1.0E-10  !!angstrom to cm^-1
at_u=1.6605390666E-27 !!atomic mass unit-kilogram relationship
hessian_factor=REAL(const_charge/(at_u*ang*ang),KIND=8)
framecount_rtp_pade=80000
!laser_in=9398.50d0  !cm^-1       
!laser_in=200000.0d0  !cm^-1       
!laser_in=15797.78d0 !cm^-1
!laser_in=18796.99 
END SUBROUTINE constants

!*************************************************************************************************
!*************************************************************************************************

SUBROUTINE read_input(filename,static_pol,static_dip_free_file,static_dip_x_file,static_dip_y_file,&
        static_dip_z_file,normal_freq_file,normal_displ_file,read_function,system,&
        length,box_all,box_x,box_y,box_z,dt,type_input,wannier_free,wannier_x,wannier_y,wannier_z,&
        input_mass,periodic,direction,averaging,type_dipole,cell_type,rtp_dipole_x,rtp_dipole_y,&
        rtp_dipole_z,framecount_rtp,dt_rtp,laser_in_resraman,frag_type,type_static,force_file,laser_in,check_pade)

CHARACTER(LEN=40),INTENT(OUT)                       :: filename,static_pol,read_function,length,system,type_input,periodic
CHARACTER(LEN=40),INTENT(OUT)                       :: wannier_free,wannier_x,wannier_y,wannier_z,input_mass,type_static
CHARACTER(LEN=40),INTENT(OUT)                       :: direction,averaging,type_dipole,cell_type,force_file
CHARACTER(LEN=40),INTENT(OUT)                       :: static_dip_free_file,static_dip_x_file,static_dip_y_file,static_dip_z_file
CHARACTER(LEN=40),INTENT(OUT)                       :: normal_freq_file,normal_displ_file,frag_type
CHARACTER(LEN=40),INTENT(OUT)                       :: rtp_dipole_x,rtp_dipole_y,rtp_dipole_z,check_pade
REAL(KIND=8),INTENT(OUT)                            :: dt,dt_rtp,box_all,box_x,box_y,box_z,laser_in,laser_in_resraman
INTEGER,INTENT(OUT)                                 :: framecount_rtp

laser_in=9398.50d0  !cm^-1       
DO 
 WRITE(*,*)'Enter which function you want to calculate (type "P" for Power spectrum , "MD-IR" for MD-based IR spectrum, &
         "MD-R" for MD-based Raman spectrum, "MD-RR" for MD-based resonance Raman, "NMA" for normal mode analysis, &
         "IR" for static IR spectrum, "R" for static Raman spectrum, "ABS" for absorption spectrum,"RR" for static &
         resonance Raman spectrum)'
 READ(*,*) read_function
 IF (read_function.NE.'P' .AND. read_function.NE.'MD-IR' .AND. read_function.NE.'MD-R' .AND. read_function.NE.'MD-RR' &
    .AND. read_function.NE.'NMA' .AND. read_function.NE.'IR' .AND. read_function.NE.'R' .AND. read_function.NE.'ABS' .AND. &
    read_function.NE.'RR') THEN
  WRITE(*,*) 'Please type P, MD-IR, MD-R, MD-RR, NMA, IR, R, ABS or RR!'
  CYCLE
 ENDIF
EXIT
ENDDO

!read_function='MD-R'

DO
    IF (read_function=='P') THEN
        WRITE(*,*) 'Enter the type of the input file (type 1 for positions, 2 for velocities)'
        READ(*,*) type_input 
        IF (type_input.NE.'1' .AND. type_input.NE.'2') THEN
            WRITE(*,*) 'Please type 1 or 2!'
            CYCLE
        ENDIF 
    ENDIF 
    EXIT
ENDDO

DO
    IF (read_function=='P') THEN
        WRITE(*,*) 'Do you want to apply mass weighting (y/n)?'
        READ(*,*) input_mass 
        IF (input_mass.NE.'y' .AND. input_mass.NE.'n') THEN
            WRITE(*,*) 'Please type y or n!'
            CYCLE
        ENDIF 
    ENDIF 
    EXIT
ENDDO

DO 
    IF (read_function=='IR' .OR. read_function=='R') THEN
        WRITE(*,*) 'Do you want the normal modes to be calculated (type "1") or read from an external file (type "2")?'
        READ(*,*) type_static
        IF (type_static.NE.'1' .AND. type_static.NE.'2') THEN
            WRITE(*,*) 'Please type 1 or 2!'
            CYCLE
        ENDIF
    ENDIF
    EXIT
ENDDO

DO
    IF (read_function=='MD-R' .OR. read_function=='R') THEN
        WRITE(*,*) 'Which one do you want to use: Wannier centers (1), Berry phase dipole moments (2)&
                    or DFPT polarizabilities (3)?'
        READ(*,*) type_dipole 
        IF (type_dipole.NE.'1' .AND. type_dipole.NE.'2' .AND. type_dipole.NE.'3') THEN
            WRITE(*,*) 'Please type 1, 2 or 3!!'
            CYCLE
        ENDIF 
    ENDIF 
    EXIT
ENDDO

DO
    IF (read_function=='MD-IR') THEN
        WRITE(*,*) 'Which one do you want to use: Wannier centers (1) or Berry phase dipole moments (2)?'
        READ(*,*) type_dipole 
        IF (type_dipole.NE.'1' .AND. type_dipole.NE.'2') THEN
            WRITE(*,*) 'Please type 1 or 2!!'
            CYCLE
        ENDIF 
    ENDIF 
    EXIT
ENDDO

DO
    IF (read_function=='MD-IR' .OR. read_function=='MD-R') THEN
        WRITE(*,*)'Do you want to apply the fragment approach (1) or the molecular approach? (2)'
        READ(*,*) system
        IF (system.NE.'1' .AND. system.NE.'2') THEN 
            WRITE(*,*) 'Please type 1 or 2!'
            CYCLE
        ENDIF
    ENDIF
    EXIT
ENDDO  

DO
    IF (read_function=='MD-IR' .OR. read_function=='MD-R') THEN
        WRITE(*,*)'Is it the k-point trajectory (1), supercell trajectory (2) or solvent trajectory (3)? '
        READ(*,*) cell_type
        IF (cell_type.NE.'1' .AND. cell_type.NE.'2' .AND. cell_type.NE.'3') THEN 
            WRITE(*,*) 'Please type 1, 2 or 3!'
            CYCLE
        ENDIF
    ENDIF
    EXIT
ENDDO  

DO
    IF (read_function=='MD-R') THEN
        IF (type_dipole=='3') THEN
            IF (cell_type=='1') THEN
                wannier_free='alpha_x_kp_5000.xyz'
                wannier_x='alpha_x_kp_5000.xyz'
                wannier_y='alpha_y_kp_5000.xyz'
                wannier_z='alpha_z_kp_5000.xyz'
            ELSEIF (cell_type=='2') THEN
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
                wannier_free='alpha_x_solv_5000_new.xyz'
                wannier_x='alpha_x_solv_5000_new.xyz'
                wannier_y='alpha_y_solv_5000_new.xyz'
                wannier_z='alpha_z_solv_5000_new.xyz'
            ENDIF
        ELSEIF (type_dipole=='2') THEN
            IF (cell_type=='1') THEN
                wannier_free='berry_dipole_free_5000_kp.xyz'
                wannier_x='berry_dipole_x_5000_kp_large.xyz'
                wannier_y='berry_dipole_y_5000_kp_large.xyz'
                wannier_z='berry_dipole_z_5000_kp_large.xyz'
            ELSEIF (cell_type=='2') THEN
              !  wannier_free='dipole_o-NP_free.xyz'
                ! wannier_x='dipole_o-NP_X.xyz'
                ! wannier_y='dipole_o-NP_Y.xyz'
                ! wannier_z='dipole_o-NP_Z.xyz'
              !  wannier_x='dipole_o-NP_X_smallfield.xyz'
              !  wannier_y='dipole_o-NP_Y_smallfield.xyz'
              !  wannier_z='dipole_o-NP_Z_smallfield.xyz'
                 wannier_free='berry_dipole_free_5000_sc.xyz'
                 wannier_x='berry_dipole_X_5000_sc.xyz' 
                 wannier_y='berry_dipole_Y_5000_sc.xyz'
                 wannier_z='berry_dipole_Z_5000_sc.xyz'
            ENDIF
        ELSEIF (type_dipole=='1') THEN
            IF (cell_type=='1') THEN
                wannier_free='wannier_free_COF-1_kp_5000.xyz'
                wannier_x='wannier_X_COF-1_kp_large_5000.xyz'
                wannier_y='wannier_Y_COF-1_kp_large_5000.xyz'
                wannier_z='wannier_Z_COF-1_kp_large_5000.xyz'
            ELSEIF (cell_type=='2') THEN
                wannier_free='wannier_free_COF-1_sc_5000.xyz'
                wannier_x='wannier_X_COF-1_sc_5000.xyz'
                wannier_y='wannier_Y_COF-1_sc_5000.xyz'
                wannier_z='wannier_Z_COF-1_sc_5000.xyz'
            ELSEIF (cell_type=='3') THEN
                wannier_free='COF-1_solv_wannier_free.xyz'
                wannier_x='COF-1_solv_wannier_X.xyz'
                wannier_y='COF-1_solv_wannier_Y.xyz'
                wannier_z='COF-1_solv_wannier_Z.xyz'
            ENDIF
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
        ENDIF 
    ENDIF 
    EXIT
ENDDO

DO
    IF (read_function=='MD-R') THEN
        WRITE(*,*) 'What is the type of averaging you want to apply? (type 1 for orientational, 2 for isotropic)'
        READ(*,*) averaging
        IF (averaging.NE.'1' .AND. averaging.NE.'2') THEN
            WRITE(*,*) 'Please type 1 or 2!'
            CYCLE
        ENDIF
    ENDIF
    EXIT
ENDDO

DO 
    IF (read_function=='MD-R' .AND. averaging=='2') THEN
        WRITE(*,*) 'What is the direction of the applied electric field (type 1 for x, 2 for y, 3 for z)?'
        READ(*,*) direction
        IF (direction.NE.'1' .AND. direction.NE.'2' .AND. direction.NE.'3') THEN
            WRITE(*,*) 'Please type 1, 2 or 3!'
            CYCLE
        ENDIF
    ENDIF
    EXIT
ENDDO

DO 
    IF (read_function=='MD-RR' .OR. read_function=='RR') THEN
        ! WRITE(*,*) 'What is the wavenumber (cm^-1) of the incident laser?'
        ! READ(*,*) laser_in
    ENDIF
    EXIT
ENDDO

DO
    IF (read_function=='MD-IR' .OR. read_function=='MD-R') THEN
        WRITE(*,*) 'Does the system contain more than one molecule? (y/n)'
        READ(*,*) periodic
        IF (periodic.NE.'y' .AND. periodic.NE.'n') THEN
            WRITE(*,*) 'Please type y or n!'
            CYCLE
        ENDIF 
    ENDIF
    EXIT
ENDDO


DO
    IF (system=='1') THEN
            WRITE(*,*)'Which fragments do you want to calculate? (BO (1), Ph (2) or B-C bonds (3)?'
        READ(*,*) frag_type
        IF (frag_type.NE.'1' .AND. frag_type.NE.'2' .AND. frag_type.NE.'3') THEN 
            WRITE(*,*) 'Please type 1, 2 or 3!'
            CYCLE
        ENDIF
    ENDIF
    EXIT
ENDDO  

DO 
    IF (read_function=='NMA' .OR. type_static=='1') THEN
        force_file='BDBA-force.data3'
        filename='BDBA.xyz'
    ENDIF
    EXIT
ENDDO

DO
    IF (read_function=='R' ) THEN
        ! normal_freq_file='normal_freqs_o-NP.dat'
        ! normal_displ_file='normal_displacements_o-NP.dat'
        ! normal_freq_file='normal_freqs_2cat_triplet.dat'
        ! normal_displ_file='normal_displacements_2cat_triplet.dat'
        normal_freq_file='normal_freqs_2cat_singlet_CS.dat'
        normal_displ_file='normal_displacements_2cat_singlet_CS.dat'
        !normal_freq_file='normal_freqs_r-met.dat'
        !normal_displ_file='normal_displacements_r-met.dat'
        ! filename='water.xyz'
        !filename='o-nitrophenol.xyz'
        !filename='2cat_triplet.xyz'
        filename='2cat_singlet_CS.xyz'
        !filename='r-met.xyz'
        !filename='COF-1.xyz'
        IF (type_dipole=='3') THEN
            static_pol='polarizabilities.dat'
        ELSEIF (type_dipole=='2') THEN
            static_dip_free_file='dipole_2cat_singlet_CS_free_static.xyz'
            static_dip_x_file='dipole_2cat_singlet_CS_X_static.xyz'
            static_dip_y_file='dipole_2cat_singlet_CS_Y_static.xyz'
            static_dip_z_file='dipole_2cat_singlet_CS_Z_static.xyz'
            !static_dip_free_file='dipole_r-met_free_static.xyz'
            ! static_dip_x_file='dipole_r-met_X_static.xyz'
            ! static_dip_y_file='dipole_r-met_Y_static.xyz'
            ! static_dip_z_file='dipole_r-met_Z_static.xyz'
            ! static_dip_free_file='dipole_o-NP_free_static.xyz'
            ! static_dip_x_file='dipole_o-NP_X_static.xyz'
            ! static_dip_y_file='dipole_o-NP_Y_static.xyz'
            ! static_dip_z_file='dipole_o-NP_Z_static.xyz'
         ENDIF
        !WRITE(*,*)'Enter the name of the normal frequencies file'
        !READ(*,*) norma_freq_file
        !WRITE(*,*)'Enter the name of the normal displacements file'
        !READ(*,*) norma_displ_file
        !WRITE(*,*)'Enter the name of the coordinate file'
        !READ(*,*) filename
        !WRITE(*,*)'Enter the name of the polarizability file'
        !READ(*,*) static_pol
    ENDIF
    EXIT
 ENDDO
 
DO
    IF (read_function=='MD-RR' ) THEN
        rtp_dipole_x='o-NP_RTP_dipoles_X_256.xyz'
        rtp_dipole_y='o-NP_RTP_dipoles_Y_256.xyz'
        rtp_dipole_z='o-NP_RTP_dipoles_Z_256.xyz'
        ! rtp_dipole_x='o-NP_RTP_dipoles_X.xyz'
        ! rtp_dipole_y='o-NP_RTP_dipoles_Y.xyz'
        ! rtp_dipole_z='o-NP_RTP_dipoles_Z.xyz'
        ! WRITE(*,*)'What is the number of RTP frames?'
        ! READ(*,*) framecount_rtp
        ! framecount_rtp=1280
        framecount_rtp=256
        ! dt_rtp=0.0125d0
        dt_rtp=0.0625d0
        WRITE(*,*)'What is the wavenumber of the incident laser (cm^-1)?'
        READ(*,*) laser_in_resraman
        !laser_in_resraman=15797.788309636651d0
    ENDIF
    EXIT
ENDDO

DO
    IF (read_function=='RR' .OR. read_function=='ABS') THEN
        ! normal_freq_file='normal_freqs_o-NP.dat'
        ! normal_displ_file='normal_displacements_o-NP.dat'
        normal_freq_file='normal_freqs_r-met.dat'
        normal_displ_file='normal_displacements_r-met.dat'
        ! filename='o-nitrophenol.xyz'
        filename='r-met.xyz'
        !  static_dip_x_file='o-NP_RTP_dipoles_static_X.xyz'
        !  static_dip_y_file='o-NP_RTP_dipoles_static_Y.xyz'
        !  static_dip_z_file='o-NP_RTP_dipoles_static_Z.xyz'
        static_dip_x_file='r-met_RTP_dipoles_static_X.xyz'
        static_dip_y_file='r-met_RTP_dipoles_static_Y.xyz'
        static_dip_z_file='r-met_RTP_dipoles_static_Z.xyz'
        !framecount_rtp=1280
        framecount_rtp=50000
        ! framecount_rtp=256
        !dt_rtp=0.0125d0
        dt_rtp=0.00242d0
        !  dt_rtp=0.0625d0
        !WRITE(*,*)'What is the wavenumber of the incident laser (cm^-1)?'
        !READ(*,*) laser_in_resraman
        laser_in_resraman=15797.788309636651d0
        !WRITE(*,*) 'Do you want to apply Pade approximants? (y/n)'
        !READ(*,*) check_pade
        check_pade='n'
    ENDIF
    EXIT
ENDDO

DO
    IF (read_function=='P' .OR. read_function=='MD-IR') THEN
        WRITE(*,*)'Enter the name of the trajectory'
        READ(*,*) filename
    ENDIF
    IF (read_function=='P' .OR. read_function=='MD-IR' .OR. read_function=='MD-R') THEN
        WRITE(*,*)'Enter the time step (fs)'
        READ(*,*) dt
        IF (read_function.NE.'P') THEN
            WRITE(*,*)'Are the 3 cell vectors the same length? (y/n)'
            READ(*,*) length
            IF (length=='y') THEN
                WRITE(*,*)'Enter the cell vector (in Armstrong)'
                READ(*,*) box_all
            ELSEIF (length=='n') THEN
                IF (cell_type=='1') THEN
                    box_x=15.100d0
                    box_y=15.100d0
                    box_z=13.448d0
                ELSEIF (cell_type=='2') THEN
                    box_x=15.100d0
                    box_y=15.101d0
                    box_z=13.457d0
                ELSEIF (cell_type=='3') THEN
                    box_x=15.100d0
                    box_y=15.100d0
                    box_z=20.172d0
                ENDIF    
               !         WRITE(*,*)'Enter the cell vector for X direction (in Armstrong)'
               !         READ(*,*) box_x
                !        WRITE(*,*)'Enter the cell vector for Y direction (in Armstrong)'
                 !       READ(*,*) box_y
                  !      WRITE(*,*)'Enter the cell vector for Z direction (in Armstrong)'
                   !     READ(*,*) box_z
            ENDIF
        ENDIF
    ENDIF    
    EXIT
ENDDO
 
IF (length=='y') THEN
    box_x=box_all
    box_y=box_all
    box_z=box_all
ENDIF        

END SUBROUTINE read_input

!*********************************************************************************************
!*********************************************************************************************

SUBROUTINE masses_charges(natom,mass_atom,atom_mass_inv_sqrt,mass_mat,element,mass_tot,charge)

CHARACTER(LEN=2),DIMENSION(:),ALLOCATABLE,INTENT(IN)   :: element
INTEGER,INTENT(INOUT)                                  :: natom
REAL(KIND=8),INTENT(INOUT)                             :: mass_tot
REAL(KIND=8),DIMENSION(:),ALLOCATABLE,INTENT(OUT)      :: atom_mass_inv_sqrt,mass_atom,charge
REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT)    :: mass_mat

INTEGER                                                :: i,j,stat
REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE                :: mat1,mat2

ALLOCATE(atom_mass_inv_sqrt(natom),mass_mat(natom,natom),mass_atom(natom),charge(natom))
ALLOCATE(mat1(natom,1),mat2(1,natom))

mass_atom=0.0d0
mass_tot=0.0d0
DO i=1,natom
    IF (element(i)=='O') THEN
        mass_atom(i)=15.999d0
        charge(i)=6.0d0
    ELSEIF (element(i)=='H') THEN
        mass_atom(i)=1.00784d0
        charge(i)=1.0d0
    ELSEIF (element(i)=='C') THEN
        mass_atom(i)=12.011d0
        charge(i)=4.0d0
    ELSEIF (element(i)=='B') THEN
        mass_atom(i)=10.811d0
        charge(i)=3.0d0
    ELSEIF (element(i)=='N') THEN
        mass_atom(i)=14.0067d0
        charge(i)=5.0d0
    ELSEIF (element(i)=='X') THEN
        mass_atom(i)=0.00d0
        charge(i)=-2.0d0
    ENDIF         
    mass_tot=mass_tot+mass_atom(i)
ENDDO

atom_mass_inv_sqrt(:)=SQRT(REAL(1.0d0/mass_atom(:),KIND=8))

mat1(:,:)=RESHAPE(atom_mass_inv_sqrt(:),(/natom,1/))
mat2(:,:)=RESHAPE(atom_mass_inv_sqrt(:),(/1,natom/))

mass_mat=MATMUL(mat1,mat2)

DEALLOCATE(mat1,mat2)

END SUBROUTINE masses_charges

!*********************************************************************************************
!*********************************************************************************************

SUBROUTINE conversion(dt,dom,dt_rtp,dom_rtp,speed_light,freq_range,t_cor,sinc_const)

INTEGER,INTENT(IN)                    :: t_cor
REAL(KIND=8),INTENT(OUT)              :: dom,dom_rtp,freq_range,sinc_const
REAL(KIND=8),INTENT(IN)               :: speed_light,dt,dt_rtp

INTEGER                               :: stat   ! error status of OPEN statements
INTEGER                               :: i, j, k

dom=REAL((1.0d0/(dt*1E-15))/speed_light,KIND=8)
dom_rtp=REAL((1.0d0/(dt_rtp*1E-15))/speed_light,KIND=8)

freq_range=REAL(dom/(2.0d0*t_cor),KIND=8)
sinc_const=freq_range*dt*1.883652d-4 !!for sinc function

END SUBROUTINE conversion

!*********************************************************************************************
!*********************************************************************************************

SUBROUTINE pbc_orthorombic(coord2,coord1,vec,vec_pbc,box_all,box_x,box_y,box_z)

REAL(KIND=8),DIMENSION(3),INTENT(INOUT)                :: vec,vec_pbc,coord2,coord1
REAL(KIND=8),INTENT(IN)                                :: box_all,box_x,box_y,box_z

vec(:)=coord2(:)-coord1(:)
     
vec_pbc(1)=vec(1)-box_x*ANINT((1./box_x)*vec(1))
vec_pbc(2)=vec(2)-box_y*ANINT((1./box_y)*vec(2))
vec_pbc(3)=vec(3)-box_z*ANINT((1./box_z)*vec(3))

END SUBROUTINE pbc_orthorombic

!********************************************************************************************
!********************************************************************************************

SUBROUTINE pbc_hexagonal(coord2,coord1,vec,vec_pbc,box_all,box_x,box_y,box_z)

REAL(KIND=8),DIMENSION(3),INTENT(INOUT)                :: vec,vec_pbc,coord2,coord1
REAL(KIND=8),INTENT(IN)                                :: box_all,box_x,box_y,box_z

REAL(KIND=8)                                           :: h_inv(3,3),a,s(3),hmat(3,3)
REAL(KIND=8)                                           :: acosa,asina,sqrt3,det_a

sqrt3=1.73205080756887729352744634d0

a = 0.5d0*(box_x + box_y)
acosa = 0.5d0*a
asina = sqrt3*acosa
hmat(1, 1) = a; hmat(1, 2) = acosa; hmat(1, 3) = 0.0d0
hmat(2, 1) = 0.0d0; hmat(2, 2) = asina; hmat(2, 3) = 0.0d0
hmat(3, 1) = 0.0d0; hmat(3, 2) = 0.0d0; hmat(3, 3) = box_z


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

vec(:)=coord2(:)-coord1(:)

s(1) = h_inv(1, 1)*vec(1) + h_inv(1, 2)*vec(2) + h_inv(1, 3)*vec(3)
s(2) = h_inv(2, 1)*vec(1) + h_inv(2, 2)*vec(2) + h_inv(2, 3)*vec(3)
s(3) = h_inv(3, 1)*vec(1) + h_inv(3, 2)*vec(2) + h_inv(3, 3)*vec(3)

s(1) = s(1) - ANINT(s(1))
s(2) = s(2) - ANINT(s(2))
s(3) = s(3) - ANINT(s(3))

vec_pbc(1) = hmat(1, 1)*s(1) + hmat(1, 2)*s(2) + hmat(1, 3)*s(3)
vec_pbc(2) = hmat(2, 1)*s(1) + hmat(2, 2)*s(2) + hmat(2, 3)*s(3)
vec_pbc(3) = hmat(3, 1)*s(1) + hmat(3, 2)*s(2) + hmat(3, 3)*s(3)

END SUBROUTINE pbc_hexagonal
END MODULE setup

