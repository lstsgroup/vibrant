MODULE read_traj

    USE kinds,              ONLY: dp

IMPLICIT NONE

PRIVATE 

PUBLIC :: read_coord,read_coord_frame,read_static,read_static_resraman

CONTAINS
SUBROUTINE read_coord(natom,framecount,element,coord,filename,periodic,mol_num,system,read_function,&
           framecount_rtp,type_dipole)

CHARACTER(LEN=40),INTENT(IN)                               :: filename,periodic,system,read_function,type_dipole
INTEGER,INTENT(IN)                                         :: framecount_rtp
INTEGER,INTENT(OUT)                                        :: natom,framecount,mol_num
CHARACTER(LEN=2),DIMENSION(:),ALLOCATABLE,INTENT(OUT)      :: element
REAL(kind=dp),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT)        :: coord

INTEGER                                                    :: i,j,stat

framecount=0

IF (read_function.NE.'MD-RR') THEN
    OPEN(UNIT=50,FILE=filename,STATUS='old',IOSTAT=stat)
    READ(50,*) natom
    CLOSE(50)
ELSEIF (read_function=='MD-RR') THEN
    natom=framecount_rtp
ENDIF        

ALLOCATE(element(natom),coord(natom,3))

OPEN(UNIT=51,FILE=filename,STATUS='old',IOSTAT=stat)
DO 
    READ(51,*,END=998)
    READ(51,*)
    framecount=framecount+1
    DO i=1,natom
        READ(51,*) element(i),coord(i,1),coord(i,2),coord(i,3)
    ENDDO
ENDDO
998 CONTINUE
CLOSE(51)

IF (type_dipole=='2' .OR. type_dipole=='3') THEN !!gas phase
    mol_num=1
ELSEIF ((periodic=='n' .AND. system=='1') .OR. type_dipole=='1') THEN !!fragment approach
    mol_num=44 !20 !! fix later to 20
ENDIF  
print*,mol_num,'mol num'
END SUBROUTINE read_coord

!********************************************************************************************
!********************************************************************************************

SUBROUTINE read_coord_frame(natom,framecount,element,filename,coord_v)

CHARACTER(LEN=40),INTENT(IN)                               :: filename
INTEGER,INTENT(INOUT)                                      :: natom,framecount
CHARACTER(LEN=2),DIMENSION(:),ALLOCATABLE,INTENT(INOUT)    :: element
REAL(kind=dp),DIMENSION(:,:,:),ALLOCATABLE,INTENT(OUT)      :: coord_v

INTEGER                                                    :: i,j,stat

ALLOCATE(coord_v(framecount,natom,3))
OPEN(UNIT=52,FILE=filename,STATUS='old',IOSTAT=stat)
DO 
    DO j=1,framecount
        READ(52,*,END=999) 
        READ(52,*) 
        DO i=1,natom
            READ(52,*) element(i),coord_v(j,i,1),coord_v(j,i,2),coord_v(j,i,3)   
        ENDDO
    ENDDO
ENDDO
999 CONTINUE
CLOSE(52)

END SUBROUTINE read_coord_frame

!********************************************************************************************
!********************************************************************************************

SUBROUTINE read_static(natom,element,normal_freq_file,normal_displ_file,static_pol,pol,freq,disp,&
         nmodes,static_dip_free_file,static_dip_x_file,static_dip_y_file,static_dip_z_file,&
         type_dipole,static_dipole_x,static_dipole_y,static_dipole_z,static_dipole_free,read_function,&
         type_static,force_file,force)

CHARACTER(LEN=40),INTENT(IN)                               :: static_pol,type_dipole,static_dip_z_file
CHARACTER(LEN=40),INTENT(IN)                               :: static_dip_free_file,static_dip_x_file
CHARACTER(LEN=40),INTENT(IN)                               :: static_dip_y_file,read_function,type_static
CHARACTER(LEN=40),INTENT(IN)                               :: normal_freq_file,normal_displ_file,force_file
INTEGER,INTENT(INOUT)                                      :: natom
INTEGER,INTENT(OUT)                                        :: nmodes
CHARACTER(LEN=2),DIMENSION(:),ALLOCATABLE,INTENT(INOUT)    :: element
REAL(kind=dp),DIMENSION(:,:,:,:,:),ALLOCATABLE,INTENT(OUT)  :: pol,force
REAL(kind=dp),DIMENSION(:,:,:,:),ALLOCATABLE,INTENT(OUT)    :: static_dipole_x,static_dipole_free
REAL(kind=dp),DIMENSION(:,:,:,:),ALLOCATABLE,INTENT(OUT)    :: static_dipole_y,static_dipole_z
REAL(kind=dp),DIMENSION(:),ALLOCATABLE,INTENT(OUT)          :: freq
REAL(kind=dp),DIMENSION(:,:,:),ALLOCATABLE,INTENT(OUT)      :: disp

CHARACTER(LEN=40)                                          :: chara
INTEGER                                                    :: i,j,k,m,n,data_number
INTEGER                                                    :: stat


IF (read_function=='NMA' .OR. type_static=='1') THEN
    ALLOCATE (force(2,natom,3,natom,3))
    OPEN(UNIT=49,FILE=force_file,STATUS='old',IOSTAT=stat) !Reading forces
    DO i=1,2
        DO j=1,natom
            DO m=1,3
                DO n=1,natom
                    READ(49,*) force(i,j,m,n,1),force(i,j,m,n,2),force(i,j,m,n,3)
                ENDDO
            ENDDO
        ENDDO
    ENDDO
    CLOSE(49)

print*,force(1,1,1,3,2),"force"

ELSEIF (type_static=='2') THEN
    nmodes=0
    OPEN(UNIT=50,FILE=normal_freq_file,STATUS='old',IOSTAT=stat) !Reading normal freqs/coords
    DO
        READ(50,*,END=998) chara
        nmodes=nmodes+1
    ENDDO
    998 CONTINUE
    CLOSE(50)
    
    ALLOCATE(freq(nmodes),disp(natom,nmodes,3))
    ALLOCATE(pol(natom,3,2,3,3),static_dipole_free(natom,3,2,3),static_dipole_x(natom,3,2,3))
    ALLOCATE(static_dipole_y(natom,3,2,3),static_dipole_z(natom,3,2,3))
    
    OPEN(UNIT=51,FILE=normal_freq_file,STATUS='old',IOSTAT=stat) !Reading normal freqs/coords
    DO i=1,nmodes
        READ(51,*,END=997) freq(i)
    ENDDO
    997 CONTINUE
    CLOSE(51)
    
    OPEN(UNIT=51,FILE=normal_displ_file,STATUS='old',IOSTAT=stat) !Reading normal freqs/coords
    DO i=1,nmodes       !look above!!!
        DO j=1,natom        !changed the order of these two lines for o-NP, for water it is reverse!!
            READ(51,*,END=996) disp(j,i,1),disp(j,i,2),disp(j,i,3)
        ENDDO
    ENDDO
    996 CONTINUE
    CLOSE(51)
ENDIF

IF (read_function=='IR' .OR. read_function=='R') THEN
    IF (type_dipole=='3') THEN
        OPEN(UNIT=52,FILE=static_pol,STATUS='old',IOSTAT=stat) !Reading polarizabilties
        DO
            DO k=1,2
                DO i=1,natom
                    DO j=1,3
                        READ(52,*,END=995) 
                        READ(52,*) 
                        READ(52,*) 
                        READ(52,*) 
                        READ(52,*) 
                        READ(52,*) 
                        READ(52,*) chara,chara,chara,pol(i,j,k,1,1),pol(i,j,k,2,2),pol(i,j,k,3,3)   
                        READ(52,*) chara,chara,chara,pol(i,j,k,1,2),pol(i,j,k,1,3),pol(i,j,k,2,3)   
                        READ(52,*) chara,chara,chara,pol(i,j,k,2,1),pol(i,j,k,3,1),pol(i,j,k,3,2)   
                    ENDDO
                ENDDO
            ENDDO
        ENDDO
        995 CONTINUE
        CLOSE(52)

    ELSEIF(type_dipole=='2') THEN
        OPEN(UNIT=53,FILE=static_dip_free_file,STATUS='old',IOSTAT=stat) !Reading dipoles
        DO
            DO k=1,2
                DO i=1,natom
                    DO j=1,3
                        READ(53,*,END=994) 
                        READ(53,*) 
                        READ(53,*) chara,static_dipole_free(i,j,k,1),static_dipole_free(i,j,k,2),static_dipole_free(i,j,k,3)   
                    ENDDO
                ENDDO
            ENDDO
        ENDDO
        994 CONTINUE
        CLOSE(53)

        print*,static_dipole_free(1,1,1,1),"dipole free test"
        OPEN(UNIT=54,FILE=static_dip_x_file,STATUS='old',IOSTAT=stat) !Reading dipoles
        DO
            DO k=1,2
                DO i=1,natom
                    DO j=1,3
                        READ(54,*,END=993) 
                        READ(54,*) 
                        READ(54,*) chara,static_dipole_x(i,j,k,1),static_dipole_x(i,j,k,2),static_dipole_x(i,j,k,3)   
                    ENDDO
                ENDDO
            ENDDO
        ENDDO
        993 CONTINUE
        CLOSE(54)
 
        print*,static_dipole_x(1,1,1,1),"dipole x test"
        OPEN(UNIT=55,FILE=static_dip_y_file,STATUS='old',IOSTAT=stat) !Reading dipoles

        DO
            DO k=1,2
                DO i=1,natom
                    DO j=1,3
                        READ(55,*,END=992) 
                        READ(55,*) 
                        READ(55,*) chara,static_dipole_y(i,j,k,1),static_dipole_y(i,j,k,2),static_dipole_y(i,j,k,3)   
                    ENDDO
                ENDDO
            ENDDO
        ENDDO
        992 CONTINUE
        CLOSE(55)
        print*,static_dipole_y(1,1,1,1),"dipole y test"
        OPEN(UNIT=56,FILE=static_dip_z_file,STATUS='old',IOSTAT=stat) !Reading dipoles
        DO
            DO k=1,2
                DO i=1,natom
                    DO j=1,3
                        READ(56,*,END=991) 
                        READ(56,*) 
                        READ(56,*) chara,static_dipole_z(i,j,k,1),static_dipole_z(i,j,k,2),static_dipole_z(i,j,k,3)   
                    ENDDO
                ENDDO
            ENDDO
        ENDDO
        991 CONTINUE
        CLOSE(56)
        print*,static_dipole_z(1,1,1,1),"dipole z test"
        print*,static_dipole_z(1,1,1,1),static_dipole_free(1,2,2,2),"dipoles"
    ENDIF
ENDIF

END SUBROUTINE read_static

!********************************************************************************************
!********************************************************************************************

SUBROUTINE read_static_resraman(natom,element,normal_freq_file,normal_displ_file,freq,disp,&
           nmodes,static_dip_x_file,static_dip_y_file,static_dip_z_file,framecount_rtp,&
                static_dipole_x_rtp,static_dipole_y_rtp,static_dipole_z_rtp)

CHARACTER(LEN=40),INTENT(IN)                               :: static_dip_z_file
CHARACTER(LEN=40),INTENT(IN)                               :: static_dip_y_file,static_dip_x_file
CHARACTER(LEN=40),INTENT(IN)                               :: normal_freq_file,normal_displ_file
INTEGER,INTENT(INOUT)                                      :: natom,framecount_rtp
INTEGER,INTENT(OUT)                                        :: nmodes
CHARACTER(LEN=2),DIMENSION(:),ALLOCATABLE,INTENT(INOUT)    :: element
REAL(kind=dp),DIMENSION(:,:,:,:,:),ALLOCATABLE,INTENT(OUT)  :: static_dipole_x_rtp
REAL(kind=dp),DIMENSION(:,:,:,:,:),ALLOCATABLE,INTENT(OUT)  :: static_dipole_y_rtp,static_dipole_z_rtp
REAL(kind=dp),DIMENSION(:),ALLOCATABLE,INTENT(OUT)          :: freq
REAL(kind=dp),DIMENSION(:,:,:),ALLOCATABLE,INTENT(OUT)      :: disp

CHARACTER(LEN=40)                                          :: chara
INTEGER                                                    :: i,j,k,m,data_number,step
INTEGER                                                    :: stat

nmodes=0
step=0
OPEN(UNIT=50,FILE=normal_freq_file,STATUS='old',IOSTAT=stat) !Reading normal freqs/coords
DO
 READ(50,*,END=998) chara
 nmodes=nmodes+1
ENDDO
998 CONTINUE
CLOSE(50)


print*,nmodes,'nmodes'

ALLOCATE(freq(nmodes),disp(natom,nmodes,3))
ALLOCATE(static_dipole_x_rtp(natom,3,2,3,framecount_rtp+1))
ALLOCATE(static_dipole_y_rtp(natom,3,2,3,framecount_rtp+1),static_dipole_z_rtp(natom,3,2,3,framecount_rtp+1))

OPEN(UNIT=51,FILE=normal_freq_file,STATUS='old',IOSTAT=stat) !Reading normal freqs/coords
DO i=1,nmodes
 READ(51,*,END=997) freq(i)
ENDDO
997 CONTINUE
CLOSE(51)

OPEN(UNIT=51,FILE=normal_displ_file,STATUS='old',IOSTAT=stat) !Reading normal freqs/coords
DO i=1,nmodes       !look above!!!
 DO j=1,natom        !changed the order of these two lines for o-NP, for water it is reverse!!
  READ(51,*,END=996) disp(j,i,1),disp(j,i,2),disp(j,i,3)
 ENDDO
ENDDO
996 CONTINUE
CLOSE(51)

print*,disp(1,2,1),'disp'
print*,disp(2,1,1),'disp'

PRINT*,step

OPEN(UNIT=53,FILE=static_dip_x_file,STATUS='old',IOSTAT=stat) !Reading polarizabilties
DO
 DO k=1,2
  DO i=1,natom
   DO j=1,3
    READ(53,*,END=994) 
    READ(53,*)
    DO m=1,framecount_rtp+1 
     READ(53,*) chara,static_dipole_x_rtp(i,j,k,1,m),static_dipole_x_rtp(i,j,k,2,m),static_dipole_x_rtp(i,j,k,3,m)
    ENDDO
   ENDDO 
  ENDDO
 ENDDO
ENDDO
994 CONTINUE
CLOSE(53)

PRINT*,"CHECKPOINT1"

OPEN(UNIT=54,FILE=static_dip_y_file,STATUS='old',IOSTAT=stat) !Reading polarizabilties
DO
DO k=1,2
 DO i=1,natom
  DO j=1,3
   READ(54,*,END=993) 
   READ(54,*)
   DO m=1,framecount_rtp+1 
    READ(54,*) chara,static_dipole_y_rtp(i,j,k,1,m),static_dipole_y_rtp(i,j,k,2,m),static_dipole_y_rtp(i,j,k,3,m)   
   ENDDO
 ENDDO
ENDDO
ENDDO
ENDDO
993 CONTINUE
CLOSE(54)

OPEN(UNIT=55,FILE=static_dip_z_file,STATUS='old',IOSTAT=stat) !Reading polarizabilties
DO
DO k=1,2
 DO i=1,natom
  DO j=1,3
   READ(55,*,END=992) 
   READ(55,*)
   DO m=1,framecount_rtp+1 
    READ(55,*) chara,static_dipole_z_rtp(i,j,k,1,m),static_dipole_z_rtp(i,j,k,2,m),static_dipole_z_rtp(i,j,k,3,m)   
   ENDDO
 ENDDO
ENDDO
ENDDO
ENDDO
992 CONTINUE
CLOSE(55)

print*,static_dipole_z_rtp(1,1,1,1,1),"dipole rtp z test"
print*,static_dipole_x_rtp(1,1,1,2,1),"dipole rtp x test"

END SUBROUTINE read_static_resraman
END MODULE read_traj 
