!
!   Copyright 2025 Ekin E. Winogradow, Johannes Scheffler, Moritz Leucke, Dorothea Golze
!
!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at
!
!       http://www.apache.org/licenses/LICENSE-2.0
!
!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.
!

!> @brief Module containing routines that are for reading of different types of trajectories.
MODULE read_traj

    USE kinds, ONLY: dp, str_len
    USE ISO_FORTRAN_ENV, ONLY: output_unit, error_unit
    USE vib_types, ONLY: global_settings, systems, static, dipoles, raman, molecular_dynamics, static_property
    USE output_io, ONLY: check_file_open

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: read_coord, read_coord_frame, read_normal_modes, read_static, read_static_resraman

CONTAINS
  
!****************************************************************************!
!****************************************************************************!

    !> @brief Reads atomic coordinates or the dipole moment files
    !>        and determines the number of atoms and frames.
    !>
    !> This subroutine reads atomic coordinate data from a given file (typically an `.xyz`
    !> trajectory or static geometry file) and initializes the corresponding system arrays.
    !> It can also read dipole moment files given in the format of .xyz files to determine 
    !> the number of MD frames.
    !>
    !> @param[in]     filename  -- Path to the coordinate file (e.g., `.xyz` format).  
    !> @param[in,out] gs        -- Global settings (provides `spectral_type%read_function`).  
    !> @param[in,out] sys       -- System information (allocates `element` and `coord` arrays).  
    !> @param[in,out] dips      -- Dipole data structure (optional; provides `type_dipole`). 
    !> @param[in,out] rams      -- Raman data structure (optional; provides frame count for RTP modes).  
    !>
    SUBROUTINE read_coord(filename, gs, sys, dips, rams)

        TYPE(global_settings), INTENT(INOUT)   :: gs
        TYPE(systems), INTENT(INOUT)        :: sys
        TYPE(dipoles), OPTIONAL        :: dips
        TYPE(raman), OPTIONAL        :: rams
        CHARACTER(LEN=40), INTENT(IN)                               :: filename

        CHARACTER(len=str_len)                                     :: msg  !store error message
        INTEGER                                                   :: i, j, stat, runit, ios

        !Initialize 
        sys%framecount = 0

        !!If the file contains cartesian coordinates or static dipole moments
        IF (gs%spectral_type%read_function/='MD-RR' .AND. gs%spectral_type%read_function/='MD-ABS') THEN
            OPEN (FILE=filename, STATUS='old', ACTION='read', IOSTAT=stat, IOMSG=msg, NEWUNIT=runit)
            !Check if file exists
            CALL check_file_open(stat, msg, filename, runit)
            READ (runit, *) sys%natom
            CLOSE (runit)
        !!If the file contains dynamic time-dependent dipole moments
        ELSEIF (gs%spectral_type%read_function=='MD-RR' .OR. gs%spectral_type%read_function=='MD-ABS') THEN
            rams%RR%framecount_rtp = rams%RR%framecount_rtp + 1
            sys%natom = rams%RR%framecount_rtp
        END IF
        
        !!Allocate
        ALLOCATE (sys%element(sys%natom), sys%coord(sys%natom, 3))
 
        !!Read the file to assign the coordinates and elements, determine the MD frame count
        OPEN (FILE=filename, STATUS='old', ACTION='read', IOSTAT=stat, IOMSG=msg, NEWUNIT=runit)
        !Check if file exists
        CALL check_file_open(stat, msg, filename, runit)
        DO
            READ (runit, *, END=998)
            READ (runit, *)
            sys%framecount = sys%framecount + 1
            DO i = 1, sys%natom
                READ (runit, *) sys%element(i), sys%coord(i, 1), sys%coord(i, 2), sys%coord(i, 3)
            END DO
        END DO
998     CONTINUE
        CLOSE (runit)

        !!If dipole moment of the whole system is considered (no fragment approach or molecules defined)
        IF (gs%spectral_type%read_function/='P') THEN
            IF (dips%type_dipole=='berry' .OR. dips%type_dipole=='dfpt' .OR. dips%type_dipole=='wannier') THEN !!gas phase
                sys%mol_num = 1
            END IF
        END IF

    END SUBROUTINE read_coord

!********************************************************************************************
!********************************************************************************************
    !> @brief Reads all coordinate frames from a trajectory containing coordinates, dipole moments
    !>        or polarizability tensors
    !>
    !> This subroutine reads the full set of Cartesian coordinates or dipole moment 
    !> vectors polarizability tensors for each MD frame in a trajectory file 
    !> (typically in XYZ format) and stores them in a 3D array `coord_v(framecount, natom, 3)`.
    !>
    !> @param[in,out] sys       -- System information (provides `framecount` and optional element labels).  
    !> @param[in]     filename  -- Path to the coordinate file (e.g., `.xyz` trajectory).  
    !> @param[in,out] natom     -- Number of atoms for cartesian coordinates, or just `1`
    !>                             if the quantity (e.g. dipole or polarizability) refers
    !>                             to the entire system.
    !> @param[out]    coord_v   -- 3D array of coordinates with dimensions `(framecount, natom, 3)`.  
    !>
    SUBROUTINE read_coord_frame(natom, filename, coord_v, sys)

        TYPE(systems), INTENT(INOUT)        :: sys
        CHARACTER(LEN=40), INTENT(IN)                               :: filename
        INTEGER, INTENT(INOUT)                               :: natom
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE, INTENT(OUT)      :: coord_v

        CHARACTER(len=str_len)                                     :: msg  ! store error message
        INTEGER                                                    :: i, j, stat, runit, buff, ios

        !!Allocate
        ALLOCATE (coord_v(sys%framecount, natom, 3))
        
        OPEN (FILE=filename, STATUS='old', ACTION='read', IOSTAT=stat, IOMSG=msg, NEWUNIT=runit)
        !Check if file exists
        CALL check_file_open(stat, msg, filename, runit)
        !Start reading until the end of the file is reached
        DO
            DO j = 1, sys%framecount
                READ (runit, *, END=999)
                READ (runit, *)
                DO i = 1, natom
                    IF (ALLOCATED(sys%element)) THEN
                        READ (runit, *) sys%element(i), coord_v(j, i, 1), coord_v(j, i, 2), coord_v(j, i, 3)
                    ELSE
                        READ (runit, *) buff, coord_v(j, i, 1), coord_v(j, i, 2), coord_v(j, i, 3)
                    END IF
                END DO

            END DO
        END DO
999     CONTINUE
        CLOSE (runit)

    END SUBROUTINE read_coord_frame
!********************************************************************************************
!********************************************************************************************
    !> @brief Reads normal mode information from forces or external files that contain
    !>        normal mode data.
    !> This subroutine initializes and reads the normal mode data required for
    !> static vibrational spectrum calculations. Depending on the calculation
    !> setup, it either:
    !> - Reads raw force data for Hessian diagonalization, or  
    !> - Reads precomputed normal mode frequencies and displacement vectors.  
    !>
    !> @param[in,out] gs    --  Global settings (provides spectral type selection).  
    !> @param[in,out] sys   --  System information (provides natom).  
    !> @param[in,out] stats --  Static data structure (provides normal mode, Hessian, and
    !>                          finite displacement information).  
    !>
    SUBROUTINE read_normal_modes(gs, sys, stats)

        ! Variables of your derived types:
        TYPE(global_settings), INTENT(INOUT)   :: gs
        TYPE(systems), INTENT(INOUT)        :: sys
        TYPE(static), INTENT(INOUT)        :: stats

        CHARACTER(LEN=str_len)                                          :: chara, msg
        INTEGER                                                    :: i, j, k, m, n, d, xyz, ios
        INTEGER                                                    :: stat, runit

        !!If the user requested a normal mode analysis based on diagonalizing the hessian
        IF (gs%spectral_type%read_function=='NMA' .OR. stats%diag_hessian=='y') THEN

            CALL stats%init_force(sys%natom, 1)
            
            OPEN (FILE=stats%force_file, STATUS='old', ACTION='read', IOSTAT=stat, IOMSG=msg, NEWUNIT=runit)
            !Check if file exists
            CALL check_file_open(stat, msg, stats%force_file, runit)
            !Read forces
            DO i = 1, 2
                DO j = 1, sys%natom
                    DO m = 1, 3
                        DO n = 1, sys%natom
                            READ (runit, *) stats%force(n, 1)%atom(j)%displacement(i)%XYZ(m)%frame(1), &
                                stats%force(n, 2)%atom(j)%displacement(i)%XYZ(m)%frame(1), &
                                stats%force(n, 3)%atom(j)%displacement(i)%XYZ(m)%frame(1)
                        END DO
                    END DO
                END DO
            END DO
            CLOSE (runit)

        !!If the user already provided the files containing the normal mode frequencies and displacements
        ELSEIF (stats%diag_hessian=='n') THEN
            stats%nmodes = 0
            OPEN (FILE=stats%normal_freq_file, STATUS='old', ACTION='read', IOSTAT=stat, IOMSG=msg, NEWUNIT=runit) !Reading normal freqs/coords
            !Check if file exists
            CALL check_file_open(stat, msg, stats%normal_freq_file, runit)
            !First determine the number of modes from the provided frequency file
            DO
                READ (runit, *, END=998) chara
                stats%nmodes = stats%nmodes + 1
            END DO
998         CONTINUE
            CLOSE (runit)
            !Allocate
            ALLOCATE (stats%freq(stats%nmodes), stats%disp(stats%nmodes, sys%natom, 3))

            !!Read frequencies
            OPEN (FILE=stats%normal_freq_file, STATUS='old', ACTION='read', IOSTAT=stat, IOMSG=msg, NEWUNIT=runit) !Reading normal freqs/coords
            !Check if file exists
            CALL check_file_open(stat, msg, stats%normal_freq_file, runit)
            DO i = 1, stats%nmodes
                READ (runit, *, END=997) stats%freq(i)
            END DO
997         CONTINUE
            CLOSE (runit)

            !!Read normal mode coordinates
            OPEN (FILE=stats%normal_displ_file, STATUS='old', ACTION='read', IOSTAT=stat, IOMSG=msg, NEWUNIT=runit) !Reading normal freqs/coords
            !Check if file exists
            CALL check_file_open(stat, msg, stats%normal_displ_file, runit)
            DO i = 1, stats%nmodes
                DO j = 1, sys%natom
                    READ (runit, *, END=996) stats%disp(i, j, 1), stats%disp(i, j, 2), stats%disp(i, j, 3)
                END DO
            END DO
996         CONTINUE
            CLOSE (runit)
        END IF

    END SUBROUTINE read_normal_modes

!********************************************************************************************
!********************************************************************************************
    !> @brief Reads static dipole moments or polarizabilities to be used in IR and Raman
    !>        calculations.
    !>
    !> @param[in,out] gs   -- Global settings structure (provides spectral context).  
    !> @param[in,out] sys  -- System information (provides natom).  
    !> @param[in,out] dips -- Dipole-related quantities (optional).  
    !> @param[in,out] rams -- Raman-related quantities (optional).  
    !>
    SUBROUTINE read_static(gs, sys, dips, rams)
        TYPE(global_settings), INTENT(INOUT)   :: gs
        TYPE(systems), INTENT(INOUT)        :: sys
        TYPE(dipoles), INTENT(INOUT), OPTIONAL        :: dips
        TYPE(raman), INTENT(INOUT), OPTIONAL        :: rams

        CHARACTER(LEN=str_len)                                          :: chara, msg

        INTEGER                                                    :: i, j, k, m, n, d, runit
        INTEGER                                                    :: stat, i_pol, j_pol, xyz

        !IF (PRESENT(static_dip)) ALLOCATE (static_dip(sys%natom, 3, 2, 3))
        IF (PRESENT(dips)) CALL dips%init_dip(sys%natom, 1)
        IF (PRESENT(rams)) CALL rams%init_pol(sys%natom, 1)

        !!If the user provided DFPT polarizabilities given in a specific CP2K style:
        IF (dips%type_dipole=='dfpt') THEN
            OPEN (FILE=rams%static_pol_file, STATUS='old', ACTION='read', IOSTAT=stat, IOMSG=msg, NEWUNIT=runit)
            !Check if file exists
            CALL check_file_open(stat, msg, rams%static_pol_file, runit)
            DO
                DO k = 1, 2 !+/-
                    DO i = 1, sys%natom
                        DO j = 1, 3 !x,y,z
                            READ (runit, *, END=995)
                            READ (runit, *)
                            READ (runit, *)
                            READ (runit, *)
                            READ (runit, *)
                            READ (runit, *)
                            READ (runit, *) chara, chara, chara, rams%pol(1, 1)%atom(i)%displacement(k)%XYZ(j)%frame(1), &
                                rams%pol(2, 2)%atom(i)%displacement(k)%XYZ(j)%frame(1), rams%pol(3, 3)%atom(i)%displacement(k)%XYZ(j)%frame(1)
                            READ (runit, *) chara, chara, chara, rams%pol(1, 2)%atom(i)%displacement(k)%XYZ(j)%frame(1), &
                                rams%pol(1, 3)%atom(i)%displacement(k)%XYZ(j)%frame(1), rams%pol(2, 3)%atom(i)%displacement(k)%XYZ(j)%frame(1)
                            READ (runit, *) chara, chara, chara, rams%pol(2, 1)%atom(i)%displacement(k)%XYZ(j)%frame(1), &
                                rams%pol(3, 1)%atom(i)%displacement(k)%XYZ(j)%frame(1), rams%pol(3, 2)%atom(i)%displacement(k)%XYZ(j)%frame(1)
                        END DO
                    END DO
                END DO
            END DO
995         CONTINUE
            CLOSE (runit)

        !!If the user provided berry phase dipole moments:
        ELSEIF (dips%type_dipole=='berry') THEN
            OPEN (FILE=dips%dip_file, STATUS='old', ACTION='read', IOSTAT=stat, IOMSG=msg, NEWUNIT=runit)
            !Check if file exists
            CALL check_file_open(stat, msg, dips%dip_file, runit)
            DO
                DO k = 1, 2 !+/-
                    DO i = 1, sys%natom
                        DO j = 1, 3 !x,y,z
                            READ (runit, *, END=994)
                            READ (runit, *)
                            READ (runit, *) chara, dips%static_dip(1)%atom(i)%displacement(k)%XYZ(j)%frame(1), &
                                dips%static_dip(2)%atom(i)%displacement(k)%XYZ(j)%frame(1), &
                                dips%static_dip(3)%atom(i)%displacement(k)%XYZ(j)%frame(1)
                        END DO
                    END DO
                END DO
            END DO
994         CONTINUE
            CLOSE (runit)

        END IF

    END SUBROUTINE read_static

!********************************************************************************************
!********************************************************************************************

    !> @brief Reads time-dependent dipole data for resonance Raman calculations.
    !>
    !> This subroutine reads static dipole trajectories from an RT-TDDFT calculation, 
    !> for each atomic displacement and polarization direction.
    !> The data are stored in the `static_dip_rtp` array, which is used for
    !> computing resonance Raman response functions.
    !>
    !> @param[in]  static_dip_file -- Path to the file containing time-dependent dipole data.  
    !> @param[out] static_dip_rtp  -- Array of static dipole data for x, y, and z polarizations.  
    !> @param[in,out] sys          -- System information (provides natom).  
    !> @param[in,out] rams         -- Raman calculation data (provides framecount and structure definitions).  
    !>
    SUBROUTINE read_static_resraman(static_dip_file, static_dip_rtp, sys, rams)

        TYPE(systems), INTENT(INOUT)        :: sys
        TYPE(raman), INTENT(INOUT)        :: rams
        CHARACTER(LEN=40), INTENT(IN)                               :: static_dip_file
        TYPE(static_property), DIMENSION(3), INTENT(OUT)  :: static_dip_rtp

        CHARACTER(LEN=str_len)                                          :: chara, msg
        INTEGER                                                    :: x, y, i, j, k, m, stat, stat2, runit

        !! Allocate
        CALL rams%RR%init_rr_static_dip(static_dip_rtp, sys%natom, rams%RR%framecount_rtp + 1)

        !! Read the time-dependent dipoles for each displaced structure
        OPEN (FILE=static_dip_file, STATUS='old', ACTION='read', IOSTAT=stat, IOMSG=msg, NEWUNIT=runit) !Reading polarizabilties
        !Check if file exists
        CALL check_file_open(stat, msg, static_dip_file, runit)
        DO
            DO k = 1, 2 !!+/- 
                DO i = 1, sys%natom
                    DO j = 1, 3 !!x, y, z
                        READ (runit, *, END=994)
                        READ (runit, *)
                        DO m = 1, rams%RR%framecount_rtp + 1
                            READ (runit, *) chara, static_dip_rtp(1)%atom(i)%displacement(k)%XYZ(j)%frame(m), &
                                static_dip_rtp(2)%atom(i)%displacement(k)%XYZ(j)%frame(m), &
                                static_dip_rtp(3)%atom(i)%displacement(k)%XYZ(j)%frame(m)

                        END DO
                    END DO
                END DO
            END DO
        END DO
994     CONTINUE
        CLOSE (runit)
    END SUBROUTINE read_static_resraman
END MODULE read_traj
