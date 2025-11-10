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

!> @brief Module containing all procedures involving autocorrelation functions.
MODULE vel_cor

    USE kinds, ONLY: dp, str_len
    USE ISO_FORTRAN_ENV, ONLY: output_unit, error_unit
    USE constants, ONLY: pi, am_u, ang, fs2s, at_u, bohr2ang
    USE vib_types, ONLY: global_settings, systems, molecular_dynamics, static, dipoles, raman
    USE output_io, ONLY: check_file_open

    IMPLICIT NONE
    PUBLIC :: cvv, cvv_iso, cvv_aniso, cvv_resraman

CONTAINS

    !> @brief Computes the autocorrelation function of vector quantities (velocity or dipole) over time.
    !>
    !> This routine evaluates the time-correlation function of atomic or molecular
    !> quantities stored in `coord_v` (e.g. velocities, positions, or dipole moments),
    !> over a correlation time `t_cor`. Unit conversions and a Hann window are applied
    !> depending on the selected spectral type. The function is later used to obtain vibrational
    !> spectra via Fourier transformation.
    !>
    !> @param[in,out] sys     -- System information (provides atomic masses, trajectory type, etc.).
    !> @param[in,out] gs      -- Global settings (provides `spectral_type%read_function`).
    !> @param[in,out] md      -- Molecular dynamics parameters and arrays (`t_cor`, `z`, etc.).
    !> @param[in,out] natom   -- Number of atoms for cartesian coordinates, or just `1`
    !>                           if the quantity (e.g. dipole or polarizability) refers
    !>                           to the entire system.
    !> @param[in,out] coord_v -- 3D array of vector data (positions, velocities, or dipoles),
    !>                           dimensions (nframes, natom, 3). Used as input for the correlation.
    !>
    SUBROUTINE cvv(natom, coord_v, sys, gs, md)

        TYPE(systems), INTENT(INOUT)                :: sys
        TYPE(global_settings), INTENT(INOUT)                :: gs
        TYPE(molecular_dynamics), INTENT(INOUT)     :: md
        INTEGER, INTENT(INOUT)                                    :: natom
        REAL(kind=dp), DIMENSION(:, :, :), INTENT(INOUT)  :: coord_v

        CHARACTER(LEN=40)                                        :: chara, msg
        INTEGER                                                  :: stat, i, j, k, m, t0, t1, l, runit
        INTEGER, DIMENSION(:), ALLOCATABLE                         :: norm
        REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE                  :: coord

        !! Allocate
        ALLOCATE (md%z(0:2*md%t_cor - 1), norm(0:2*md%t_cor - 1))

        !! Initialize
        norm = 0
        md%z = 0.0_dp
        k = 0
        j = 0

        !! Start the time autocorrelations for power of IR spectra
        DO t0 = 1, sys%framecount - 2
            t1 = MIN(sys%framecount - 2, t0 + md%t_cor)
            DO j = 1, natom
                k = j
                !! If mass-weighting is selected for power spectrum
                IF (sys%input_mass=='y') THEN
                    md%z(0:t1 - t0) = md%z(0:t1 - t0) + (coord_v(t0, k, 1)*coord_v(t0:t1, j, 1) + coord_v(t0, j, 2)* &
                                                         coord_v(t0:t1, j, 2) + coord_v(t0, j, 3)*coord_v(t0:t1, j, 3))*sys%mass_atom(j)
                !! If mass-weighting is not selected for power spectrum
                ELSE
                    DO m = 1, 3
                        md%z(0:t1 - t0) = md%z(0:t1 - t0) + coord_v(t0, k, m)*coord_v(t0:t1, k, m)
                    END DO
                END IF
            END DO
            norm(0:t1 - t0) = norm(0:t1 - t0) + 1
        END DO

        !! Normalization
        md%z(:) = md%z(:)/norm(:)

      !! Unit conversion of dipole autocorrelation function to debye^2/s^2
        IF (gs%spectral_type%read_function=='MD-IR') THEN
            md%z(:) = md%z(:)/(fs2s*fs2s)
      !! Unit conversion of velocity autocorrelation function to m^2/s^2, or m^2 /(s^2 kg) if not mass-weighted
        ELSE IF (gs%spectral_type%read_function=='P') THEN
            IF (sys%type_traj=='pos') THEN
                md%z(:) = md%z(:)*ang*ang/(fs2s*fs2s)
            ELSE IF (sys%type_traj=='vel') THEN
                md%z(:) = md%z(:)*ang*ang*bohr2ang*bohr2ang/(at_u*at_u)
            END IF
            IF (sys%input_mass=='y') THEN
                md%z(:) = md%z(:)*am_u
            END IF
        END IF

        !! Data mirroring
        md%z(md%t_cor) = 0.0_dp
        DO i = 1, md%t_cor - 1
            md%z(md%t_cor + i) = md%z(md%t_cor - i)
        END DO

        !! Hann Window function
        DO i = 0, 2*md%t_cor - 1
            md%z(i) = md%z(i)*((COS(i/(md%t_cor - 1.0_dp)/2.0_dp*pi))**2)
        END DO

        DEALLOCATE (norm)

    END SUBROUTINE cvv

!***********************************************************************************!
!***********************************************************************************!

    !> @brief Computes the isotropic autocorrelation function of Raman polarizabilities.
    !>
    !> This routine evaluates the time autocorrelation of the isotropic part of the
    !> polarizability tensor using the time-dependent differences of polarizabilities
    !> along the x, y, and z directions. The resulting correlation function `z_iso` is normalized,
    !> windowed, and later Fourier transformed to be used for computing the Raman intensities.
    !>
    !> @param[in,out] sys          -- System information (provides framecount, etc.).
    !> @param[in,out] md           -- Molecular dynamics data (provides correlation time `t_cor`).
    !> @param[in,out] mol_num      -- Number of molecules, for cartesian coordinates it is equal to natom,
    !>                                if the dipole moment of the whole system is considered, it is `1`.
    !> @param[out]    z_iso        -- Output of isotropic autocorrelation function array,
    !>                                dimensions (0 : 2*md%t_cor - 1).
    !> @param[in,out] alpha_diff_x -- 3D array of time-dependent polarizability differences
    !>                                along x, dimensions (nframes, mol_num, 3).
    !> @param[in,out] alpha_diff_y -- 3D array of time-dependent polarizability differences
    !>                                along y, dimensions (nframes, mol_num, 3).
    !> @param[in,out] alpha_diff_z -- 3D array of time-dependent polarizability differences
    !>                                along z, dimensions (nframes, mol_num, 3).
    !>
    SUBROUTINE cvv_iso(mol_num, z_iso, alpha_diff_x, alpha_diff_y, alpha_diff_z, sys, md)

        TYPE(systems), INTENT(INOUT)                :: sys
        TYPE(molecular_dynamics), INTENT(INOUT)     :: md
        INTEGER, INTENT(INOUT)                                    :: mol_num
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(OUT)        :: z_iso
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE, INTENT(INOUT)  :: alpha_diff_x, alpha_diff_y, alpha_diff_z

        INTEGER                                                  :: stat, i, j, k, t0, t1, runit
        CHARACTER(len=str_len)                                  :: msg
        INTEGER, DIMENSION(:), ALLOCATABLE                         :: norm

        !! Allocate
        ALLOCATE (z_iso(0:2*md%t_cor - 1), norm(0:2*md%t_cor - 1))

        !! Initialize
        norm = 0
        z_iso = 0.0_dp

        !! Start the time autocorrelations for the isotropic Raman polarizabilities
        DO t0 = 1, sys%framecount - 2
            t1 = MIN(sys%framecount - 2, t0 + md%t_cor)
            DO j = 1, mol_num
                k = j
                z_iso(0:t1 - t0) = z_iso(0:t1 - t0) + (alpha_diff_x(t0, k, 1) + alpha_diff_y(t0, k, 2) + alpha_diff_z(t0, k, 3))* &
                                   (alpha_diff_x(t0:t1, k, 1) + alpha_diff_y(t0:t1, k, 2) + alpha_diff_z(t0:t1, k, 3))
            END DO
            norm(0:t1 - t0) = norm(0:t1 - t0) + 1
        END DO

        !! Normalization
        z_iso(:) = z_iso(:)/(norm(:)*9._dp)

      !!Unit conversion of Debye^2/(E^2*fs^2) into Debye^2/(E^2*s^2)
        z_iso(:) = z_iso(:)/(fs2s*fs2s)

        !! Hann Window function
        DO i = 0, md%t_cor - 1
            z_iso(i) = z_iso(i)*((COS(i/(md%t_cor - 1.0_dp)/2.0_dp*pi))**2)
        END DO

        !! Data mirroring
        z_iso(md%t_cor) = 0.0_dp
        DO i = 1, md%t_cor - 1
            z_iso(md%t_cor + i) = z_iso(md%t_cor - i)
        END DO

        DEALLOCATE (norm)

    END SUBROUTINE cvv_iso

!***********************************************************************************!
!***********************************************************************************!

    !> @brief Computes the anisotropic autocorrelation function of Raman polarizabilities.
    !>
    !> This routine evaluates the time autocorrelation of the anisotropic part of the
    !> polarizability tensor using the time-dependent differences of polarizabilities
    !> along the x, y, and z directions. The resulting correlation function `z_aniso` is normalized,
    !> windowed, and later Fourier transformed to be used for computing the Raman intensities.
    !>
    !> @param[in,out] sys          -- System information (provides framecount, etc.).
    !> @param[in,out] md           -- Molecular dynamics data (provides correlation time `t_cor`).
    !> @param[in,out] mol_num      -- Number of molecules, for cartesian coordinates it is equal to natom,
    !>                                if the dipole moment of the whole system is considered, it is `1`.
    !> @param[out]    z_aniso      -- Output of anisotropic autocorrelation function array,
    !>                                dimensions (0 : 2*md%t_cor - 1).
    !> @param[in,out] alpha_diff_x -- 3D array of time-dependent polarizability differences
    !>                                along x, dimensions (nframes, mol_num, 3).
    !> @param[in,out] alpha_diff_y -- 3D array of time-dependent polarizability differences
    !>                                along y, dimensions (nframes, mol_num, 3).
    !> @param[in,out] alpha_diff_z -- 3D array of time-dependent polarizability differences
    !>                                along z, dimensions (nframes, mol_num, 3).
    !>
    SUBROUTINE cvv_aniso(mol_num, z_aniso, alpha_diff_x, alpha_diff_y, alpha_diff_z, sys, md)

        TYPE(systems), INTENT(INOUT)                :: sys
        TYPE(molecular_dynamics), INTENT(INOUT)     :: md
        INTEGER, INTENT(INOUT)                                    :: mol_num
        REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(OUT)        :: z_aniso
        REAL(kind=dp), DIMENSION(:, :, :), ALLOCATABLE, INTENT(INOUT)  :: alpha_diff_x, alpha_diff_y, alpha_diff_z

        INTEGER                                                  :: stat, i, j, k, t0, t1, runit
        CHARACTER(len=str_len)                                  :: msg
        INTEGER, DIMENSION(:), ALLOCATABLE                         :: norm

        !! Allocate
        ALLOCATE (z_aniso(0:2*md%t_cor - 1), norm(0:2*md%t_cor - 1))

        !! Initialize
        norm = 0
        z_aniso = 0.0_dp

        !! Start the time autocorrelations for the anisotropic Raman polarizabilities
        DO t0 = 1, sys%framecount - 2
            t1 = MIN(sys%framecount - 2, t0 + md%t_cor)
            DO j = 1, mol_num
                k = j
                z_aniso(0:t1 - t0) = z_aniso(0:t1 - t0) + (alpha_diff_x(t0, k, 1) - alpha_diff_y(t0, k, 2)) &
                                     *(alpha_diff_x(t0:t1, k, 1) - alpha_diff_y(t0:t1, k, 2))/2.0_dp
                z_aniso(0:t1 - t0) = z_aniso(0:t1 - t0) + (alpha_diff_y(t0, k, 2) - alpha_diff_z(t0, k, 3)) &
                                     *(alpha_diff_y(t0:t1, k, 2) - alpha_diff_z(t0:t1, k, 3))/2.0_dp
                z_aniso(0:t1 - t0) = z_aniso(0:t1 - t0) + (alpha_diff_z(t0, k, 3) - alpha_diff_x(t0, k, 1)) &
                                     *(alpha_diff_z(t0:t1, k, 3) - alpha_diff_x(t0:t1, k, 1))/2.0_dp
                z_aniso(0:t1 - t0) = z_aniso(0:t1 - t0) + (alpha_diff_x(t0, k, 2)*0.50_dp + alpha_diff_y(t0, k, 1)*0.50_dp) &
                                     *(alpha_diff_x(t0:t1, k, 2)*0.50_dp + alpha_diff_y(t0:t1, k, 1)*0.50_dp)*3.0_dp
                z_aniso(0:t1 - t0) = z_aniso(0:t1 - t0) + (alpha_diff_y(t0, k, 3)*0.50_dp + alpha_diff_z(t0, k, 2)*0.50_dp) &
                                     *(alpha_diff_y(t0:t1, k, 3)*0.50_dp + alpha_diff_z(t0:t1, k, 2)*0.50_dp)*3.0_dp
                z_aniso(0:t1 - t0) = z_aniso(0:t1 - t0) + (alpha_diff_z(t0, k, 1)*0.50_dp + alpha_diff_x(t0, k, 3)*0.50_dp) &
                                     *(alpha_diff_z(t0:t1, k, 1)*0.50_dp + alpha_diff_x(t0:t1, k, 3)*0.50_dp)*3.0_dp
            END DO
            norm(0:t1 - t0) = norm(0:t1 - t0) + 1
        END DO

        !! Normalization
        z_aniso(:) = z_aniso(:)/norm(:)

      !!Unit conversion of Debye^2/(E^2*fs^2) into C^4*s^2/kg^2
        z_aniso(:) = z_aniso(:)/(fs2s*fs2s)

        !! Hann Window function
        DO i = 0, md%t_cor - 1
            z_aniso(i) = z_aniso(i)*((COS(i/(md%t_cor - 1.0_dp)/2.0_dp*pi))**2)
        END DO

        !! Data mirroring
        z_aniso(md%t_cor) = 0.0_dp
        DO i = 1, md%t_cor - 1
            z_aniso(md%t_cor + i) = z_aniso(md%t_cor - i)
        END DO

        DEALLOCATE (norm)

    END SUBROUTINE cvv_aniso

!**********************************************************************************!
!***********************************************************************************!

    !> @brief Computes isotropic and anisotropic complex autocorrelation functions
    !>        of time-dependent polarizability derivatives for resonance Raman spectra.
    !>
    !> This routine constructs both the isotropic (`z_iso_resraman`) and anisotropic
    !> (`z_aniso_resraman`) time autocorrelation functions of the real-time
    !> polarizability response functions obtained from RT-TDDFT simulations. The calculated
    !> isotropic and anisotropic polarizabilities are later Fourier transformed to yield
    !> the resonance Raman spectra.
    !>
    !> @param[in,out] sys  -- System information (natom, framecount, etc.).
    !> @param[in,out] md   -- Molecular dynamics data (provides correlation time `t_cor`).
    !> @param[in,out] rams -- Raman data structure containing:
    !>                        - `rams%RR%alpha_resraman_[xyz]_diff_re/im`:
    !>                           real and imaginary polarizability derivatives,
    !>                        - `rams%RR%z_iso_resraman`, `z_aniso_resraman`:
    !>                           output complex correlation functions,
    !>                        - `rams%RR%check_pade`: flag for Padé interpolation,
    !>                           used to determine frequency grid size `Nw`.
    !>
    SUBROUTINE cvv_resraman(sys, md, rams)

        TYPE(systems), INTENT(INOUT)     :: sys
        TYPE(molecular_dynamics), INTENT(INOUT)     :: md
        TYPE(raman), INTENT(INOUT)     :: rams

        CHARACTER(len=str_len)                                  :: msg
        INTEGER                                                  :: stat, i, j, k, m, t0, t1, runit, Nw
        INTEGER, DIMENSION(:, :), ALLOCATABLE                       :: norm_iso, norm_aniso
        COMPLEX(kind=dp)                                          :: im_unit

        !! Determine Nw based on the flag for Padé interpolation
        IF (rams%RR%check_pade=='n') THEN
            Nw = rams%RR%framecount_rtp - 1
        ELSEIF (rams%RR%check_pade=='y') THEN
            Nw = rams%RR%framecount_rtp_pade
        END IF

        !! Allocate
        ALLOCATE (rams%RR%z_iso_resraman(0:2*md%t_cor, Nw), norm_iso(0:2*md%t_cor, Nw))
        ALLOCATE (rams%RR%z_aniso_resraman(0:2*md%t_cor, Nw), norm_aniso(0:2*md%t_cor, Nw))

        !! Initialize
        im_unit = (0.0_dp, 1.0_dp)
        rams%RR%z_iso_resraman = (0.0_dp, 0.0_dp)
        rams%RR%z_aniso_resraman = (0.0_dp, 0.0_dp)
        norm_iso = 0.0_dp
        norm_aniso = 0.0_dp
        sys%framecount = sys%framecount - 2

        !! Start complex autocorrelations for the isotropic part (alpha_xx + alpha_yy + alpha_zz)
        DO t0 = 2, sys%framecount
            t1 = MIN(sys%framecount, t0 + md%t_cor)
            DO k = 1, Nw
                !!Real * Real
                rams%RR%z_iso_resraman(0:t1 - t0, k) = rams%RR%z_iso_resraman(0:t1 - t0, k) &
                                                       + (rams%RR%alpha_resraman_x_diff_re(t0, k, 1) + rams%RR%alpha_resraman_y_diff_re(t0, k, 2) &
                                                          + rams%RR%alpha_resraman_z_diff_re(t0, k, 3)) &
                                                       *(rams%RR%alpha_resraman_x_diff_re(t0:t1, k, 1) &
                                                         + rams%RR%alpha_resraman_y_diff_re(t0:t1, k, 2) + rams%RR%alpha_resraman_z_diff_re(t0:t1, k, 3))
                !!Imaginary * Imaginary
                rams%RR%z_iso_resraman(0:t1 - t0, k) = rams%RR%z_iso_resraman(0:t1 - t0, k) &
                                                       + (rams%RR%alpha_resraman_x_diff_im(t0, k, 1) + rams%RR%alpha_resraman_y_diff_im(t0, k, 2) + &
                                                          rams%RR%alpha_resraman_z_diff_im(t0, k, 3)) &
                                                       *(rams%RR%alpha_resraman_x_diff_im(t0:t1, k, 1) &
                                                         + rams%RR%alpha_resraman_y_diff_im(t0:t1, k, 2) + rams%RR%alpha_resraman_z_diff_im(t0:t1, k, 3))
                !!Real * Imaginary
                rams%RR%z_iso_resraman(0:t1 - t0, k) = rams%RR%z_iso_resraman(0:t1 - t0, k) &
                                                       + ((rams%RR%alpha_resraman_x_diff_re(t0, k, 1) + rams%RR%alpha_resraman_y_diff_re(t0, k, 2) &
                                                           + rams%RR%alpha_resraman_z_diff_re(t0, k, 3)) &
                                                          *(rams%RR%alpha_resraman_x_diff_im(t0:t1, k, 1) &
                                                            + rams%RR%alpha_resraman_y_diff_im(t0:t1, k, 2) &
                                                            + rams%RR%alpha_resraman_z_diff_im(t0:t1, k, 3)))*im_unit
                !!Real * Imaginary (SUBTRACT)
                rams%RR%z_iso_resraman(0:t1 - t0, k) = rams%RR%z_iso_resraman(0:t1 - t0, k) &
                                                       - ((rams%RR%alpha_resraman_x_diff_im(t0, k, 1) + rams%RR%alpha_resraman_y_diff_im(t0, k, 2) &
                                                           + rams%RR%alpha_resraman_z_diff_im(t0, k, 3)) &
                                                          *(rams%RR%alpha_resraman_x_diff_re(t0:t1, k, 1) + rams%RR%alpha_resraman_y_diff_re(t0:t1, k, 2) &
                                                            + rams%RR%alpha_resraman_z_diff_re(t0:t1, k, 3)))*im_unit
                !!Normalize
                norm_iso(0:t1 - t0, k) = norm_iso(0:t1 - t0, k) + 1.0_dp
            END DO
        END DO

        !! Normalize
        rams%RR%z_iso_resraman(:, :) = rams%RR%z_iso_resraman(:, :)/REAL(norm_iso(:, :), kind=dp)
        rams%RR%z_iso_resraman(:, :) = rams%RR%z_iso_resraman(:, :)/9._dp

        !!Unit conversion of Debye^2/(E^2*fs^2) into Debye^2/(E^2*s^2)
        rams%RR%z_iso_resraman(:, :) = rams%RR%z_iso_resraman(:, :)/(fs2s*fs2s)

        !!!Hann window function
        DO i = 0, md%t_cor - 1
            DO j = 1, Nw
                rams%RR%z_iso_resraman(i, j) = rams%RR%z_iso_resraman(i, j)*((COS(i/(md%t_cor - 1.0_dp)/2.0_dp*pi))**2)
            END DO
        END DO

        !!Data mirroring
        DO i = 1, Nw
            rams%RR%z_iso_resraman(md%t_cor, i) = 0.0_dp
        END DO
        DO i = 1, md%t_cor - 1
            DO j = 1, Nw
                rams%RR%z_iso_resraman(md%t_cor + i, j) = CONJG(rams%RR%z_iso_resraman(md%t_cor - i, j))
            END DO
        END DO

        DEALLOCATE (norm_iso)

        !! Start complex autocorrelations for the anisotropic part

        DO t0 = 2, sys%framecount
            t1 = MIN(sys%framecount, t0 + md%t_cor)
            DO k = 1, Nw
                !!Real * Real
                !! (alpha_xx - alpha_yy)
                rams%RR%z_aniso_resraman(0:t1 - t0, k) = rams%RR%z_aniso_resraman(0:t1 - t0, k) + (rams%RR%alpha_resraman_x_diff_re(t0, k, 1) &
                                                                                                   - rams%RR%alpha_resraman_y_diff_re(t0, k, 2)) &
                                                         *(rams%RR%alpha_resraman_x_diff_re(t0:t1, k, 1) &
                                                           - rams%RR%alpha_resraman_y_diff_re(t0:t1, k, 2))/2.0_dp

               !! (alpha_yy - alpha_zz)
                rams%RR%z_aniso_resraman(0:t1 - t0, k) = rams%RR%z_aniso_resraman(0:t1 - t0, k) + (rams%RR%alpha_resraman_y_diff_re(t0, k, 2) &
                                                                                                   - rams%RR%alpha_resraman_z_diff_re(t0, k, 3)) &
                                                         *(rams%RR%alpha_resraman_y_diff_re(t0:t1, k, 2) &
                                                           - rams%RR%alpha_resraman_z_diff_re(t0:t1, k, 3))/2.0_dp
               !! (alpha_zz - alpha_xx)
                rams%RR%z_aniso_resraman(0:t1 - t0, k) = rams%RR%z_aniso_resraman(0:t1 - t0, k) + (rams%RR%alpha_resraman_z_diff_re(t0, k, 3) &
                                                                                                   - rams%RR%alpha_resraman_x_diff_re(t0, k, 1)) &
                                                         *(rams%RR%alpha_resraman_z_diff_re(t0:t1, k, 3) &
                                                           - rams%RR%alpha_resraman_x_diff_re(t0:t1, k, 1))/2.0_dp
               !! (alpha_xy + alpha yx)/2
                rams%RR%z_aniso_resraman(0:t1 - t0, k) = rams%RR%z_aniso_resraman(0:t1 - t0, k) &
                                                         + (rams%RR%alpha_resraman_x_diff_re(t0, k, 2)*0.50_dp + &
                                                            rams%RR%alpha_resraman_y_diff_re(t0, k, 1)*0.50_dp) &
                                                         *(rams%RR%alpha_resraman_x_diff_re(t0:t1, k, 2)*0.50_dp + &
                                                           rams%RR%alpha_resraman_y_diff_re(t0:t1, k, 1)*0.50_dp)*3.0_dp

               !! (alpha_yz + alpha_zy)/2
                rams%RR%z_aniso_resraman(0:t1 - t0, k) = rams%RR%z_aniso_resraman(0:t1 - t0, k) &
                                                         + (rams%RR%alpha_resraman_y_diff_re(t0, k, 3)*0.50_dp + &
                                                            rams%RR%alpha_resraman_z_diff_re(t0, k, 2)*0.50_dp) &
                                                         *(rams%RR%alpha_resraman_y_diff_re(t0:t1, k, 3)*0.50_dp + &
                                                           rams%RR%alpha_resraman_z_diff_re(t0:t1, k, 2)*0.50_dp)*3.0_dp

               !! (alpha_zx + alpha_xz)/2
                rams%RR%z_aniso_resraman(0:t1 - t0, k) = rams%RR%z_aniso_resraman(0:t1 - t0, k) &
                                                         + (rams%RR%alpha_resraman_z_diff_re(t0, k, 1)*0.50_dp + &
                                                            rams%RR%alpha_resraman_x_diff_re(t0, k, 3)*0.50_dp) &
                                                         *(rams%RR%alpha_resraman_z_diff_re(t0:t1, k, 1)*0.50_dp + &
                                                           rams%RR%alpha_resraman_x_diff_re(t0:t1, k, 3)*0.50_dp)*3.0_dp

                !!Imaginary * Imaginary
                !! (alpha_xx - alpha_yy)
                rams%RR%z_aniso_resraman(0:t1 - t0, k) = rams%RR%z_aniso_resraman(0:t1 - t0, k) &
                                                         + (rams%RR%alpha_resraman_x_diff_im(t0, k, 1) - rams%RR%alpha_resraman_y_diff_im(t0, k, 2)) &
                                                         *(rams%RR%alpha_resraman_x_diff_im(t0:t1, k, 1) &
                                                           - rams%RR%alpha_resraman_y_diff_im(t0:t1, k, 2))/2.0_dp
               !! (alpha_yy - alpha_zz)
                rams%RR%z_aniso_resraman(0:t1 - t0, k) = rams%RR%z_aniso_resraman(0:t1 - t0, k) &
                                                         + (rams%RR%alpha_resraman_y_diff_im(t0, k, 2) - rams%RR%alpha_resraman_z_diff_im(t0, k, 3)) &
                                                         *(rams%RR%alpha_resraman_y_diff_im(t0:t1, k, 2) &
                                                           - rams%RR%alpha_resraman_z_diff_im(t0:t1, k, 3))/2.0_dp
               !! (alpha_zz - alpha_xx)
                rams%RR%z_aniso_resraman(0:t1 - t0, k) = rams%RR%z_aniso_resraman(0:t1 - t0, k) &
                                                         + (rams%RR%alpha_resraman_z_diff_im(t0, k, 3) - rams%RR%alpha_resraman_x_diff_im(t0, k, 1)) &
                                                         *(rams%RR%alpha_resraman_z_diff_im(t0:t1, k, 3) &
                                                           - rams%RR%alpha_resraman_x_diff_im(t0:t1, k, 1))/2.0_dp

               !! (alpha_xy + alpha yx)/2
                rams%RR%z_aniso_resraman(0:t1 - t0, k) = rams%RR%z_aniso_resraman(0:t1 - t0, k) &
                                                         + (rams%RR%alpha_resraman_x_diff_im(t0, k, 2)*0.50_dp + &
                                                            rams%RR%alpha_resraman_y_diff_im(t0, k, 1)*0.50_dp) &
                                                         *(rams%RR%alpha_resraman_x_diff_im(t0:t1, k, 2)*0.50_dp + &
                                                           rams%RR%alpha_resraman_y_diff_im(t0:t1, k, 1)*0.50_dp)*3.0_dp

               !! (alpha_yz + alpha_zy)/2
                rams%RR%z_aniso_resraman(0:t1 - t0, k) = rams%RR%z_aniso_resraman(0:t1 - t0, k) &
                                                         + (rams%RR%alpha_resraman_y_diff_im(t0, k, 3)*0.50_dp + &
                                                            rams%RR%alpha_resraman_z_diff_im(t0, k, 2)*0.50_dp) &
                                                         *(rams%RR%alpha_resraman_y_diff_im(t0:t1, k, 3)*0.50_dp + &
                                                           rams%RR%alpha_resraman_z_diff_im(t0:t1, k, 2)*0.50_dp)*3.0_dp

               !! (alpha_zx + alpha_xz)/2
                rams%RR%z_aniso_resraman(0:t1 - t0, k) = rams%RR%z_aniso_resraman(0:t1 - t0, k) &
                                                         + (rams%RR%alpha_resraman_z_diff_im(t0, k, 1)*0.50_dp + &
                                                            rams%RR%alpha_resraman_x_diff_im(t0, k, 3)*0.50_dp) &
                                                         *(rams%RR%alpha_resraman_z_diff_im(t0:t1, k, 1)*0.50_dp + &
                                                           rams%RR%alpha_resraman_x_diff_im(t0:t1, k, 3)*0.50_dp)*3.0_dp

                !!Real * Imaginary
                !! (alpha_xx - alpha_yy)
                rams%RR%z_aniso_resraman(0:t1 - t0, k) = rams%RR%z_aniso_resraman(0:t1 - t0, k) &
                                                         + ((rams%RR%alpha_resraman_x_diff_re(t0, k, 1) - rams%RR%alpha_resraman_y_diff_re(t0, k, 2)) &
                                                            *(rams%RR%alpha_resraman_x_diff_im(t0:t1, k, 1) &
                                                              - rams%RR%alpha_resraman_y_diff_im(t0:t1, k, 2)))/2.0_dp*im_unit
               !! (alpha_yy - alpha_zz)
                rams%RR%z_aniso_resraman(0:t1 - t0, k) = rams%RR%z_aniso_resraman(0:t1 - t0, k) &
                                                         + ((rams%RR%alpha_resraman_y_diff_re(t0, k, 2) - rams%RR%alpha_resraman_z_diff_re(t0, k, 3)) &
                                                            *(rams%RR%alpha_resraman_y_diff_im(t0:t1, k, 2) &
                                                              - rams%RR%alpha_resraman_z_diff_im(t0:t1, k, 3)))/2.0_dp*im_unit
               !! (alpha_zz - alpha_xx)
                rams%RR%z_aniso_resraman(0:t1 - t0, k) = rams%RR%z_aniso_resraman(0:t1 - t0, k) &
                                                         + ((rams%RR%alpha_resraman_z_diff_re(t0, k, 3) - rams%RR%alpha_resraman_x_diff_re(t0, k, 1)) &
                                                            *(rams%RR%alpha_resraman_z_diff_im(t0:t1, k, 3) &
                                                              - rams%RR%alpha_resraman_x_diff_im(t0:t1, k, 1)))/2.0_dp*im_unit

               !! (alpha_xy + alpha yx)/2
                rams%RR%z_aniso_resraman(0:t1 - t0, k) = rams%RR%z_aniso_resraman(0:t1 - t0, k) &
                                                         + ((rams%RR%alpha_resraman_x_diff_re(t0, k, 2)* &
                                                             0.50_dp + rams%RR%alpha_resraman_y_diff_re(t0, k, 1) &
                                                             *0.50_dp)*(rams%RR%alpha_resraman_x_diff_im(t0:t1, k, 2)*0.50_dp &
                                                                        + rams%RR%alpha_resraman_y_diff_im(t0:t1, k, 1)*0.50_dp))*3.0_dp*im_unit
               !! (alpha_yz + alpha_zy)/2
                rams%RR%z_aniso_resraman(0:t1 - t0, k) = rams%RR%z_aniso_resraman(0:t1 - t0, k) &
                                                         + ((rams%RR%alpha_resraman_y_diff_re(t0, k, 3)* &
                                                             0.50_dp + rams%RR%alpha_resraman_z_diff_re(t0, k, 2) &
                                                             *0.50_dp)*(rams%RR%alpha_resraman_y_diff_im(t0:t1, k, 3)*0.50_dp &
                                                                        + rams%RR%alpha_resraman_z_diff_im(t0:t1, k, 2)*0.50_dp))*3.0_dp*im_unit
               !! (alpha_zx + alpha_xz)/2
                rams%RR%z_aniso_resraman(0:t1 - t0, k) = rams%RR%z_aniso_resraman(0:t1 - t0, k) &
                                                         + ((rams%RR%alpha_resraman_z_diff_re(t0, k, 1)* &
                                                             0.50_dp + rams%RR%alpha_resraman_x_diff_re(t0, k, 3) &
                                                             *0.50_dp)*(rams%RR%alpha_resraman_z_diff_im(t0:t1, k, 1)*0.50_dp &
                                                                        + rams%RR%alpha_resraman_x_diff_im(t0:t1, k, 3)*0.50_dp))*3.0_dp*im_unit

                !!Real * Imaginary (SUBTRACT)
                !! (alpha_xx - alpha_yy)
                rams%RR%z_aniso_resraman(0:t1 - t0, k) = rams%RR%z_aniso_resraman(0:t1 - t0, k) &
                                                         - ((rams%RR%alpha_resraman_x_diff_im(t0, k, 1) - rams%RR%alpha_resraman_y_diff_im(t0, k, 2)) &
                                                            *(rams%RR%alpha_resraman_x_diff_re(t0:t1, k, 1) &
                                                              - rams%RR%alpha_resraman_y_diff_re(t0:t1, k, 2)))/2.0_dp*im_unit
               !! (alpha_yy - alpha_zz)
                rams%RR%z_aniso_resraman(0:t1 - t0, k) = rams%RR%z_aniso_resraman(0:t1 - t0, k) &
                                                         - ((rams%RR%alpha_resraman_y_diff_im(t0, k, 2) - rams%RR%alpha_resraman_z_diff_im(t0, k, 3)) &
                                                            *(rams%RR%alpha_resraman_y_diff_re(t0:t1, k, 2) &
                                                              - rams%RR%alpha_resraman_z_diff_re(t0:t1, k, 3)))/2.0_dp*im_unit
               !! (alpha_zz - alpha_xx)
                rams%RR%z_aniso_resraman(0:t1 - t0, k) = rams%RR%z_aniso_resraman(0:t1 - t0, k) &
                                                         - ((rams%RR%alpha_resraman_z_diff_im(t0, k, 3) - rams%RR%alpha_resraman_x_diff_im(t0, k, 1)) &
                                                            *(rams%RR%alpha_resraman_z_diff_re(t0:t1, k, 3) &
                                                              - rams%RR%alpha_resraman_x_diff_re(t0:t1, k, 1)))/2.0_dp*im_unit
               !! (alpha_xy + alpha yx)/2
                rams%RR%z_aniso_resraman(0:t1 - t0, k) = rams%RR%z_aniso_resraman(0:t1 - t0, k) &
                                                         - ((rams%RR%alpha_resraman_x_diff_im(t0, k, 2)*0.50_dp &
                                                             + rams%RR%alpha_resraman_y_diff_im(t0, k, 1) &
                                                             *0.50_dp)*(rams%RR%alpha_resraman_x_diff_re(t0:t1, k, 2)*0.50_dp &
                                                                        + rams%RR%alpha_resraman_y_diff_re(t0:t1, k, 1)*0.50_dp))*3.0_dp*im_unit
               !! (alpha_yz + alpha_zy)/2
                rams%RR%z_aniso_resraman(0:t1 - t0, k) = rams%RR%z_aniso_resraman(0:t1 - t0, k) &
                                                         - ((rams%RR%alpha_resraman_y_diff_im(t0, k, 3)*0.50_dp &
                                                             + rams%RR%alpha_resraman_z_diff_im(t0, k, 2) &
                                                             *0.50_dp)*(rams%RR%alpha_resraman_y_diff_re(t0:t1, k, 3)*0.50_dp &
                                                                        + rams%RR%alpha_resraman_z_diff_re(t0:t1, k, 2)*0.50_dp))*3.0_dp*im_unit
               !! (alpha_zx + alpha_xz)/2
                rams%RR%z_aniso_resraman(0:t1 - t0, k) = rams%RR%z_aniso_resraman(0:t1 - t0, k) &
                                                         - ((rams%RR%alpha_resraman_z_diff_im(t0, k, 1)*0.50_dp &
                                                             + rams%RR%alpha_resraman_x_diff_im(t0, k, 3)*0.50_dp) &
                                                            *(rams%RR%alpha_resraman_z_diff_re(t0:t1, k, 1)*0.50_dp &
                                                              + rams%RR%alpha_resraman_x_diff_re(t0:t1, k, 3)*0.50_dp))*3.0_dp*im_unit
                !! Normalize
                norm_aniso(0:t1 - t0, k) = norm_aniso(0:t1 - t0, k) + 1.0_dp
            END DO
        END DO

        !! Normalize
        rams%RR%z_aniso_resraman(:, :) = rams%RR%z_aniso_resraman(:, :)/REAL(norm_aniso(:, :), kind=dp)

        !!Unit conversion of Debye^2/(E^2*fs^2) into Debye^2/(E^2*s^2)
        rams%RR%z_aniso_resraman(:, :) = rams%RR%z_aniso_resraman(:, :)/(fs2s*fs2s)

        !!!Hann window function
        DO i = 0, md%t_cor - 1
            DO j = 1, Nw
                rams%RR%z_aniso_resraman(i, j) = rams%RR%z_aniso_resraman(i, j)*0.5_dp*(1 + COS(2.0_dp*pi*i/(2.0_dp*(md%t_cor - 1))))
            END DO
        END DO

        !!Data mirroing
        DO i = 1, Nw
            rams%RR%z_aniso_resraman(md%t_cor, i) = 0.0_dp
        END DO
        DO i = 1, md%t_cor - 1
            DO j = 1, Nw
                rams%RR%z_aniso_resraman(md%t_cor + i, j) = CONJG(rams%RR%z_aniso_resraman(md%t_cor - i, j))
            END DO
        END DO

    END SUBROUTINE cvv_resraman

END MODULE vel_cor
