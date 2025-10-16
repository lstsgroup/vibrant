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

!> @brief Module containing procedures involving the formatting and outputting the results.
!>
MODULE output_io

    USE kinds, ONLY: dp, str_len
    USE vib_types, ONLY: global_settings, systems, static, dipoles, raman
    USE constants, ONLY: bohr2ang
    USE ISO_FORTRAN_ENV, ONLY: output_unit, error_unit

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: write_spectra_data, write_mol_file, append_column, check_file_open
CONTAINS

    !> @brief Checks the status of a file open operation and aborts on failure.
    !>
    !> @param[in] stat      -- I/O status flag returned by an OPEN statement
    !> @param[in] msg       -- Custom error message or context description
    !> @param[in] filename  -- Name of the file that failed to open
    !>
    SUBROUTINE check_file_open(stat, msg, filename)

        INTEGER, INTENT(IN) :: stat
        CHARACTER(*), INTENT(IN) :: msg
        CHARACTER(*), INTENT(IN) :: filename

        IF (stat/=0) THEN
            WRITE (error_unit, '(4X,"[ERROR] could not open file ",A)') TRIM(filename)
            WRITE (error_unit, '(4X,"I/O error message: ",A)') TRIM(msg)
            STOP
        END IF

    END SUBROUTINE check_file_open

!***************************************************************************************!
!***************************************************************************************!

    !> @brief Writes spectral data (frequency and intensity) to a text file.
    !>
    !> @param[in] outfile               --  Name of the output file
    !> @param[in] header                --  Header label for the intensity column
    !> @param[in] freq                  --  Array of frequencies
    !> @param[in] intensities           --  Array of intensities corresponding to each
    !>                                      frequency
    !> @param[in, optional] freq_cutoff --  Optional upper frequency limit for output
    !>
    SUBROUTINE write_spectra_data(outfile, header, freq, intensities, freq_cutoff)
        CHARACTER(*), INTENT(in) :: outfile
        CHARACTER(*), INTENT(in) :: header
        REAL(dp), INTENT(in) :: freq(:)
        REAL(dp), INTENT(in) :: intensities(:)
        REAL(dp), INTENT(in), OPTIONAL :: freq_cutoff

        INTEGER :: runit, stat, i
        CHARACTER(len=256) :: msg

        OPEN (file=outfile, status='replace', action='write', iostat=stat, iomsg=msg, newunit=runit)
        CALL check_file_open(stat, msg, outfile)

        WRITE (runit, '(A22,1X,A25)') 'FREQ', TRIM(header)
        !WRITE (runit, '(A6,6X,A20)') 'FREQ', TRIM(header)
        IF (PRESENT(freq_cutoff)) THEN

            DO i = 1, SIZE(freq)
                IF (freq(i).GE.freq_cutoff) EXIT
                WRITE (runit, '(F22.16,1X,ES25.16E3)') freq(i), intensities(i)
            END DO
        ELSE
            DO i = 1, SIZE(freq)
                WRITE (runit, '(I22,1X,ES25.16E3)') i, intensities(i)
            END DO
        END IF

        CLOSE (runit)
    END SUBROUTINE write_spectra_data

!***************************************************************************************!
!***************************************************************************************!

    !> @brief Writes a Molden-compatible .mol file containing molecular geometry,
    !>        vibrational frequencies, and optionally intensities.
    !>
    !> @param[inout] sys         -- System data (atomic elements and coordinates)
    !> @param[inout] stats       -- Vibrational analysis data (frequencies, displacements)
    !> @param[inout] gs          -- Global settings (determines spectral type)
    !> @param[in]    outfile     -- Name of the output .mol file
    !> @param[inout] intensities -- Optional array of intensities corresponding to each
    !>                              vibrational mode
    !>
    SUBROUTINE write_mol_file(sys, stats, gs, outfile, intensities)
        TYPE(systems), INTENT(INOUT)        :: sys
        TYPE(static), INTENT(INOUT)        :: stats
        TYPE(global_settings), INTENT(INOUT)        :: gs
        CHARACTER(*), INTENT(IN) :: outfile
        REAL(dp), INTENT(INOUT), OPTIONAL :: intensities(:)

        INTEGER :: runit, stat, i, j
        CHARACTER(len=256) :: msg

        OPEN (FILE=outfile, STATUS='unknown', ACTION='write', IOSTAT=stat, IOMSG=msg, NEWUNIT=runit)
        CALL check_file_open(stat, msg, outfile)
        WRITE (runit, *) "[Molden Format]"
        WRITE (runit, *) "[GEOMETRIES] XYZ"
        WRITE (runit, *) sys%natom
        WRITE (runit, *)
        DO i = 1, sys%natom
            WRITE (runit, *) sys%element(i), sys%coord(i, 1), sys%coord(i, 2), sys%coord(i, 3)
        END DO
        WRITE (runit, *) "[FREQ]"
        DO i = 1, stats%nmodes
            WRITE (runit, *) stats%freq(i)
        END DO
        IF (gs%spectral_type%read_function.NE.'NMA') THEN
            WRITE (runit, *) "[INT]"
            intensities = REAL(intensities/MINVAL(intensities), kind=dp)
            DO i = 1, stats%nmodes
                WRITE (runit, *) intensities(i)
            END DO
        END IF
        WRITE (runit, *) "[FR-COORD]"
        WRITE (runit, *) sys%natom
        WRITE (runit, *)
        DO i = 1, sys%natom
            WRITE (runit, *) sys%element(i), sys%coord(i, 1)/bohr2ang, sys%coord(i, 2)/bohr2ang, sys%coord(i, 3)/bohr2ang
        END DO
        WRITE (runit, *) "[FR-NORM-COORD]"
        DO i = 1, stats%nmodes
            WRITE (runit, *) "vibration", i
            DO j = 1, sys%natom
                WRITE (runit, *) stats%disp(i, j, 1)/bohr2ang, stats%disp(i, j, 2)/bohr2ang, stats%disp(i, j, 3)/bohr2ang
            END DO
        END DO
        CLOSE (runit)
    END SUBROUTINE write_mol_file

!***************************************************************************************!
!***************************************************************************************!
    !> @brief Appends a new data column to an existing spectra file.
    !>
    !> This subroutine adds a column of intensity data to an existing text file
    !> containing frequency values and previously written data columns.
    !> If the file is newly created, the provided `header` is written as the
    !> column label on the first line.
    !>
    !> @param[in] filename         -- Name of the target file to append to
    !> @param[in] header           -- Header label for the new column
    !> @param[in] intensities      -- Array of intensity values to append
    !> @param[in] freq             -- Array of frequencies corresponding to the intensity values
    !> @param[in, optional] cutoff -- Optional upper frequency limit for writing
    !>
    SUBROUTINE append_column(filename, header, intensities, freq, cutoff)

        CHARACTER(len=*), INTENT(in) :: filename
        CHARACTER(len=*), INTENT(in) :: header
        REAL(kind=dp), INTENT(in) :: intensities(:)
        REAL(dp), INTENT(in) :: freq(:)
        REAL(dp), INTENT(in), OPTIONAL :: cutoff

        INTEGER :: runit, ios, n, i, stat, nlines
        LOGICAL :: exists
        CHARACTER(len=1024) :: line
        CHARACTER(len=30)   :: VALUE, msg
        CHARACTER(len=1024) :: buf
        CHARACTER(len=str_len), ALLOCATABLE :: lines(:)
        REAL(dp), ALLOCATABLE :: freq_in_file(:)

        ! ---------- loop 1: read all data rows into lines(:) ----------
        OPEN (file=filename, status="old", action="readwrite", IOMSG=msg, IOSTAT=stat, newunit=runit)
        CALL check_file_open(stat, msg, filename)

        nlines = 0
        DO
            READ (runit, '(A)', iostat=ios) buf
            IF (ios/=0) EXIT
            IF (LEN_TRIM(buf)==0) EXIT   ! stop at empty line if you want
            nlines = nlines + 1
        END DO
        REWIND (runit)

        ! --- allocate after counting ---
        ALLOCATE (lines(nlines))
        ALLOCATE (freq_in_file(nlines - 1))

        ! --- pass 2: fill  lines and extract freq data---
        DO i = 1, nlines
            READ (runit, '(A)', iostat=ios) lines(i)
            IF (i.GE.2) THEN
                READ (lines(i), *) freq_in_file(i - 1)
            END IF
        END DO
        REWIND (runit)

        ! ---------- loop 2: write the updated file ----------
        WRITE (runit, '(A,1X,A25)') TRIM(lines(1)), TRIM(header)
        IF (PRESENT(cutoff)) THEN
            DO i = 1, SIZE(lines) - 1
                IF (freq(i).GE.cutoff) EXIT
                IF (freq_in_file(i)==freq(i)) THEN
                    !PRINT *,freq(i)!, TRIM(lines(i+1))
                    WRITE (runit, '(A,1X,ES25.16E3)') TRIM(lines(i + 1)), intensities(i)
                ELSE
                    WRITE (error_unit, '(4X,"[ERROR] ",A)') "append_column: freq in file not equal to given freq."
                    RETURN
                END IF
            END DO
        ELSE
            DO i = 1, SIZE(freq)
                IF (freq_in_file(i)==freq(i)) THEN
                    !PRINT *,freq(i)!, TRIM(lines(i+1))
                    WRITE (runit, '(A,1X,ES25.16E3)') TRIM(lines(i + 1)), intensities(i)
                ELSE
                    WRITE (error_unit, '(4X,"[ERROR] ",A)') "append_column: freq in file not equal to given freq."
                    RETURN
                END IF
                !WRITE (runit, '(A,1X,ES25.16E3)') TRIM(lines(i+1)), intensities(i)
            END DO
        END IF
        CLOSE (runit)
    END SUBROUTINE append_column

END MODULE output_io
