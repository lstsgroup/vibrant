MODULE output_io

    USE kinds, ONLY: dp, str_len
    USE ISO_FORTRAN_ENV, ONLY: output_unit, error_unit

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: write_spectra_data, append_column, check_file_open
    CONTAINS


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

    SUBROUTINE write_spectra_data(outfile, header, freq, intensities, freq_cutoff)
        CHARACTER(*), INTENT(in) :: outfile
        CHARACTER(*), INTENT(in) :: header
        !INTEGER, INTENT(in) :: start_freq, end_freq
        REAL(dp), INTENT(in) :: freq(:)
        REAL(dp), INTENT(in) :: intensities(:)
        REAL(dp),     INTENT(in), OPTIONAL :: freq_cutoff

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

    SUBROUTINE append_column(filename, header, intensities, freq, cutoff)
        !! Fügt eine Spalte `data` in die Datei `filename` hinzu.
        !! Falls die Datei neu ist, wird `header` als Kopfzeile geschrieben.
        CHARACTER(len=*), INTENT(in) :: filename
        CHARACTER(len=*), INTENT(in) :: header
        REAL(kind=dp), INTENT(in) :: intensities(:)
        REAL(dp),     INTENT(in) :: freq(:)
        REAL(dp),     INTENT(in), OPTIONAL :: cutoff
    
        INTEGER :: runit, ios, n, i, stat, nlines
        LOGICAL :: exists
        CHARACTER(len=1024) :: line
        CHARACTER(len=30)   :: VALUE, msg
        CHARACTER(len=1024) :: buf
        CHARACTER(len=str_len), ALLOCATABLE :: lines(:)
        REAL(dp),    ALLOCATABLE :: freq_in_file(:)
    
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
        ALLOCATE (freq_in_file(nlines-1))
    
        ! --- pass 2: fill  lines and extract freq data---
        DO i = 1, nlines
            READ (runit, '(A)', iostat=ios) lines(i)
            IF (i.GE.2) THEN
            READ (lines(i), *) freq_in_file(i-1)
            END IF
        END DO
        REWIND (runit)

        ! ---------- loop 2: write the updated file ----------
        WRITE (runit, '(A,1X,A25)') TRIM(lines(1)), TRIM(header)
        IF (PRESENT(cutoff)) THEN
            DO i = 1, SIZE(lines) -1
                IF (freq(i).GE.cutoff) EXIT
                IF (freq_in_file(i) == freq(i))THEN
                !PRINT *,freq(i)!, TRIM(lines(i+1))
                    WRITE (runit, '(A,1X,ES25.16E3)') TRIM(lines(i+1)), intensities(i)
                ELSE
                    WRITE (error_unit, '(4X,"[ERROR] ",A)') "append_column: freq in file not equal to given freq."
                    RETURN
                END IF
            END DO
        ELSE
            DO i = 1, SIZE(freq)
                IF (freq_in_file(i) == freq(i))THEN
                    !PRINT *,freq(i)!, TRIM(lines(i+1))
                        WRITE (runit, '(A,1X,ES25.16E3)') TRIM(lines(i+1)), intensities(i)
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