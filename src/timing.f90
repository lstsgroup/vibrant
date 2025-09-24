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

MODULE timing

    USE kinds, ONLY: dp, str_len
    USE iso_fortran_env, ONLY: output_unit, error_unit

    IMPLICIT NONE
    PRIVATE
    PUBLIC :: timings

    TYPE :: event
        REAL(kind=8) :: start
        REAL(kind=8) :: STOP
        REAL(kind=8) :: elapsed
        LOGICAL :: running
        LOGICAL :: is_over
        CHARACTER(len=100) :: name
    CONTAINS
        PROCEDURE :: start_event
        PROCEDURE :: stop_event
        PROCEDURE :: report_event
    END TYPE event

    TYPE :: time_table
        INTEGER :: n_events = 0
        TYPE(event), DIMENSION(100) :: events
    CONTAINS
        PROCEDURE :: register
        PROCEDURE :: report_all
    END TYPE time_table

    !> global time table instance
    TYPE(time_table) :: timings

CONTAINS

    !> @brief record the timing of an event and stop the last one
    !> @param name -- the name of the new event
    SUBROUTINE register(this, name)
        CLASS(time_table), INTENT(inout) :: this
        CHARACTER(len=*), INTENT(in) :: name
        IF (this%n_events>=SIZE(this%events)) THEN
            WRITE(error_unit,'(4X,"[ERROR] ",A)') 'Maximum number of events reached.'
            RETURN
        END IF
        IF (this%n_events>0) THEN
            CALL this%events(this%n_events)%stop_event()
        END IF
        this%n_events = this%n_events + 1
        CALL this%events(this%n_events)%start_event(name)
    END SUBROUTINE register

    !> @brief report the timing of all events
    SUBROUTINE report_all(this)
        CLASS(time_table), INTENT(inout) :: this
        INTEGER :: i
        CHARACTER(len=str_len) :: time
        REAL(kind=8) :: total_time
        IF (this%n_events>0) THEN
            CALL this%events(this%n_events)%stop_event()
        END IF
        WRITE(*,'(/,90A)') REPEAT("-",90)
        WRITE(*,'(2X, A)') "Timing Report:"
        total_time = 0.0_dp
        DO i = 1, this%n_events
            total_time = total_time + this%events(i)%elapsed
            CALL this%events(i)%report_event()
        END DO
        !WRITE(*,'(/,T17,A52, A, T60,F12.4, A)') "TOTAL", ": ", total_time, " s"
        WRITE(time,'(F12.4, " s")') total_time  
        WRITE(*,'(/,T17,A,":",T60,A)') "TOTAL", TRIM(ADJUSTL(time))
        WRITE(*,'(90A,/)') REPEAT("-",90)
    END SUBROUTINE report_all

    !> @brief start timing an event
    !> @param name -- the name of the event
    SUBROUTINE start_event(this, name)
        CLASS(event), INTENT(inout) :: this
        CHARACTER(len=*), INTENT(in) :: name
        INTEGER :: count, rate

        IF (this%running .OR. this%is_over) THEN
            WRITE(error_unit,'(4X,"[WARN]  ",A)') 'Event is already running or over.'
            RETURN
        END IF

        this%name = TRIM(name)
        CALL SYSTEM_CLOCK(count, rate)
        this%start = REAL(count, kind=8)/REAL(rate, kind=8)  ! seconds
        this%running = .TRUE.
    END SUBROUTINE start_event

    !> @brief stop timing of this event
    SUBROUTINE stop_event(this)
        CLASS(event), INTENT(inout) :: this
        INTEGER :: count, rate

        IF (.NOT. this%running .OR. this%is_over) THEN
            WRITE(error_unit,'(4X,"[WARN]  ",A)') 'Event is already running or over.'
            RETURN
        END IF

        CALL SYSTEM_CLOCK(count, rate)
        this%STOP = REAL(count, kind=8)/REAL(rate, kind=8)   ! seconds
        this%elapsed = this%STOP - this%start
        this%running = .FALSE.
        this%is_over = .TRUE.
    END SUBROUTINE stop_event

    !> @brief report the timing of this event
    SUBROUTINE report_event(this)
        CLASS(event), INTENT(in) :: this
        CHARACTER(len=str_len) :: time

        IF (this%running .OR. .NOT. this%is_over) THEN
            WRITE(error_unit,'(4X,"[INFO]  ",A,A)') TRIM(this%name), "' is still running or didn't start."
        ELSE
            WRITE(time,'(F12.4, " s")') this%elapsed    
            WRITE(*,'(T17,A,":",T60,A)') TRIM(this%name), TRIM(ADJUSTL(time))

        END IF
    END SUBROUTINE report_event

END MODULE timing
