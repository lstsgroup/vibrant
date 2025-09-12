MODULE timing

    USE kinds, ONLY: dp

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
            PRINT *, "Error: Maximum number of events reached."
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
        REAL(kind=8) :: total_time
        IF (this%n_events>0) THEN
            CALL this%events(this%n_events)%stop_event()
        END IF
        PRINT *, ""
        WRITE (*, '(A)') REPEAT('-', 80)
        PRINT *, "Timing Report:"
        total_time = 0.0_DP
        DO i = 1, this%n_events
            total_time = total_time + this%events(i)%elapsed
            CALL this%events(i)%report_event()
        END DO
        PRINT *, ""
        PRINT '(A50, A, F12.4, A)', "TOTAL", ": ", total_time, " s"
        WRITE (*, '(A)') REPEAT('-', 80)
        PRINT *, ""
    END SUBROUTINE report_all

    !> @brief start timing an event
    !> @param name -- the name of the event
    SUBROUTINE start_event(this, name)
        CLASS(event), INTENT(inout) :: this
        CHARACTER(len=*), INTENT(in) :: name
        INTEGER :: count, rate

        IF (this%running .OR. this%is_over) THEN
            PRINT *, "Warning: Event is already running or over."
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
            PRINT *, "Warning: Event is not running or over."
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
        IF (this%running .OR. .NOT. this%is_over) THEN
            PRINT *, "Event '", TRIM(this%name), "' is still running or didn't start."
        ELSE
            PRINT '(A50, A, F12.4, A)', TRIM(this%name), ": ", this%elapsed, " s"
        END IF
    END SUBROUTINE report_event

END MODULE timing
