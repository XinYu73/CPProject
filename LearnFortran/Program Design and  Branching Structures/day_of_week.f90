program dayOfWeek
    implicit none
    CHARACTER(len=11)::cDay
    CHARACTER(len=11)::cType
    WRITE(*,*)"Enter the name of the day :  "
    READ(*,*)cDay
    SELECT CASE (cDay)
    case ('Monday','Tuesday','Wednesday','Thursday','Friday')
        cType="Weekday"
    case ('Saturday','Sunday')
        cType="Weekend"
    case DEFAULT
        cType="Wrong day"
    end select
    WRITE (*,*) 'Day Type = ', cType
end program dayOfWeek