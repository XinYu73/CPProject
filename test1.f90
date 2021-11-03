program name
    use MyModule
    implicit none
    REAL(4) :: x1,y1,z1
    WRITE(*,*)"Enter the position"
    READ(*,*)x1,y1,z1
    CALL GetPosition(x1,y1,z1)
    CALL spline()
end program name