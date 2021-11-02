program namee
    use MyModule
    use AlgebraMethod
    implicit none
    REAL :: x,y,z
    COMPLEX :: Ex,Ey,Ez,Hx,Hy,Hz
    WRITE(*,*)"Enter the position"
    READ(*,*)x,y,z
    CALL GetEx(x,y,z,Ex)
    CALL GetEy(x,y,z,Ey)
    CALL GetEz(x,y,z,Ez)
    WRITE(*,*)"Ex:  ",Ex,"Ey:   ",Ey,"Ez    ",Ez
    CALL GetH(x,y,z,Hx,Hy,Hz)
    WRITE(*,*)"H:",Hx," ",Hy," ",Hz
    CALL spline()
end program namee