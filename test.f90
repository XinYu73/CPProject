program name
    use AlgebraMethod
    implicit none
    REAL :: x1,y1,z1
    COMPLEX :: Ex,Ey,Ez,Hx,Hy,Hz
    WRITE(*,*)"Enter the position"
    READ(*,*)x1,y1,z1
    CALL GetPosition(x1,y1,z1)
    CALL GetEx(Ex)
    CALL GetEy(Ey)
    CALL GetEz(Ez)
    WRITE(*,*)"E:",Ex,Ey,Ez
    CALL GetH(Hx,Hy,Hz)
    WRITE(*,*)"H:",Hx,Hy,Hz
end program name