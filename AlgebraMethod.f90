program AlgebraMethod
    implicit none
    !data dictionary
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
end program AlgebraMethod

subroutine GetEx(x,y,z,Ex)
    implicit none
    REAL,PARAMETER:: pi=3.14159265
    REAL,PARAMETER:: Mz =1 !1A m**2
    REAL,PARAMETER:: omega =2 * pi * 20
    REAL,PARAMETER:: sigma=0.1
    REAL,PARAMETER:: mu = 4 * pi * 10**(-7.0)
    REAL,intent(in) :: x
    REAL,intent(in) :: y
    REAL,intent(in) :: z
    REAL :: r 
    COMPLEX,INTENT(OUT)::Ex
    COMPLEX,PARAMETER :: i=(0,1)
    COMPLEX :: k=sqrt(i*omega*mu*sigma)
    r=sqrt( x**2+y**2+z**2)
    Ex = ( Mz / (4 * pi)) * i *  omega *  mu * (y * exp(i *  k *  r) /  r**3) * (i *  k *  r - 1)
end subroutine GetEx

subroutine GetEy(x,y,z,Ey)
    implicit none
    REAL,PARAMETER:: pi=3.14159265
    REAL,PARAMETER:: Mz =1 !1A m**2
    REAL,PARAMETER:: omega =2 * pi * 20000
    COMPLEX:: sigma
    REAL,PARAMETER:: mu = 4 * pi * 10**(-7.0)
    REAL,intent(in) :: x
    REAL,intent(in) :: y
    REAL,intent(in) :: z
    REAL :: r 
    COMPLEX :: k
    COMPLEX,INTENT(OUT)::Ey
    COMPLEX,PARAMETER :: i=(0,1)
    WRITE(*,*)mu
    r=sqrt( x**2+y**2+z**2)
    sigma=0.1
    k=sqrt(i*omega*mu*sigma)
    Ey = -( Mz / (4 * pi)) * i *  omega *  mu * ( x * exp(i *  k *  r) /  r**3) * (i *  k *  r - 1)
end subroutine GetEy

subroutine GetEz(x,y,z,Ez)
    implicit none
    REAL,intent(in) :: x
    REAL,intent(in) :: y
    REAL,intent(in) :: z
    COMPLEX,INTENT(OUT)::Ez
    Ez=0*x*y*z
end subroutine GetEz

subroutine GetH(x,y,z,Hx,Hy,Hz)
    implicit none
    REAL,PARAMETER:: pi=3.14159265
    REAL,PARAMETER:: Mz =1 !1A m**2
    REAL,PARAMETER:: omega =2 * pi * 20000
    COMPLEX:: sigma
    REAL,PARAMETER:: mu = 4 * pi * 10**(-7.0)
    REAL,intent(in) :: x
    REAL,intent(in) :: y
    REAL,intent(in) :: z
    REAL :: r 
    COMPLEX :: k
    COMPLEX,INTENT(OUT)::Hx,Hy,Hz
    COMPLEX,PARAMETER :: i=(0,1)
    r=sqrt( x**2+y**2+z**2)
    sigma=0.1
    k=sqrt(i*omega*mu*sigma)
    Hx=-((Mz*x*z*exp(i*k*r))/(4*pi))*(k**2/r**3+3*i*k/r**4-3/r**5)
    Hy=-((Mz*y*z*exp(i*k*r))/(4*pi))*(k**2/r**3+3*i*k/r**4-3/r**5)
    Hz=((Mz*exp(i*k*r))/(4*pi*r))*(k**2+i*k/r-(k**2*z**2+1)/r**2-(3*i*k*z**2)/r**3+3*z**2/r**4)
end subroutine GetH