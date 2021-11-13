module AlgebraMethod
    implicit none
    save
    REAL,PARAMETER:: pi=3.14159265359
    REAL,PARAMETER:: Mz =1 !1A m**2
    REAL,PARAMETER:: omega =2 * pi * 20000000
    COMPLEX,PARAMETER :: i=(0,1)
    real,PARAMETER::epsilon=8.854*10**(-12.0)
    complex,PARAMETER:: sigma=0.1-i * omega *epsilon
    REAL,PARAMETER:: mu = 4 * pi * 10**(-7.0)
    COMPLEX :: k=sqrt(i*omega*mu*sigma)
    real::x,y,z
    !!!!!!!!!!
    contains
    subroutine GetPosition(arg1, arg2 , arg3)
        implicit none
        real,intent(in) :: arg1
        real,intent(in) ::  arg2
        real,intent(in) ::  arg3
        x=arg1
        y=arg2
        z=arg3
    end subroutine GetPosition
    subroutine GetEx(Ex)
    implicit none
    REAL :: r 
    COMPLEX,INTENT(OUT)::Ex
    r=sqrt( x**2+y**2+z**2)
    Ex = ( Mz / (4 * pi)) * i *  omega *  mu * (y * exp(i *  k *  r) /  r**3) * (i *  k *  r - 1)
    end subroutine GetEx
    !!!!!!
    subroutine GetEy(Ey)
        implicit none
        REAL :: r 
        COMPLEX,INTENT(OUT)::Ey
        r=sqrt( x**2+y**2+z**2)
        Ey = -( Mz / (4 * pi)) * i *  omega *  mu * ( x * exp(i *  k *  r) /  r**3) * (i *  k *  r - 1)
    end subroutine GetEy
    !!!!!
    subroutine GetEz(Ez)
        COMPLEX,INTENT(OUT)::Ez
        Ez=0*x*y*z
    end subroutine GetEz
    !!!
    subroutine GetH(Hx,Hy,Hz)
        REAL :: r 
        COMPLEX,INTENT(OUT)::Hx,Hy,Hz
        r=sqrt( x**2+y**2+z**2)
        Hx=-((Mz*x*z*exp(i*k*r))/(4*pi))*(k**2/r**3+3*i*k/r**4-3/r**5)
        Hy=-((Mz*y*z*exp(i*k*r))/(4*pi))*(k**2/r**3+3*i*k/r**4-3/r**5)
        Hz=((Mz*exp(i*k*r))/(4*pi*r))*(k**2+i*k/r-(k**2*z**2+1)/r**2-(3*i*k*z**2)/r**3+3*z**2/r**4)
    end subroutine GetH
end module AlgebraMethod




