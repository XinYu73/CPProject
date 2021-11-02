module AlgebraMethod
    implicit none
    save
    REAL,PARAMETER:: pi=3.14159265359
    REAL,PARAMETER:: Mz =1 !1A m**2
    REAL,PARAMETER:: omega =2 * pi * 20000
    COMPLEX,PARAMETER :: i=(0,1)
    real,PARAMETER::epsilon=8.854*10**(-12.0)
    complex,PARAMETER:: sigma=0.1-i * omega *epsilon
    REAL,PARAMETER:: mu = 4 * pi * 10**(-7.0)
    COMPLEX :: k=sqrt(i*omega*mu*sigma)
    !!!!!!!!!!
    contains
    subroutine GetEx(x,y,z,Ex)
    implicit none
    REAL,intent(in) :: x
    REAL,intent(in) :: y
    REAL,intent(in) :: z
    REAL :: r 
    COMPLEX,INTENT(OUT)::Ex
    r=sqrt( x**2+y**2+z**2)
    Ex = ( Mz / (4 * pi)) * i *  omega *  mu * (y * exp(i *  k *  r) /  r**3) * (i *  k *  r - 1)
    end subroutine GetEx
    !!!!!!
    subroutine GetEy(x,y,z,Ey)
        implicit none
        REAL,intent(in) :: x
        REAL,intent(in) :: y
        REAL,intent(in) :: z
        REAL :: r 
        COMPLEX,INTENT(OUT)::Ey
        r=sqrt( x**2+y**2+z**2)
        Ey = -( Mz / (4 * pi)) * i *  omega *  mu * ( x * exp(i *  k *  r) /  r**3) * (i *  k *  r - 1)
    end subroutine GetEy
    !!!!!
    subroutine GetEz(x,y,z,Ez)
        REAL,intent(in) :: x
        REAL,intent(in) :: y
        REAL,intent(in) :: z
        COMPLEX,INTENT(OUT)::Ez
        Ez=0*x*y*z
    end subroutine GetEz
    !!!
    subroutine GetH(x,y,z,Hx,Hy,Hz)
        REAL,intent(in) :: x
        REAL,intent(in) :: y
        REAL,intent(in) :: z
        REAL :: r 
        COMPLEX,INTENT(OUT)::Hx,Hy,Hz
        r=sqrt( x**2+y**2+z**2)
        Hx=-((Mz*x*z*exp(i*k*r))/(4*pi))*(k**2/r**3+3*i*k/r**4-3/r**5)
        Hy=-((Mz*y*z*exp(i*k*r))/(4*pi))*(k**2/r**3+3*i*k/r**4-3/r**5)
        Hz=((Mz*exp(i*k*r))/(4*pi*r))*(k**2+i*k/r-(k**2*z**2+1)/r**2-(3*i*k*z**2)/r**3+3*z**2/r**4)
    end subroutine GetH
end module AlgebraMethod




