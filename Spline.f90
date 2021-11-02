module MyModule
    implicit none
    save
    REAL,PARAMETER:: pi=3.14159265359
    REAL,PARAMETER:: Mz =1 !1A m**2
    REAL,PARAMETER:: omega =2 * pi * 20000
    COMPLEX,PARAMETER :: i=(0,1)
    real,PARAMETER::epsilon=8.854*10**(-12.0)
    complex,PARAMETER:: sigma=0.1-i * omega *epsilon
    REAL,PARAMETER:: mu = 4 * pi * 10**(-7.0)
    !!!!
    contains
    !test Function
    complex function cubic(arg,j) result(retval)
    implicit none
    INTEGER,INTENT(IN):: j
    complex,INTENT(IN) :: arg
    complex::kz
    real:: rho
    real::x,y,z
    x=0.5
    y=0.0
    z=1.0
    SELECT CASE (j)
    case (0)
        retval = arg**3+1.0
    case (1) ! E phi

        rho = sqrt(x**2+y**2)
        kz=sqrt(i*omega*mu*sigma-arg**2)
        retval = -(omega *mu*Mz/(4.0*pi))*(BESSEL_JN(1,real(arg*rho))*arg**2*exp(i*abs(z)*kz))/kz
    case DEFAULT
        WRITE(*,*)"Wrong day" 
    end select
    end function cubic
    !spline and others
    subroutine spline()
        implicit none
        external sgesv
        INTEGER ,PARAMETER :: N = 5000
        complex,DIMENSION(N,4)::Cmatrix
        INTEGER :: iter
        integer :: v(N+1),iflag
        complex::Integate=0
        complex,DIMENSION(N+1) :: Dd = 0.0
        complex,DIMENSION(N+1) :: Xd = 0.0
        complex,DIMENSION(N+1,N+1) :: Mmatrix = 0.0
        !
        complex::temp
        real,DIMENSION(5)::t=[-0.9061798459,0.9061798459,0.5384693101,-0.5384693101,0.0]
        real,DIMENSION(5)::W=[0.2369268851,0.2369268851,0.4786286705,0.4786286705,0.568888889]
        INTEGER :: initer
        !
        real::hj
        real::a,b ! start and end point of the intval
        INTEGER::caseSelect=1
        WRITE(*,*)"Enter the interval a,b ,b should be greater than a"
        READ(*,*)a,b
        hj=(b-a)/real(N) ! 均匀采点
        do iter = 2 , N
            Xd(iter) = real(iter)*(b-a)/ real(N)
            !WRITE(*,*)Xd(iter)
        end do
        Xd(1)=a
        Xd(N+1)=b  
        !Initializing Matrix for Mj
        Mmatrix(1,1)=2.0
        Mmatrix(1,2) = 1.0
        Mmatrix(N+1,N)=1.0
        Mmatrix(N+1,N+1) = 2.0
        do iter = 2, N
            Mmatrix(iter,iter-1) = 0.5
            Mmatrix(iter,iter) = 2
            Mmatrix(iter,iter+1) = 0.5
        end do
        !initializing done
        !initializing vector for Mj
        do iter = 2, N
            Dd(iter) = 6.0*(cubic(Xd(iter+1),caseSelect)-cubic(Xd(iter),caseSelect))/(2.0*hj) &
                        -6.0*(cubic(Xd(iter),caseSelect)-cubic(Xd(iter-1),caseSelect))/(2.0*hj)
        end do
        Dd(1)=Dd(2)
        Dd(N+1)= Dd(N)
        !WRITE(*,*)Xd
        !WRITE(*,*)Dd
        !WRITE(*,*)Mmatrix(1,:)
        call sgesv(N+1,1,Mmatrix,N+1,v,Dd,N+1,iflag)
        !WRITE(*,*)Dd   ! Dd 现在是spline所需的插值
        !********
        do iter  = 1,N
            Cmatrix(iter,4)=(Dd(iter+1)-Dd(iter))/(6.0*hj)
            Cmatrix(iter,3)=(Xd(iter+1)*Dd(iter)-Xd(iter)*Dd(iter+1))/(2.0*hj)
            Cmatrix(iter,2)=(hj/6.0)*(Dd(iter)-Dd(iter+1))+(-Xd(iter+1)**2*Dd(iter)&
                            +Xd(iter)**2*Dd(iter+1))/(2.0*hj)+(cubic(Xd(iter+1),caseSelect)-cubic(Xd(iter),caseSelect))/(hj)!cubic
            Cmatrix(iter,1)=(hj/6.0)*(Xd(iter)*Dd(iter+1)-Xd(iter+1)*Dd(iter))+&
                            (Xd(iter+1)**3 * Dd(iter)-Xd(iter)**3 * Dd(iter+1))/(6.0*hj)&
                            -(Xd(iter)*cubic(Xd(iter+1),caseSelect)-Xd(iter+1)*cubic(Xd(iter),caseSelect))/(hj)
        end do
        !WRITE(*,*)Cmatrix(N,:)
        !WRITE(*,*)1.0/real(N)
        ! the integrate
        do iter = 1,N
            Integate=Integate+(Cmatrix(iter,4)*(Xd(iter+1)**4-Xd(iter)**4))/4.0+(Cmatrix(iter,3)*(Xd(iter+1)**3-Xd(iter)**3))/3.0 &
                    +(Cmatrix(iter,2)*(Xd(iter+1)**2-Xd(iter)**2))/2.0+(Cmatrix(iter,1)*(Xd(iter+1)-Xd(iter)))
        end do
        WRITE(*,*) "Result:     ",Integate
    
        !GaussLegendre
        Integate=0
        do iter = 1,N
            do initer = 1,5
                temp=(Xd(iter+1)+Xd(iter))/2.0+(Xd(iter+1)-Xd(iter))*t(initer)/2.0
                Integate=Integate+((Xd(iter+1)-Xd(iter))/2.0)*W(initer)*((Cmatrix(iter,4)*(temp)**3)&
                        +(Cmatrix(iter,3)*(temp)**2) &
                        +(Cmatrix(iter,2)*(temp)**1) &
                        +(Cmatrix(iter,1)))
            end do
        end do
        WRITE(*,*) "Result:     ",Integate
    end subroutine spline
end module MyModule