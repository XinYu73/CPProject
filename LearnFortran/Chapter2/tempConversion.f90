program tempConversion
    implicit none
    REAL :: temp_f,temp_k
    WRITE(*,*)"Enter the temperature in degrees Fahrenheit"
    READ(*,*) temp_f
    CALL convert(temp_f,temp_k)
    WRITE(*,*)temp_f, ' degrees Fahrenheit = ', temp_k, ' kelvins'
end program tempConversion

subroutine convert(f,  k)
implicit none
REAL,intent(in) :: f
REAL,intent(out) ::  k
k=(5./9.)*(f-32.)+273.15
end subroutine convert