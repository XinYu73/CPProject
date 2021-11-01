program roots
    implicit none
    REAL::a,b,c,discriminant,imagPart,realPart,x1,x2
    WRITE(*,*)"Enter the coefs of aX^2 + BX+C =0"
    READ(*,*)a,b,c
    WRITE(*,*)"The coefs are :  ",a,b,c
    discriminant = b**2.0 - 4. * a * c
    if (discriminant>0 ) then
        x1 = (-b - sqrt(discriminant))/(2.*a)
        x2 = (-b + sqrt(discriminant))/(2.*a)
        WRITE(*,*)"There are two different roots :  ",x1,x2
    else if (discriminant<0 ) then
        realPart = ( -b ) / ( 2. * a )
        imagPart = sqrt ( abs ( discriminant ) ) / ( 2. * a )
        WRITE (*,*) 'X1 = ', realPart, ' +i ', imagPart
        WRITE (*,*) 'X2 = ', realPart, ' -i ', imagPart 
    end if
end program roots
