!  Fixed point method to locate roots of 
!   
!  f(x)=x-cos^3 x  near 0.6.
!
!  Exercise 4b, chapter 2.

PROGRAM fixed_point
  IMPLICIT NONE

  INTERFACE 
     FUNCTION f(x)
       IMPLICIT NONE
       DOUBLE PRECISION :: x,f
     END FUNCTION f

     FUNCTION g(x)
       IMPLICIT NONE
       DOUBLE PRECISION :: x,g
     END FUNCTION g
  END INTERFACE

  DOUBLE PRECISION :: x
  DOUBLE PRECISION :: fnval    !  function value 
  DOUBLE PRECISION, PARAMETER :: toler = 1D-8   !  a small constant number
  DOUBLE PRECISION, PARAMETER :: vfar =  1D5    !  a very large number

  X = 0.7D0 !     initial guess for x

  DO
     fnval = f(x)
     
     PRINT *, "The current approximation is ", x 
     PRINT *, "The current value of F(x) is ", fnval
     PRINT *
     
     IF (ABS(x) > vfar) THEN
        PRINT *, "X is very large; probably g(x) is not appropriate"
        EXIT
     ENDIF
     
     IF (ABS(fnval) < toler) EXIT
     
     X = G(X)
  END DO

END PROGRAM fixed_point


FUNCTION f(x)
  IMPLICIT NONE
  DOUBLE PRECISION :: x,f

  f = x-COS(x)**3
END FUNCTION f

!     if f(x) above then 
!     f(x) = 0 => x = g(x) 
!     with g as follows:

FUNCTION g(x)
  IMPLICIT NONE
  DOUBLE PRECISION :: x,g
  
  G = ACOS(X**(1.D0/3.D0))
  
  !     g(x)=cos^3(x) is not appropriate.
END FUNCTION g
