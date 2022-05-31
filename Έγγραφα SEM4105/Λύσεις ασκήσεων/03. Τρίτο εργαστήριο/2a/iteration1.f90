!  Fixed point method to locate a root of 
!   
!   f(x)=x^2-6x+5

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

  X = 1.3D0 !     initial guess for x

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
  DOUBLE PRECISION, INTENT (in) :: x
  DOUBLE PRECISION :: f

  f = x**2 - 6.0D0 * x + 5.0D0
END FUNCTION f


!     if f(x) above then 
!     f(x) = 0 => x = g(x) 
!     with g as follows:

FUNCTION g(x)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: x
  DOUBLE PRECISION :: g

  g = (x**2 + 5.0D0) / 6.0D0
END FUNCTION g
