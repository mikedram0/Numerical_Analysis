!  Bisection method to locate roots of 
!
!  f(x)=sqrt(x)-cos(x) in [0,1].
!
!  Exercise 1b, chapter 2.

PROGRAM bisection
  IMPLICIT NONE
  DOUBLE PRECISION :: a,b    ! left and right limits bracketing the root 
  DOUBLE PRECISION :: x      ! approximation for the root
  DOUBLE PRECISION :: fnval  ! function value
  DOUBLE PRECISION :: fa, fb  ! function values

  INTERFACE 
     FUNCTION f(x)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: x
       DOUBLE PRECISION :: f
     END FUNCTION f

       FUNCTION samesign(a,b)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: a,b
       LOGICAL :: samesign
     END FUNCTION samesign
  END INTERFACE

  DOUBLE PRECISION, PARAMETER :: toler = 1D-8   ! a small constant number
  
  !  initial limits
  a=0.0D0
  b=1.0D0

  fa = f(a)
  fb = f(b)

  IF (samesign(fa, fb)) STOP

  DO 
     x=(a+b)/2.0D0
     
     fnval = f(x)
     
     PRINT *, "The current approximation is ", x 
     PRINT *, "The current value of F(x) is ", fnval
     PRINT *
     
     !    Check if root is found
     IF (ABS(fnval) < toler) EXIT

     IF (samesign(fnval,fa)) THEN 
        a = x
        fa = fnval
     ELSE
        b = x
        fb = fnval
     ENDIF
  END DO

END PROGRAM bisection


FUNCTION f(x)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: x
  DOUBLE PRECISION :: f

  f = SQRT(x)-COS(x)
END FUNCTION f



FUNCTION samesign(a,b)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: a,b
  LOGICAL :: samesign

  !     samesign = a*b > 0.0D0
  samesign = SIGN(a,b) == a
END FUNCTION samesign
