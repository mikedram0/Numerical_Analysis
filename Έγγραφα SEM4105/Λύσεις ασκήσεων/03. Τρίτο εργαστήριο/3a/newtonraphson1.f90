!     Epilisi tis f(x)=0 me f(x)=sin(x)-x**2
!     me th methodo Newton-Raphson.

PROGRAM newton_raphson
  IMPLICIT NONE
  DOUBLE PRECISION ::  X

  INTERFACE 
     FUNCTION F(X)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: X
       DOUBLE PRECISION :: F
     END FUNCTION F

     FUNCTION DF(X)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: X
       DOUBLE PRECISION :: DF
     END FUNCTION DF
  END INTERFACE

  DOUBLE PRECISION, PARAMETER :: TOLER = 1D-8  ! A small constant number

  X = 1.0D0    !     Initial approximation 
  DO
     X = X - F(X) / DF(X)
     IF ( ABS(F(X)) < TOLER ) EXIT
  ENDDO

  PRINT *, "The root is approximately ", X
  PRINT *, "The function value is ", F(X)

END PROGRAM newton_raphson

!     The function
FUNCTION F(X)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: X
  DOUBLE PRECISION :: F

  F = SIN(X)-X**2
END FUNCTION F

!     The first derivative 
FUNCTION DF(X)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: X
  DOUBLE PRECISION :: DF

  DF = COS(X) - 2.0D0 * X
END FUNCTION DF
