SUBROUTINE SECANT(X1, X2, X, TOLER, FUNC)
  IMPLICIT NONE
  INTEGER, PARAMETER :: prc = KIND(1d0)

  COMPLEX (prc), INTENT (inout) :: X1, X2
  COMPLEX (prc), INTENT (out)   :: X
  DOUBLE PRECISION, INTENT (in)    :: TOLER
  INTERFACE
     FUNCTION FUNC(X)
       IMPLICIT NONE
       INTEGER, PARAMETER :: prc = KIND(1d0)
       complex (prc), INTENT (in) :: X
       COMPLEX (prc) :: FUNC
     END FUNCTION FUNC
  END INTERFACE

  COMPLEX (prc) :: F1, F2

  F1 = FUNC(X1)
  F2 = FUNC(X2)

  DO 
     X = X2 - F2 * (X2 - X1) / (F2 - F1)

     X1 = X2
     F1 = F2

     X2 = X
     F2 = FUNC(X)

     IF (ABS(F2) < TOLER) EXIT
  END DO

END SUBROUTINE SECANT


FUNCTION g(x,y)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: x,y
  DOUBLE PRECISION :: g

  g =  4.0d0 * x*x - y**3 + 28.0d0
END FUNCTION g

FUNCTION h(x,y)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: x,y
  DOUBLE PRECISION :: h

  h = 3.0d0 * x**3 + 4.0d0 * y*y - 145.0d0
END FUNCTION h


FUNCTION F(Z)
  IMPLICIT NONE
  INTEGER, PARAMETER :: prc = KIND(1d0)
  COMPLEX (prc), INTENT (in) :: Z
  COMPLEX (prc) :: F

  DOUBLE PRECISION :: x,y

  INTERFACE
     FUNCTION g(x,y)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: x,y
       DOUBLE PRECISION :: g
     END FUNCTION g

     FUNCTION h(x,y)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: x,y
       DOUBLE PRECISION :: h
     END FUNCTION h

  END INTERFACE
  x = REAL(z,prc)
  y = AIMAG(z)

  f = CMPLX(g(x,y), h(x,y), prc)
END FUNCTION F

PROGRAM nonlin
  IMPLICIT NONE
  INTEGER, PARAMETER :: prc = KIND(1d0)

  INTERFACE
     FUNCTION F(Z)
       IMPLICIT NONE
       INTEGER, PARAMETER :: prc = KIND(1d0)
       COMPLEX (prc), INTENT (in) :: Z
       COMPLEX (prc) :: F
     END FUNCTION F

     SUBROUTINE SECANT(X1, X2, X, TOLER, FUNC)
       IMPLICIT NONE
       INTEGER, PARAMETER :: prc = KIND(1d0)

       COMPLEX (prc), INTENT (inout) :: X1, X2
       COMPLEX (prc), INTENT (out)   :: X
       DOUBLE PRECISION, INTENT (in)    :: TOLER
       INTERFACE
          FUNCTION FUNC(X)
            IMPLICIT NONE
            INTEGER, PARAMETER :: prc = KIND(1d0)
            complex (prc), INTENT (in) :: X
            COMPLEX (prc) :: FUNC
          END FUNCTION FUNC
       END INTERFACE
     END SUBROUTINE SECANT

  END INTERFACE

  COMPLEX(prc) :: z1,z2, z
  z1 = (0d0, 0d0)
  z2 = (1d0, 1d0)

  CALL secant(z1,z2,z,1d-8, f)

  PRINT *, REAL(z,prc), AIMAG(z)     
END PROGRAM nonlin
