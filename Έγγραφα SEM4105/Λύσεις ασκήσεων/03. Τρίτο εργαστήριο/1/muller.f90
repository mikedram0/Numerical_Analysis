! Find the roots of f(x)=0 with f(x)=sin(x) - x^2
! using the Muller's method.
! use complex arithmetic

PROGRAM MULLERS
  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: epsi = 1d-6, delta = 1d-6

  COMPLEX (KIND(1d0)) :: X,X0,X1,X2, W1, W0,A,B,C, D, D1, D2, P
  COMPLEX (KIND(1d0)) :: F0,F1,F2, T

  INTERFACE
     FUNCTION F(X)
       IMPLICIT NONE
       COMPLEX (KIND(1d0)), INTENT (in) :: X
       COMPLEX (KIND(1d0)) :: F
     END FUNCTION F
  END INTERFACE

  X0 = (0.1D0, 0.0D0)
  X1 = (0.2D0, 0.0D0)
  X2 = (0.5D0, 0.0D0)

  F0 = F(X0)
  F1 = F(X1)
  F2 = F(X2)

  DO
     W0 = (F2-F0)/(X2-X0)
     W1 = (F2-F1)/(X2-X1)

     A = (W1-W0)/(X1-X0)
     B = W0+A*(X2-X0)
     C = F2

     P = B*B-4*A*C

     D1 = B + SQRT(P)
     D2 = B - SQRT(P)

     IF (ABS(D1) > ABS(D2)) THEN
        D = D1
     ELSE
        D = D2
     ENDIF

     X = X2 - 2*C / D

     T = F(X)
     !     termination condition
     IF((ABS(X - X2) < DELTA) .OR. (ABS(T) < EPSI)) EXIT

     X0 = X1
     X1 = X2
     X2 = X
     F0 = F1
     F1 = F2
     F2 = T
  ENDDO

  WRITE (*,*) "Η ρίζα είναι", X
  WRITE (*,*) "Η τιμή της συνάρτησης είναι", F(X)

END PROGRAM MULLERS


!     Define f(x)
FUNCTION F(X)
  IMPLICIT NONE
  COMPLEX(KIND(1D0)), INTENT (in) :: X
  COMPLEX(KIND(1D0)) :: F

  F = SIN(X) - X**2
END FUNCTION F
