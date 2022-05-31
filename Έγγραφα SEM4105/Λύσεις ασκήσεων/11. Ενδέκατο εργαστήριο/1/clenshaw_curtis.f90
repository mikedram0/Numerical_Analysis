FUNCTION c(i,n)
  IMPLICIT NONE
  INTEGER, INTENT (in) :: i, n
  DOUBLE PRECISION :: c

  c = 2.0d0
  IF (MOD(i,n) == 0) c = 1.0d0
END FUNCTION c

FUNCTION b(i,n)
  IMPLICIT NONE
  INTEGER, INTENT (in) :: i, n
  DOUBLE PRECISION :: b

  b = 2.0d0
  IF (MOD(i,n/2) == 0) b = 1.0d0
END FUNCTION b


SUBROUTINE coefs(w)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (out) :: w(0:)

  INTEGER :: n, i,j
  DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979323846d0

  INTERFACE
     FUNCTION c(i,n)
       IMPLICIT NONE
       INTEGER, INTENT (in) :: i, n
       DOUBLE PRECISION :: c
     END FUNCTION c

     FUNCTION b(i,n)
       IMPLICIT NONE
       INTEGER, INTENT (in) :: i, n
       DOUBLE PRECISION :: b
     END FUNCTION b
  END INTERFACE

  n = SIZE(w)-1

  DO i=0,n
     w(i) = 0.0d0
     DO j =0, n/2
        w(i) = w(i) + b(j,n) / (1.d0-4d0*j*j) * COS(2d0*i*j*pi/n)
     END DO
     w(i) = w(i) * c(i,n)/n
  END DO

END SUBROUTINE coefs


!  \int_{-2}^{2} \frac{1}{1+x^2} \D x = \int_{-1}^{1} \frac{2}{1+4x^2} \D x
FUNCTION f(x)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(in) :: x
  DOUBLE PRECISION :: f

  f = 2.0d0 / (1.0d0+ 4.0d0*x*x)
END FUNCTION f



FUNCTION calc(n)
  IMPLICIT NONE
  INTEGER, INTENT (in) :: n
  DOUBLE PRECISION :: calc

  INTEGER :: i
  DOUBLE PRECISION :: w(0:n), x
  DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979323846d0

  INTERFACE
     FUNCTION f(x)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT(in) :: x
       DOUBLE PRECISION :: f
     END FUNCTION f


     SUBROUTINE coefs(w)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (out) :: w(0:)
     END SUBROUTINE coefs
  END INTERFACE

  CALL coefs(w)

  calc = 0d0
  DO i=0,n
     x= COS(i*pi/n)
     calc = calc + w(i) * f(x)
  END DO

END FUNCTION calc


PROGRAM clenshaw_curtis
  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: correct = 2d0 * ATAN(2d0)

  INTEGER :: n
  DOUBLE PRECISION :: s, diff
  INTERFACE
     FUNCTION calc(n)
       IMPLICIT NONE
       INTEGER, INTENT (in) :: n
       DOUBLE PRECISION :: calc
     END FUNCTION calc
  END INTERFACE

  n = 3
  DO
     s = calc(n)
     diff = s - correct
     PRINT *, s, correct, diff
     IF (ABS(diff) < 1d-12) EXIT
     n=n+1
  END DO

  PRINT *, n
END PROGRAM clenshaw_curtis
