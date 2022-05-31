PROGRAM gauss_hermite
  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: pi=3.14159265358979323846d0
  INTEGER, PARAMETER :: n=4
  INTEGER :: i
  DOUBLE PRECISION :: x(n), w, fval, sum
  INTERFACE
     FUNCTION f(x)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) ::x
       DOUBLE PRECISION :: f
     END FUNCTION f

     FUNCTION weight(n, x)
       IMPLICIT NONE
       INTEGER, INTENT (in) :: n
       DOUBLE PRECISION, INTENT (in):: x
       DOUBLE PRECISION :: weight
     END FUNCTION weight
  END INTERFACE

  x(1) = SQRT ( (3.0d0-SQRT(6.0d0)) / 2.0d0 )
  x(2) = -x(1)
  x(3) = SQRT ( (3.0d0+SQRT(6.0d0)) / 2.0d0 )
  x(4) = -x(3) 

  sum = 0.0d0
  DO i=1,n
     w = weight(n, x(i))
     fval = f(x(i))
     sum = sum + w * fval
  ENDDO

  PRINT *, "the integral is", sum
  PRINT *, "the correct value is", SQRT(pi) / 2.0d0
  PRINT *, "exactly the same (why?)"
END PROGRAM gauss_hermite

FUNCTION f(x)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) ::x
  DOUBLE PRECISION :: f
  f = x*x
END FUNCTION f

FUNCTION parag(n)
  IMPLICIT NONE
  INTEGER, INTENT (in) :: n
  INTEGER :: parag

  INTEGER :: temp, i

  temp = 1
  DO i=2,n
     temp = temp * i
  ENDDO
  parag = temp
END FUNCTION parag

FUNCTION weight(n, x)
  IMPLICIT NONE
  INTEGER, INTENT (in) :: n
  DOUBLE PRECISION, INTENT (in):: x
  DOUBLE PRECISION :: weight

  DOUBLE PRECISION, PARAMETER :: pi=3.14159265358979323846d0

  INTERFACE
     FUNCTION parag(n)
       IMPLICIT NONE
       INTEGER, INTENT (in) :: n
       INTEGER :: parag
     END FUNCTION parag
     FUNCTION hermite(n, x)
       IMPLICIT NONE
       INTEGER, INTENT (in) :: n
       DOUBLE PRECISION, INTENT (in) :: x
       DOUBLE PRECISION :: hermite
     END FUNCTION hermite
  END INTERFACE

  weight = 2.0d0**(n-1) * parag(n) * SQRT(pi) / (n * hermite(n-1,x))**2

END FUNCTION weight

FUNCTION hermite(n, x)
  IMPLICIT NONE
  INTEGER, INTENT (in) :: n
  DOUBLE PRECISION, INTENT (in) :: x
  DOUBLE PRECISION :: hermite

  SELECT CASE(n)
  CASE (0)
     hermite = 1.0d0
  CASE (1) 
     hermite = 2.0d0 * x
  CASE (2)
     hermite = 4.0d0 * x * x - 2.0d0
  CASE (3)
     hermite = x * (8.0d0 * x * x - 12.0d0)
  CASE (4) 
     hermite = x * x * (16.0d0 * x * x - 48.0d0) + 12.0d0
  CASE default
     STOP "not implemented"
  END SELECT
END FUNCTION hermite
