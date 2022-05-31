PROGRAM simpson_13_38
  IMPLICIT NONE

  DOUBLE PRECISION, ALLOCATABLE :: x(:), y(:)

  INTERFACE
     SUBROUTINE read_data(fname, x, y)
       IMPLICIT NONE
       CHARACTER(*), INTENT (in) :: fname
       DOUBLE PRECISION, ALLOCATABLE, INTENT (out) :: x(:), y(:)
     END SUBROUTINE read_data

     FUNCTION simpson_integral(x, y)
       IMPLICIT NONE

       DOUBLE PRECISION, INTENT (in) :: x(0:), y(0:)
       DOUBLE PRECISION :: simpson_integral
     END FUNCTION simpson_integral

     FUNCTION fint(x)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: x
       DOUBLE PRECISION :: fint
     END FUNCTION fint
  END INTERFACE

  CALL read_data("points.dat", X, Y)

  PRINT *, simpson_integral(x, y), fint(3.0d0) - fint(0.0d0)

END PROGRAM simpson_13_38

FUNCTION fint(x)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: x
  DOUBLE PRECISION :: fint

  fint = (SIN(x) - COS(x)) / 2.0d0 * EXP(x)
END FUNCTION fint

!     Read data file
SUBROUTINE read_data(fname, x, y)
  IMPLICIT NONE
  CHARACTER(*), INTENT (in) :: fname
  DOUBLE PRECISION, ALLOCATABLE, INTENT (out) :: x(:), y(:)

  INTEGER, PARAMETER ::  UNITNO = 54
  INTEGER :: n, i

  OPEN (UNIT = UNITNO, FILE = FNAME, status = 'old', action = 'read')
  READ (UNITNO, *) N

  n = n-1

  ALLOCATE(x(0:n), y(0:n))

  DO I=0,n
     READ (UNITNO, *) X(I), Y(I)
  ENDDO

  CLOSE (UNITNO)
END SUBROUTINE read_data

FUNCTION simpson_integral(x, y)
  IMPLICIT NONE

  DOUBLE PRECISION, INTENT (in) :: x(0:), y(0:)
  DOUBLE PRECISION :: simpson_integral

  INTEGER :: n

  INTERFACE
     FUNCTION simpson(a, b, y)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) ::  a, b, y(0:)
       DOUBLE PRECISION :: simpson
     END FUNCTION simpson

     FUNCTION simpson38(a, b, y)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) ::  a, b, y(0:3)
       DOUBLE PRECISION :: simpson38
     END FUNCTION simpson38
  END INTERFACE

  n = SIZE(x) - 1

  IF (MOD(n,2) == 1) THEN 
     simpson_integral = simpson38(x(0), x(3), y(0:3)) &
          + simpson(x(3), x(n), y(3:n)) 
  ELSE
     simpson_integral = simpson(x(0), x(n), y(0:n)) 
  ENDIF

END FUNCTION simpson_integral


FUNCTION simpson(a, b, y)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) ::  a, b, y(0:)
  DOUBLE PRECISION :: simpson

  DOUBLE PRECISION :: h
  INTEGER :: n

  n = SIZE(y) - 1

  h = (b-a)/n

  IF (MOD(n, 2) /= 0) STOP "odd size in simpson"

  simpson = y(0) + 4.0d0 * SUM(y(1:n-1:2)) + 2.0d0 * SUM(y(2:n-2:2)) + y(n)
  simpson = h * simpson / 3.0d0

END FUNCTION simpson

FUNCTION simpson38(a, b, y)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) ::  a, b, y(0:3)
  DOUBLE PRECISION :: simpson38

  DOUBLE PRECISION :: h

  h = (b-a) / 3d0

  simpson38 = 3.0d0 * h * (y(0) + 3d0 * (y(1) + y(2)) + y(3)) / 8d0

END FUNCTION simpson38



