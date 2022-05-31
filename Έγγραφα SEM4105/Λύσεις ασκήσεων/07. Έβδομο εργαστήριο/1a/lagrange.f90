SUBROUTINE create_data(fname)
  IMPLICIT NONE
  CHARACTER(*), INTENT (in) :: fname

  INTEGER, PARAMETER :: n = 15
  DOUBLE PRECISION, PARAMETER :: a = 2.0d0, b = 4.0d0
  DOUBLE PRECISION, PARAMETER :: step = (b-a) / (n-1)
  DOUBLE PRECISION :: x, y
  INTEGER :: i

  OPEN(22, file=fname)

  WRITE (22, *) n
  DO i= 1, n
     x = a + step * (i-1)
     y = SIN(x)

     WRITE(22,*) x, y
  ENDDO

  CLOSE(22)
END SUBROUTINE create_data

SUBROUTINE read_data(fname, x, y)
  IMPLICIT NONE
  CHARACTER(*), INTENT (in) :: fname
  DOUBLE PRECISION, ALLOCATABLE, INTENT (out) :: x(:), y(:)

  INTEGER :: n, i

  OPEN(22, file=fname)

  READ (22, *) n

  ALLOCATE(x(n), y(n))

  DO i= 1, n
     READ (22,*) x(i), y(i)
  ENDDO

  CLOSE(22)
END SUBROUTINE read_data


FUNCTION ell(i, z, x)
  IMPLICIT NONE
  INTEGER, INTENT (in) :: i
  DOUBLE PRECISION, INTENT (in) :: z, x(:)
  DOUBLE PRECISION :: ell

  DOUBLE PRECISION :: p
  INTEGER :: j

  p = 1.0d0
  DO j = 1, SIZE(x)
     IF (i==j) CYCLE
     p = p * (z-x(j)) / (x(i)-x(j))
  ENDDO

  ell = p 

END FUNCTION ell


FUNCTION poly(x, y, z)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: x(:), y(:), z
  DOUBLE PRECISION :: poly

  INTEGER :: i
  DOUBLE PRECISION :: p

  INTERFACE
     FUNCTION ell(i, z, x)
       IMPLICIT NONE
       INTEGER, INTENT (in) :: i
       DOUBLE PRECISION, INTENT (in) :: z, x(:)
       DOUBLE PRECISION :: ell
     END FUNCTION ell
  END INTERFACE

  p = 0.0d0
  DO i= 1, SIZE(x)
     p = p + ell(i,z,x) * y(i)
  END DO

  poly = p
END FUNCTION poly



PROGRAM lagrange
  IMPLICIT NONE

  DOUBLE PRECISION, ALLOCATABLE :: x(:), y(:)
  INTEGER, PARAMETER :: m = 100
  DOUBLE PRECISION :: min_x, max_x, step, xout, yout
  INTEGER :: i

  INTERFACE
     SUBROUTINE create_data(fname)
       IMPLICIT NONE
       CHARACTER(*), INTENT (in) :: fname
     END SUBROUTINE create_data

     SUBROUTINE read_data(fname, x, y)
       IMPLICIT NONE
       CHARACTER(*), INTENT (in) :: fname
       DOUBLE PRECISION, ALLOCATABLE, INTENT (out) :: x(:), y(:)
     END SUBROUTINE read_data

     FUNCTION poly(x, y, z)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: x(:), y(:), z
       DOUBLE PRECISION :: poly
     END FUNCTION poly
  END INTERFACE

  CALL create_data("points.dat")
  CALL read_data("points.dat", x,y)

  min_x = MINVAL(x)
  max_x = MAXVAL(x)
  step = (max_x - min_x) / (m-1)

  xout = min_x
  DO i = 1,m
     yout = poly(x, y, xout)

     PRINT *, xout, yout, SIN(xout)
     xout = xout + step
  ENDDO

END PROGRAM lagrange

