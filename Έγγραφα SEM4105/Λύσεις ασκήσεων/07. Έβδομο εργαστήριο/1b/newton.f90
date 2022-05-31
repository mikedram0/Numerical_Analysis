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

  OPEN(22, file=fname, status="old", action="read")

  READ (22, *) n

  ALLOCATE(x(0:n-1), y(0:n-1))

  DO i= 0,n-1
     READ (22,*) x(i), y(i)
  ENDDO

  CLOSE(22)
END SUBROUTINE read_data


FUNCTION q(i,x,z)
  IMPLICIT NONE
  INTEGER, INTENT (in) :: i
  DOUBLE PRECISION, INTENT (in) :: x(0:)
  DOUBLE PRECISION, INTENT (in) :: z
  DOUBLE PRECISION :: q

  INTEGER :: j

  q = 1d0
  DO j=0,i-1
     q = q * (z-x(j))
  END DO

END FUNCTION q


SUBROUTINE coefficients(x, y, a)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in)  :: x(0:), y(0:)
  DOUBLE PRECISION, ALLOCATABLE, INTENT (out) :: a(:)

  INTEGER :: n, i, j
  DOUBLE PRECISION :: s

  INTERFACE
     FUNCTION q(i,x,z)
       IMPLICIT NONE
       INTEGER, INTENT (in) :: i
       DOUBLE PRECISION, INTENT (in) :: x(0:)
       DOUBLE PRECISION, INTENT (in) :: z
       DOUBLE PRECISION :: q
     END FUNCTION q
  END INTERFACE

  n = SIZE(x)-1

  ALLOCATE(a(0:n))

  DO  i = 0, n        
     s = y(i)
     DO j = 0, i-1
        s = s - a(j) * q(j, x, x(i)) 
     ENDDO
     
     a(i) = s / q(i, x, x(i))
  ENDDO
  
END SUBROUTINE coefficients


FUNCTION poly(a, x, z)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: a(0:), x(0:), z
  DOUBLE PRECISION :: poly

  INTEGER :: i, n

  INTERFACE
     FUNCTION q(i,x,z)
       IMPLICIT NONE
       INTEGER, INTENT (in) :: i
       DOUBLE PRECISION, INTENT (in) :: x(0:)
       DOUBLE PRECISION, INTENT (in) :: z
       DOUBLE PRECISION :: q
     END FUNCTION q
  END INTERFACE

  n = SIZE(x)-1

  poly = 0d0
  DO i = 0,n
     poly = poly + a(i) * q(i, x, z)
  ENDDO

END FUNCTION poly


PROGRAM newton
  IMPLICIT NONE

  DOUBLE PRECISION, ALLOCATABLE :: x(:), y(:), a(:)
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

     SUBROUTINE coefficients(x, y, a)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in)  :: x(0:), y(0:)
       DOUBLE PRECISION, ALLOCATABLE, INTENT (out) :: a(:)
     END SUBROUTINE coefficients

     FUNCTION poly(a, x, z)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: a(0:), x(0:), z
       DOUBLE PRECISION :: poly
     END FUNCTION poly
  END INTERFACE

  CALL create_data("points.dat")
  CALL read_data("points.dat", x,y)
  CALL coefficients(x,y,a)


  min_x = MINVAL(x)
  max_x = MAXVAL(x)
  step = (max_x - min_x) / (m-1)

  DO i = 1,m
     xout = min_x + step * (i-1)
     yout = poly(a, x, xout)

     PRINT *, xout, yout
  ENDDO

END PROGRAM newton
