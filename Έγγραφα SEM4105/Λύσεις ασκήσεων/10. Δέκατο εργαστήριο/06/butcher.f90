FUNCTION runge_kutta(A, b, c, x0, y0, x1, f)
  IMPLICIT NONE

  DOUBLE PRECISION, INTENT (in) :: A(:,:), b(:), c(:)
  DOUBLE PRECISION, INTENT (in) :: x0, y0, x1
  INTERFACE
     FUNCTION f(x,y)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: x,y
       DOUBLE PRECISION :: f
     END FUNCTION f
  END INTERFACE
  DOUBLE PRECISION :: runge_kutta


  DOUBLE PRECISION :: h, sum
  DOUBLE PRECISION k(SIZE(b))
  INTEGER :: i

  h = x1-x0

  DO i = 1, SIZE(k)
     sum = y0 + DOT_PRODUCT(A(i,1:i), k(1:i))
     k(i) = h * f(x0+c(i) * h, sum)
  END DO

  runge_kutta = y0 + DOT_PRODUCT(b,k)
END FUNCTION runge_kutta


FUNCTION f(x,y)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: x,y
  DOUBLE PRECISION :: f

  f = (y+x)/(y-x)
END FUNCTION f




FUNCTION correct(x)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: x
  DOUBLE PRECISION :: correct

  correct = x + SQRT(1d0+2.0d0*x*x)

END FUNCTION correct



PROGRAM butcher
  IMPLICIT NONE

  DOUBLE PRECISION, PARAMETER :: a = 0.0d0
  DOUBLE PRECISION, PARAMETER :: yinit = 1.0d0
  DOUBLE PRECISION, PARAMETER :: b = 2.0d0
  DOUBLE PRECISION :: h, x,y

  INTEGER, PARAMETER :: s = 4
  DOUBLE PRECISION, PARAMETER :: rkA(s,s) = RESHAPE( (/   & ! column-wise !!
       0.0d0, 0.5d0, 0.0d0, 0.0d0, &
       0.0d0, 0.0d0, 0.5d0, 0.0d0, &
       0.0d0, 0.0d0, 0.0d0, 1.0d0, &
       0.0d0, 0.0d0, 0.0d0, 0.0d0 /), SHAPE(rkA))

  DOUBLE PRECISION, PARAMETER :: rkb(s) = &
       (/ 1.0d0/6.0d0, 1.0d0/3.0d0, 1.0d0/3.0d0, 1.0d0/6.0d0 /)

  DOUBLE PRECISION, PARAMETER :: rkc(s) = (/ 0.0d0, 0.5d0, 0.5d0, 1.0d0 /)


  INTERFACE
     FUNCTION runge_kutta(A, b, c, x0, y0, x1, f)
       IMPLICIT NONE

       DOUBLE PRECISION, INTENT (in) :: A(:,:), b(:), c(:)
       DOUBLE PRECISION, INTENT (in) :: x0, y0, x1
       INTERFACE
          FUNCTION f(x,y)
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT (in) :: x,y
            DOUBLE PRECISION :: f
          END FUNCTION f
       END INTERFACE
       DOUBLE PRECISION :: runge_kutta
     END FUNCTION runge_kutta

     FUNCTION f(x,y)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: x,y
       DOUBLE PRECISION :: f
     END FUNCTION f

     FUNCTION correct(x)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: x
       DOUBLE PRECISION :: correct
     END FUNCTION correct
  END INTERFACE

  h = 0.1d0

  x = a
  y = yinit

  DO
     IF (x >= b) EXIT
     IF (x+h>b) h = b-x

     y = runge_kutta(rkA, rkb, rkc, x, y, x+h, f)

     x = x + h

     PRINT *, x, y, correct(x)
  END DO

END PROGRAM butcher
