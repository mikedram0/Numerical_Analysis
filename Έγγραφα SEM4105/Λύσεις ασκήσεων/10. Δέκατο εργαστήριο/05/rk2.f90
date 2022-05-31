FUNCTION correct(x)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: x
  DOUBLE PRECISION :: correct

  correct = 1.0d0 - EXP(-x) + x*x -x
END FUNCTION correct


FUNCTION f(x,y)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: x, y
  DOUBLE PRECISION :: f

  f = x*x + x-y
END FUNCTION f


FUNCTION heun(x0, y0, x1, f)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: x0, y0, x1
  INTERFACE
     FUNCTION f(x,y)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: x, y
       DOUBLE PRECISION :: f
     END FUNCTION f
  END INTERFACE
  DOUBLE PRECISION :: heun

  DOUBLE PRECISION :: h, k1, k2
  h = x1-x0

  k1 = h * f(x0,y0)
  k2 = h * f(x1, y0 + k1)

  heun = y0 + (k1 + k2)/2.0
END FUNCTION heun


FUNCTION ralston(x0, y0, x1, f)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: x0, y0, x1
  INTERFACE
     FUNCTION f(x,y)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: x, y
       DOUBLE PRECISION :: f
     END FUNCTION f
  END INTERFACE
  DOUBLE PRECISION :: ralston

  DOUBLE PRECISION :: h, k1, k2

  h = x1-x0

  k1 = h * f(x0,y0)
  k2 = h * f(x0 +2.0d0/3.0d0 * h, y0 + 2.0d0/3.0d0*k1)

  ralston = y0 + (k1 + 3.0d0*k2)/4.0d0
END FUNCTION ralston



PROGRAM rk2
  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: a = 0.0d0
  DOUBLE PRECISION, PARAMETER :: b = 0.6d0
  DOUBLE PRECISION, PARAMETER :: yinit = 0.0d0
  DOUBLE PRECISION :: h, x, yh, yr

  INTERFACE
     FUNCTION f(x,y)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: x, y
       DOUBLE PRECISION :: f
     END FUNCTION f

     FUNCTION correct(x)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: x
       DOUBLE PRECISION :: correct
     END FUNCTION correct

     FUNCTION heun(x0, y0, x1, f)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: x0, y0, x1
       INTERFACE
          FUNCTION f(x,y)
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT (in) :: x, y
            DOUBLE PRECISION :: f
          END FUNCTION f
       END INTERFACE
       DOUBLE PRECISION :: heun
     END FUNCTION heun

     FUNCTION ralston(x0, y0, x1, f)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: x0, y0, x1
       INTERFACE
          FUNCTION f(x,y)
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT (in) :: x, y
            DOUBLE PRECISION :: f
          END FUNCTION f
       END INTERFACE
       DOUBLE PRECISION :: ralston
     END FUNCTION ralston
  END INTERFACE

  h = 1d-2
  x = a
  yh = yinit
  yr = yinit

  DO
     IF (x >= b) EXIT
     IF (x+h > b) h = b-x


     yh = heun(x, yh, x+h, f)
     yr = ralston(x, yr, x+h, f)

     x = x+h

     PRINT *, x, yh, yr, correct(x)
  END DO

END PROGRAM rk2
