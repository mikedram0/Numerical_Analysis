FUNCTION f(x,y)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: x,y
  DOUBLE PRECISION :: f

  f = cos(x) - sin(y) + x*x
END FUNCTION f

FUNCTION forward_euler(x0, y0, x1, f)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: x0, y0, x1
  INTERFACE
     FUNCTION f(x,y)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: x,y
       DOUBLE PRECISION :: f
     END FUNCTION f
  END INTERFACE
  DOUBLE PRECISION :: forward_euler

  forward_euler = y0 + (x1-x0) * f(x0,y0)
END FUNCTION forward_euler


! y1 = y0 + (x1-x0) f(x1, y1)
FUNCTION g(y1, x0, y0, x1, f)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: x0, y0, x1, y1
  INTERFACE
     FUNCTION f(x,y)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: x,y
       DOUBLE PRECISION :: f
     END FUNCTION f
  END INTERFACE
  DOUBLE PRECISION :: g

  g = y1 - y0 - (x1-x0) * f(x1,y1)

END FUNCTION g


FUNCTION backward_euler(x0, y0, x1, f)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: x0, y0, x1
  INTERFACE
     FUNCTION f(x,y)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: x,y
       DOUBLE PRECISION :: f
     END FUNCTION f
  END INTERFACE
  DOUBLE PRECISION :: backward_euler

  INTERFACE
     FUNCTION g(y1, x0, y0, x1, f)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: x0, y0, x1, y1
       INTERFACE
          FUNCTION f(x,y)
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT (in) :: x,y
            DOUBLE PRECISION :: f
          END FUNCTION f
       END INTERFACE
       DOUBLE PRECISION :: g
     END FUNCTION g

     FUNCTION forward_euler(x0, y0, x1, f)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: x0, y0, x1
       INTERFACE
          FUNCTION f(x,y)
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT (in) :: x,y
            DOUBLE PRECISION :: f
          END FUNCTION f
       END INTERFACE
       DOUBLE PRECISION :: forward_euler
     END FUNCTION forward_euler
  END INTERFACE

  DOUBLE PRECISION :: y1, y2, y, g1, g2

  ! y will be the approximation to the unknown (y1)

  ! secant for g(y, x0,y0,x1,f) = 0

  ! need two points: one is y1 from forward_euler

  y1 = forward_euler(x0, y0, x1, f)

  !the other is y2, close to y1
  y2 = y1 * 1.1d0

  DO 
     g1 = g(y1, x0, y0, x1, f)
     g2 = g(y2, x0, y0, x1, f)

     y = y2 - g2 * (y2-y1) / (g2 - g1)

     IF (ABS(g(y, x0, y0, x1, f)) < 1d-8) EXIT

     y1 = y2
     y2 = y
  END DO

  backward_euler = y
END FUNCTION backward_euler

PROGRAM fb_euler
  IMPLICIT NONE

  INTERFACE
     FUNCTION f(x,y)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: x,y
       DOUBLE PRECISION :: f
     END FUNCTION f

     FUNCTION backward_euler(x0, y0, x1, f)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: x0, y0, x1
       INTERFACE
          FUNCTION f(x,y)
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT (in) :: x,y
            DOUBLE PRECISION :: f
          END FUNCTION f
       END INTERFACE
       DOUBLE PRECISION :: backward_euler
     END FUNCTION backward_euler

     FUNCTION forward_euler(x0, y0, x1, f)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: x0, y0, x1
       INTERFACE
          FUNCTION f(x,y)
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT (in) :: x,y
            DOUBLE PRECISION :: f
          END FUNCTION f
       END INTERFACE
       DOUBLE PRECISION :: forward_euler
     END FUNCTION forward_euler
  END INTERFACE

  DOUBLE PRECISION, PARAMETER :: a = -1.0d0
  DOUBLE PRECISION, PARAMETER :: yinit = 3.0d0
  DOUBLE PRECISION, PARAMETER :: b = 1.0d0
  DOUBLE PRECISION :: h 


  DOUBLE PRECISION :: x0, y0f, y0b
  INTEGER :: n, i

  x0 = a
  y0f = yinit
  y0b = yinit
  h = 1d-2

  n = NINT((b-a)/h)

  PRINT *, "X   Forward   Backward"

  DO i = 1, n
     ! change h if x0 is beyond b:
     IF (x0+h>b) h = b-x0

     y0f = forward_euler(x0, y0f, x0+h, f)
     y0b = backward_euler(x0, y0b, x0+h, f)
     x0 = x0 + h

     PRINT *, x0, y0f, y0b
  END DO

END PROGRAM fb_euler
