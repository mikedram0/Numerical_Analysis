FUNCTION f(x,y)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: x,y
  DOUBLE PRECISION :: f

  f = COS(x)-x*SIN(x)
END FUNCTION f


FUNCTION correct(x)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: x
  DOUBLE PRECISION :: f

  f = 2.0d0 + x * COS(x)
END FUNCTION correct

SUBROUTINE euler(a, b, h, ya, f)
  IMPLICIT NONE

  DOUBLE PRECISION, INTENT (in) :: a, b, h, ya
  INTERFACE
     FUNCTION f(x,y)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: x,y
       DOUBLE PRECISION :: f
     END FUNCTION f
  END INTERFACE

  DOUBLE PRECISION :: xold, yold, xnew, ynew

  xold = a
  yold = ya
  DO
     xnew = xold + h
     ynew = yold + f(xold, yold) * h

     PRINT *, xnew, ynew

     xold = xnew
     yold = ynew

     IF (xnew > b) EXIT
  ENDDO
END SUBROUTINE euler

PROGRAM euler_integration
  IMPLICIT NONE

  INTERFACE
     SUBROUTINE euler(a, b, h, ya, f)
       IMPLICIT NONE

       DOUBLE PRECISION, INTENT (in) :: a, b, h, ya
       INTERFACE
          FUNCTION f(x,y)
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT (in) :: x,y
            DOUBLE PRECISION :: f
          END FUNCTION f
       END INTERFACE
     END SUBROUTINE euler

     FUNCTION f(x,y)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: x,y
       DOUBLE PRECISION :: f
     END FUNCTION f


     FUNCTION correct(x)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: x
       DOUBLE PRECISION :: f
     END FUNCTION correct
  END INTERFACE

  CALL euler(0.0d0, 3.0d0, 0.01d0, 2.0d0, f)

END PROGRAM euler_integration
