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

SUBROUTINE difeq(x,y,dy)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in)  :: x,y
  DOUBLE PRECISION, INTENT (out) :: dy(5)
  
  dy(1) = COS(x)-SIN(y) + x*x
  dy(2) = 2.0d0 * x - SIN(x) - COS(y) * dy(1)
  dy(3) = 2.0d0 - COS(x) - COS(y) * dy(2) + SIN(y) * dy(1)**2 
  dy(4) = SIN(x) + COS(y) * (dy(1)**3 - dy(3)) + 3.0d0 * dy(1) * dy(2) * SIN(y)
  dy(5) = COS(x) + COS(y) * (6.0d0 * dy(1)**2 * dy(2) - dy(4)) &
       + (3.0d0 * dy(2)**2 + dy(1) * (4.0d0 * dy(3) - dy(1)**3)) * SIN(y)  
END SUBROUTINE difeq


SUBROUTINE taylor(a,b,h,ya)
  IMPLICIT NONE

  DOUBLE PRECISION, INTENT (in) :: a, b, h, ya
  INTERFACE
     SUBROUTINE difeq(x,y,dy)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in)  :: x,y
       DOUBLE PRECISION, INTENT (out) :: dy(5)
     END SUBROUTINE difeq

     FUNCTION parag(n)
       IMPLICIT NONE
       INTEGER, INTENT (in) :: n
       INTEGER :: parag
     end FUNCTION parag
  END INTERFACE

  DOUBLE PRECISION :: xold, yold, xnew, ynew, dy(5)
  INTEGER :: i

  xold = a
  yold = ya
  DO
     CALL difeq(xold, yold, dy)

     xnew = xold + h

     ynew = yold
     DO i=1,5
        ynew = ynew + dy(i) * h**i / parag(i)
     END DO

     PRINT *, xnew, ynew

     xold = xnew
     yold = ynew

     IF (xnew > b) EXIT
  ENDDO
END SUBROUTINE taylor


PROGRAM taylor_prog
  IMPLICIT NONE

  INTERFACE
     SUBROUTINE taylor(a,b,h,ya)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: a, b, h, ya
     END SUBROUTINE taylor
  END INTERFACE

  CALL taylor(-1.0d0, 1.0d0, 0.01d0, 3.0d0)
END PROGRAM taylor_prog
