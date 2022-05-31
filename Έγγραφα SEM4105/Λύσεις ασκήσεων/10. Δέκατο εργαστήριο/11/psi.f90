!  y' = z
!  z' = (x*x-5) * y   
!
!
SUBROUTINE f(x, y, dy)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in)  :: x, y(:)
  DOUBLE PRECISION, INTENT (out) :: dy(:)

  dy(1) = y(2)
  dy(2) = (x*x-5d0) * y(1)
END SUBROUTINE f


SUBROUTINE ralston(x0, y0, x1, y1, f)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: x0, y0(:), x1
  DOUBLE PRECISION, INTENT (out) :: y1(:)
  INTERFACE
     SUBROUTINE f(x, y, dy)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in)  :: x, y(:)
       DOUBLE PRECISION, INTENT (out) :: dy(:)
     END SUBROUTINE f
  END INTERFACE

  DOUBLE PRECISION, DIMENSION(SIZE(y0)) :: k1, k2, p
  DOUBLE PRECISION :: h

  h = x1-x0

  CALL f(x0, y0, k1)
  k1 = k1 * h

  p = y0 + 2.0d0/3.0d0 * k1

  CALL f(x0+2.0d0/3.0d0*h, p, k2)
  k2 = k2 * h

  y1 = y0 + (k1 + 3d0 * k2) / 4d0

END SUBROUTINE ralston



PROGRAM psi
  IMPLICIT NONE

  DOUBLE PRECISION, PARAMETER :: pi=3.14159265358979323846d0
  
  INTEGER, PARAMETER :: nsteps = 50
  DOUBLE PRECISION, PARAMETER :: a=0d0, b=2.0d0, h = (b-a)/nsteps

  DOUBLE PRECISION ::  x, y(2,-nsteps:nsteps), y0(2), y1(2)
  INTEGER :: i

  INTERFACE
     SUBROUTINE f(x, y, dy)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in)  :: x, y(:)
       DOUBLE PRECISION, INTENT (out) :: dy(:)
     END SUBROUTINE f
     SUBROUTINE ralston(x0, y0, x1, y1, f)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: x0, y0(:), x1
       DOUBLE PRECISION, INTENT (out) :: y1(:)
       INTERFACE
          SUBROUTINE f(x, y, dy)
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT (in)  :: x, y(:)
            DOUBLE PRECISION, INTENT (out) :: dy(:)
          END SUBROUTINE f
       END INTERFACE
     END SUBROUTINE ralston
  END INTERFACE


  x = 0.0d0
  
  y(:, 0) = (/ -1.0d0/SQRT(2d0*SQRT(pi)), 0.0d0 /)

  DO i = 0,nsteps-1
     y0(:) = y(:,i)
     CALL ralston(x, y0, x+h, y1, f)
     x = x + h

     y(:,i+1) = y1
  END DO

  x = 0.0d0
  
  DO i = 0,-(nsteps-1), -1
     y0(:) = y(:,i)
     CALL ralston(x, y0, x-h, y1, f)
     x = x - h

     y(:,i-1) = y1
  END DO


  OPEN(14, file="psi.txt", action="write", status = "replace")

  DO i = -nsteps, nsteps
     x = i * h
     
     WRITE(14, *) x, y(:,i)
  END DO

  CLOSE(14)
END PROGRAM psi


