FUNCTION force(x)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: x
  DOUBLE PRECISION :: force

  force = x * (1.0d0 - 0.01d0*x*x)
END FUNCTION force




!
!  m x'' = F(x) =>    x' = v     =>    y0 = x, y1 = v
!                     v' = F/m
!
SUBROUTINE f(t, y, dy)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in)  :: t, y(:)
  DOUBLE PRECISION, INTENT (out) :: dy(:)

  DOUBLE PRECISION, PARAMETER :: mass = 2.0d0
  INTERFACE
     FUNCTION force(x)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: x
       DOUBLE PRECISION :: force
     END FUNCTION force
  END INTERFACE

  dy(1) = y(2)
  dy(2) = force(y(1)) / mass
END SUBROUTINE f


SUBROUTINE ralston(t0, y0, t1, y1, f)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: t0, y0(:), t1
  DOUBLE PRECISION, INTENT (out) :: y1(:)
  INTERFACE
     SUBROUTINE f(t, y, dy)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in)  :: t, y(:)
       DOUBLE PRECISION, INTENT (out) :: dy(:)
     END SUBROUTINE f
  END INTERFACE

  DOUBLE PRECISION, DIMENSION(SIZE(y0)) :: k1, k2, p
  DOUBLE PRECISION :: h

  h = t1-t0

  CALL f(t0, y0, k1)
  k1 = k1 * h

  p = y0 + 2.0d0/3.0d0 * k1

  CALL f(t0+2.0d0/3.0d0*h, p, k2)
  k2 = k2 * h

  y1 = y0 + (k1 + 3d0 * k2) / 4d0

END SUBROUTINE ralston


PROGRAM sysrk
  IMPLICIT NONE
  INTERFACE
     SUBROUTINE f(t, y, dy)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in)  :: t, y(:)
       DOUBLE PRECISION, INTENT (out) :: dy(:)
     END SUBROUTINE f
     SUBROUTINE ralston(t0, y0, t1, y1, f)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: t0, y0(:), t1
       DOUBLE PRECISION, INTENT (out) :: y1(:)
       INTERFACE
          SUBROUTINE f(t, y, dy)
            IMPLICIT NONE
            
            DOUBLE PRECISION, INTENT (in)  :: t, y(:)
            DOUBLE PRECISION, INTENT (out) :: dy(:)
          END SUBROUTINE f
       END INTERFACE
     END SUBROUTINE ralston
  END INTERFACE


  INTEGER, PARAMETER :: nsteps = 100000
  DOUBLE PRECISION, PARAMETER :: h = 1d-3

  DOUBLE PRECISION :: t0, x0, v0

  DOUBLE PRECISION ::  y0(2), y1(2)
  INTEGER :: i

  t0 = 0.0d0
  x0 = 2.5d-2
  v0 = 0.0d0

  y0 = (/ x0, v0 /)

  OPEN(14, file="sysrk.txt", action="write", status = "replace")

  DO i = 1, nsteps

     WRITE(14, *) t0, y0

     CALL ralston(t0, y0, t0+h, y1, f)

     t0 = t0 + h

     y0 = y1
  END DO

  CLOSE(14)
END PROGRAM sysrk


