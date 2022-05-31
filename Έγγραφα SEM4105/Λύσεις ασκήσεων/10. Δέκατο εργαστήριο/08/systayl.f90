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

SUBROUTINE deriv(x,y,z, dy, dz)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: x, y, z
  DOUBLE PRECISION, INTENT (out) :: dy(4), dz(4)

  dy(1) = y + z*z - x**3
  dz(1) = y**3 + z + COS(x)

  dy(2) =-3.0d0 * x*x + dy(1) + 2.0d0 * z * dz(1)
  dz(2) = 3.0d0 * y*y * dy(1) + dz(1) - SIN(x)

  dy(3) = -6.0d0 * x + dy(2) + 2.0d0 * (dz(1)**2 + z * dz(2))
  dz(3) = 3.0d0 * y * ( 2.0d0 * dy(1)**2 + y * dy(2)) + dz(2) - COS(x)

  dy(4) =-6.0d0 + dy(3) + 2.0d0 * (3.0d0*dz(1) * dz(2)+ z*dz(3)) 
  dz(4) = dz(3) + SIN(x) + 6.0d0 * dy(1)**3 + &
       3.0d0 * y * ( 6.0d0 * dy(1) * dy(2) + y * dy(3))
END SUBROUTINE deriv

SUBROUTINE system_taylor(a, b, h, ya, za)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: a, b, h, ya, za

  DOUBLE PRECISION :: xold, yold, zold, xnew, ynew, znew, dy(4), dz(4)
  INTERFACE
     FUNCTION parag(n)
       IMPLICIT NONE
       INTEGER, INTENT (in) :: n
       INTEGER :: parag
     END FUNCTION parag

     SUBROUTINE deriv(x,y,z, dy, dz)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: x, y, z
       DOUBLE PRECISION, INTENT (out) :: dy(4), dz(4)
     END SUBROUTINE deriv
  END INTERFACE

  INTEGER :: i, nsteps, k

  xold = a
  yold = ya
  zold = za

  nsteps = NINT((b-a)/h)
  
  DO k = 1,nsteps
     CALL deriv(xold, yold, zold, dy, dz)

     xnew = xold + h

     ynew = yold
     znew = zold
     DO i=1,4
        ynew = ynew + dy(i) * h**i / parag(i)
        znew = znew + dz(i) * h**i / parag(i)
     END DO

     PRINT *, xnew, ynew, znew

     xold = xnew
     yold = ynew
     zold = znew
  ENDDO
  
END SUBROUTINE system_taylor

PROGRAM systaylor
  IMPLICIT NONE

  INTERFACE
     SUBROUTINE system_taylor(a, b, h, ya, za)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: a, b, h, ya, za
     END SUBROUTINE system_taylor
  END INTERFACE

  CALL system_taylor(0.0d0, 1.0d0, 0.1d0, 0.3d0, 0.1d0)
END PROGRAM systaylor
