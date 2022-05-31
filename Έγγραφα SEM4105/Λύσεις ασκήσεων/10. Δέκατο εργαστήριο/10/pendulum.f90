SUBROUTINE rk4(x, y, h, f, xnew, ynew)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in)  ::  x, y(:), h
  DOUBLE PRECISION, INTENT (out) :: xnew, ynew(:)

  INTERFACE
     FUNCTION f(i,x,y)
       IMPLICIT NONE
       INTEGER, INTENT (in) :: i
       DOUBLE PRECISION, INTENT (in) :: x, y(:)
       DOUBLE PRECISION :: f
     END FUNCTION f
  END INTERFACE

  DOUBLE PRECISION, DIMENSION (SIZE(y)) :: k1, k2, k3, k4
  INTEGER :: i

  DO i =1,SIZE(y)
     k1(i) = h * f(i,x,y)
  ENDDO
  DO i =1,SIZE(y)
     k2(i) = h * f(i,x + h/2.0d0, y + k1/2.0d0)
  ENDDO
  DO i =1,SIZE(y)
     k3(i) = h * f(i,x + h/2.0d0, y + k2/2.0d0)
  ENDDO
  DO i =1,SIZE(y)
     k4(i) = h * f(i,x + h, y + k3)
  ENDDO
  
  ynew = y + (k1 + 2.0d0 * (k2 + k3) + k4) / 6.0d0

  xnew = x + h
END SUBROUTINE rk4



! The system is 
!   D theta / D t  = z
!   D z / D t      = -sin(theta)
!

! theta -> y(1), z-> y(2).  D theta / D t -> f(1), D z / D t -> f(2). t -> x
FUNCTION f(i, x, y)
  IMPLICIT NONE

  INTEGER, INTENT (in) :: i
  DOUBLE PRECISION, INTENT (in) ::  x, y(:)
  DOUBLE PRECISION :: f

  IF (i == 1) f = y(2)
  IF (i == 2) f = -SIN(y(1)) 
END FUNCTION f


! solution of theta'' = - theta:   theta = A cos(t) + B sin(t)
! theta(0) = 45 deg. , theta'(0) = 0 => 
SUBROUTINE solut(t,y)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: t
  DOUBLE PRECISION, INTENT (out) :: y(:)

  DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979323846d0
  DOUBLE PRECISION, PARAMETER :: A = 45.0d0 / 180.0d0 * pi

  y(1) =  A * COS(t)
  y(2) = -A * SIN(t)
END SUBROUTINE solut

PROGRAM pendulum
  IMPLICIT NONE

  INTERFACE
     SUBROUTINE rk4(x, y, h, f, xnew, ynew)
       DOUBLE PRECISION, INTENT (in)  ::  x, y(:), h
       DOUBLE PRECISION, INTENT (out) :: xnew, ynew(:)

       INTERFACE
          FUNCTION f(i,x,y)
            IMPLICIT NONE
            INTEGER, INTENT (in) :: i
            DOUBLE PRECISION, INTENT (in) :: x, y(:)
            DOUBLE PRECISION :: f
          END FUNCTION f
       END INTERFACE
     END SUBROUTINE rk4
  END INTERFACE


  INTERFACE
     FUNCTION f(i,x,y)
       IMPLICIT NONE
       INTEGER, INTENT (in) :: i
       DOUBLE PRECISION, INTENT (in) :: x, y(:)
       DOUBLE PRECISION :: f
     END FUNCTION f

     SUBROUTINE solut(t,y)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: t
       DOUBLE PRECISION, INTENT (out) :: y(:)
     END SUBROUTINE solut
  END INTERFACE

  DOUBLE PRECISION :: tinit, tfin, h
  DOUBLE PRECISION :: yinit(2), yold(2), ynew(2), yapprox(2), told, tnew
  DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979323846d0
  INTEGER :: k, niter
  
  tinit = 0.0d0
  tfin = 10.0d0

  h= 0.1d0
  yinit(1) = 45.0d0 / 180.0d0 * pi
  yinit(2) = 0.0d0

  told = tinit
  yold = yinit

  niter = NINT((tfin-tinit)/h)

  DO k = 1, niter
     CALL rk4(told,yold,h,f,tnew,ynew)
     told = tnew
     yold = ynew
     
     CALL solut(tnew, yapprox)
     
     PRINT "(5G16.8)", tnew, ynew, yapprox
  ENDDO

END PROGRAM pendulum
