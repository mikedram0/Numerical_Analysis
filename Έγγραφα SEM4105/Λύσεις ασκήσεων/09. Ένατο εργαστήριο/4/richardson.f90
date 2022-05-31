FUNCTION g(x)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: x
  DOUBLE PRECISION :: g
  g = SIN(x)
END FUNCTION g

FUNCTION trapez(f, a, b, n)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: a, b
  INTEGER, INTENT (in) :: n
  DOUBLE PRECISION :: trapez

  INTERFACE
     FUNCTION f(x)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) ::  x
       DOUBLE PRECISION :: f
     END FUNCTION f
  END INTERFACE

  INTEGER :: i
  DOUBLE PRECISION :: x, sum, h

  h = (b-a) / n

  sum = 0.5d0 * (f(a) + f(b)) 

  DO i= 1,n-1
     x = a + i * h

     sum = sum + f(x)
  ENDDO

  trapez = sum * h
END FUNCTION trapez


FUNCTION richardson(f, a, b, n)
  IMPLICIT NONE

  DOUBLE PRECISION, INTENT (in) :: a, b
  INTEGER, INTENT (in) :: n
  DOUBLE PRECISION :: richardson

  INTERFACE 
     FUNCTION f(x)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: x
       DOUBLE PRECISION :: f
     END FUNCTION f
  END INTERFACE

  INTERFACE
     FUNCTION trapez(f, a, b, n)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: a, b
       INTEGER, INTENT (in) :: n
       DOUBLE PRECISION :: trapez
       INTERFACE
          FUNCTION f(x)
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT (in) ::  x
            DOUBLE PRECISION :: f
          END FUNCTION f
       END INTERFACE
     END FUNCTION trapez

     SUBROUTINE triang(a, b)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (inout) :: a(:,:), b(:)
     END SUBROUTINE triang

     SUBROUTINE backsu(a, b, x)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in)  :: a(:,:), b(:)
       DOUBLE PRECISION, INTENT (out) :: x(:)
     END SUBROUTINE backsu
  END INTERFACE

  DOUBLE PRECISION :: lhs(3,3), rhs(3), x(3), h

  h = (b-a)/n
  lhs(1,1) = 1.0d0
  lhs(2,1) = 1.0d0
  lhs(3,1) = 1.0d0

  lhs(1,2) = h**2 
  lhs(2,2) = (h/2)**2
  lhs(3,2) = (h/4)**2

  lhs(1,3) = h**4
  lhs(2,3) = (h/2)**4
  lhs(3,3) = (h/4)**4


  rhs(1) = trapez(f, a, b,   n)
  rhs(2) = trapez(f, a, b, 2 * n)
  rhs(3) = trapez(f, a, b, 4 * n)


  CALL triang(lhs, rhs)
  CALL backsu(lhs, rhs, x)

  richardson = x(1)

END FUNCTION richardson



!     make lower triangular
SUBROUTINE triang(a, b)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (inout) :: a(:,:), b(:)

  INTEGER :: n, i, k
  DOUBLE PRECISION :: g

  n = SIZE(b)

  DO i=1, n
     DO k=i+1, n
        g = -a(k,i) / a(i,i)
        a(k,1:n) = a(k,1:n) + a(i,1:n) * g
        b(k) = b(k) + b(i) * g
     ENDDO
  ENDDO

END SUBROUTINE triang


!     backsubstitution
SUBROUTINE backsu(a, b, x)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in)  :: a(:,:), b(:)
  DOUBLE PRECISION, INTENT (out) :: x(:)

  INTEGER :: k, n

  n = SIZE(b)

  DO k = n,1,-1
     x(k) = (b(k) - DOT_PRODUCT(a(k,k+1:n), x(k+1:n)) ) / a(k,k)
  ENDDO

END SUBROUTINE backsu


PROGRAM richardson_integration
  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979323846d0
  DOUBLE PRECISION, PARAMETER :: a = 0.0d0, b = pi
  INTEGER, PARAMETER :: n  = 16
  DOUBLE PRECISION, PARAMETER :: correct = 2.0d0

  INTERFACE 
     FUNCTION g(x)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: x
       DOUBLE PRECISION :: g
     END FUNCTION g

     FUNCTION trapez(f, a, b, n)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: a, b
       INTEGER, INTENT (in) :: n
       DOUBLE PRECISION :: trapez

       INTERFACE
          FUNCTION f(x)
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT (in) ::  x
            DOUBLE PRECISION :: f
          END FUNCTION f
       END INTERFACE

     END FUNCTION trapez

     FUNCTION richardson(f, a, b, n)
       IMPLICIT NONE

       DOUBLE PRECISION, INTENT (in) :: a, b
       INTEGER, INTENT (in) :: n
       DOUBLE PRECISION :: richardson

       INTERFACE 
          FUNCTION f(x)
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT (in) :: x
            DOUBLE PRECISION :: f
          END FUNCTION f
       END INTERFACE
     END FUNCTION richardson
  END INTERFACE

  PRINT *, "     spaces      error with trapez"
  PRINT *,   n, correct - trapez(g, a, b,   n)
  PRINT *, 2*n, correct - trapez(g, a, b, 2*n)
  PRINT *, 4*n, correct - trapez(g, a, b, 4*n)

  PRINT *
  PRINT *, "error with richardson: ", correct - richardson(g, a, b, n)

END PROGRAM richardson_integration
