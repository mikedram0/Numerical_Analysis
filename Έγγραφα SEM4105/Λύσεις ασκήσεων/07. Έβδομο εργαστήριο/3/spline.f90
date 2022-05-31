!     MAKE LOWER TRIANGULAR
SUBROUTINE TRIANG(A, B)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (inout) :: A(:,:), B(:)

  INTERFACE
     SUBROUTINE pivot(a, b, i)
       IMPLICIT NONE

       DOUBLE PRECISION, INTENT (inout) :: A(:,:), B(:)
       INTEGER, INTENT (in) :: i
     END SUBROUTINE pivot
  END INTERFACE

  INTEGER :: N, I, K
  DOUBLE PRECISION :: G

  N = SIZE(B)


  DO K = 1, N
     CALL pivot(a, b, k)                  
     DO I=K+1, N
        G = -A(I,K) / A(K,K)
        A(I,K:N) = A(I,K:N) + A(K,K:N) * G
        B(I) = B(I) + B(K) * G
     ENDDO
  ENDDO

END SUBROUTINE TRIANG





!     BACKSUBSTITUTION
SUBROUTINE BACKSU(A, B, X)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in)  :: A(:,:), B(:)
  DOUBLE PRECISION, INTENT (out) :: X(:)

  INTEGER :: K, N

  N = SIZE(b)

  DO K = N,1,-1
     X(K) = (B(K) - DOT_PRODUCT(A(K,K+1:N), X(K+1:N)) ) / A(K,K)
  ENDDO

END SUBROUTINE BACKSU


!meriki odigisi (partial pivoting)
SUBROUTINE pivot(a, b, i)
  IMPLICIT NONE

  DOUBLE PRECISION, INTENT (inout) :: A(:,:), B(:)
  INTEGER, INTENT (in) :: i
  INTEGER :: j, k

  INTERFACE
     ELEMENTAL SUBROUTINE swap(a,b)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (inout) :: a, b
     END SUBROUTINE swap
  END INTERFACE

  ! find position of max element in i column below the diagonal.
  j = i
  DO k = i+1, SIZE(a,1)
     IF ( ABS(a(k,i)) > ABS(a(j,i)) ) j = k
  ENDDO
  ! in i column a(j,i) is maximum.

  IF (j > i) THEN              ! exchange i, j rows of A and B
     CALL swap(a(i,:), a(j,:))
     CALL swap(b(i), b(j))     
  END IF

END SUBROUTINE pivot

! a<->b
ELEMENTAL SUBROUTINE swap(a,b)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (inout) :: a, b

  DOUBLE PRECISION :: temp

  temp = a
  a    = b
  b    = temp

END SUBROUTINE swap


SUBROUTINE spline(x,y,a,b,c,d)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in)  :: x(0:), y(0:)
  DOUBLE PRECISION, INTENT (out) :: a(0:), b(0:), c(0:), d(0:)

  DOUBLE PRECISION :: matA(4*(SIZE(x)-1), 4*(SIZE(x)-1))
  DOUBLE PRECISION :: matB(4*(SIZE(x)-1)), matX(4*(SIZE(x)-1))

  INTEGER :: n, i

  INTERFACE 
     SUBROUTINE TRIANG(A, B)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (inout) :: A(:,:), B(:)
     END SUBROUTINE TRIANG

     SUBROUTINE BACKSU(A, B, X)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in)  :: A(:,:), B(:)
       DOUBLE PRECISION, INTENT (out) :: X(:)
     END SUBROUTINE BACKSU
  END INTERFACE

  n = SIZE(x)-1

  matA = 0.0d0
  matB = 0.0d0

  DO i=0,n-1
     matA(i+1,4*i+4) = 1.0d0

     matB(i+1) = y(i)
  ENDDO

  DO i=0,n-1
     matA(i+1+n,4*i+1) = (x(i+1)-x(i))**3
     matA(i+1+n,4*i+2) = (x(i+1)-x(i))**2
     matA(i+1+n,4*i+3) = (x(i+1)-x(i))
     matA(i+1+n,4*i+4) = 1.0d0

     matB(i+1+n) = y(i+1)
  ENDDO


  DO i=1,n-1
     matA(i+2*n, 4*(i-1)+1) = 3.0d0 * (x(i)-x(i-1))**2
     matA(i+2*n, 4*(i-1)+2) = 2.0d0 * (x(i)-x(i-1))
     matA(i+2*n, 4*(i-1)+3) = 1.0d0

     matA(i+2*n, 4*i+3) = -1.0d0
  ENDDO

  DO i=1,n-1
     matA(i-1+3*n, 4*(i-1)+1) = 3.0d0 * (x(i)-x(i-1))
     matA(i-1+3*n, 4*(i-1)+2) = 1.0d0

     matA(i-1+3*n, 4*i+2) = -1.0d0
  ENDDO

  matA(4*n-1,2) = 1.0d0

  matA(4*n,4*n-3) = 3.0d0 * (x(n)-x(n-1))
  matA(4*n,4*n-2) = 1.0d0


  CALL triang(matA, matB)
  CALL backsu(matA, matB, matX)
  
  DO i=0,n-1
     a(i) = matX(4*i+1)
     b(i) = matX(4*i+2)
     c(i) = matX(4*i+3)
     d(i) = matX(4*i+4)
  ENDDO

END SUBROUTINE spline


SUBROUTINE create_data(fname)
  IMPLICIT NONE
  CHARACTER(*), INTENT (in) :: fname

  INTEGER, PARAMETER :: n = 15
  DOUBLE PRECISION, PARAMETER :: a = 2.0d0, b = 4.0d0
  DOUBLE PRECISION, PARAMETER :: step = (b-a) / (n-1)
  DOUBLE PRECISION :: x, y
  INTEGER :: i

  OPEN(22, file=fname)

  WRITE (22, *) n
  DO i= 1, n
     x = a + step * (i-1)
     y = SIN(x)

     WRITE(22,*) x, y
  ENDDO

  CLOSE(22)
END SUBROUTINE create_data

SUBROUTINE read_data(fname, x, y)
  IMPLICIT NONE
  CHARACTER(*), INTENT (in) :: fname
  DOUBLE PRECISION, ALLOCATABLE, INTENT (out) :: x(:), y(:)

  INTEGER :: n, i

  OPEN(22, file=fname)

  READ (22, *) n

  ALLOCATE(x(n), y(n))

  DO i= 1, n
     READ (22,*) x(i), y(i)
  ENDDO

  CLOSE(22)
END SUBROUTINE read_data


FUNCTION poly(x, a, b, c, d, z)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: x(0:), a(0:), b(0:), c(0:), d(0:), z
  DOUBLE PRECISION :: poly

  INTEGER :: i, n
  DOUBLE PRECISION :: delta

  n = SIZE(x) - 1

  DO i=0,n-1
     IF (z <= x(i+1)) EXIT
  ENDDO

  delta = z-x(i)
  poly = a(i) * delta**3 + b(i) * delta**2 + c(i) * delta + d(i) 
END FUNCTION poly



PROGRAM sp
  IMPLICIT NONE

  DOUBLE PRECISION, ALLOCATABLE :: x(:), y(:), a(:), b(:), c(:), d(:)
  INTEGER, PARAMETER :: m = 100
  DOUBLE PRECISION :: min_x, max_x, step, xout, yout
  INTEGER :: i, n

  INTERFACE 
     SUBROUTINE spline(x,y,a,b,c,d)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in)  :: x(0:), y(0:)
       DOUBLE PRECISION, INTENT (out) :: a(0:), b(0:), c(0:), d(0:)
     END SUBROUTINE spline

     SUBROUTINE create_data(fname)
       IMPLICIT NONE
       CHARACTER(*), INTENT (in) :: fname
     END SUBROUTINE create_data

     SUBROUTINE read_data(fname, x, y)
       IMPLICIT NONE
       CHARACTER(*), INTENT (in) :: fname
       DOUBLE PRECISION, ALLOCATABLE, INTENT (out) :: x(:), y(:)
     END SUBROUTINE read_data

     FUNCTION poly(x, a, b, c, d, z)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: x(0:), a(0:), b(0:), c(0:), d(0:), z
       DOUBLE PRECISION :: poly
     END FUNCTION poly

  END INTERFACE

  CALL create_data("points.dat")
  CALL read_data("points.dat", x,y)

  n = SIZE(x)-1
  ALLOCATE(a(0:n-1), b(0:n-1), c(0:n-1), d(0:n-1))

  CALL spline(x,y,a,b,c,d)

  min_x = MINVAL(x)
  max_x = MAXVAL(x)
  step = (max_x - min_x) / (m-1)

  DO i = 1,m
     xout = min_x + step * (i-1)
     yout = poly(x, a,b,c,d, xout)

     PRINT *, xout, yout
  ENDDO

END PROGRAM SP
