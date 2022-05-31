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
SUBROUTINE pivot(A, B, k)
  IMPLICIT NONE

  DOUBLE PRECISION, INTENT (inout) :: A(:,:), B(:)
  INTEGER, INTENT (in) :: k
  INTEGER :: i, j, n
  DOUBLE PRECISION :: maxv

  INTERFACE
     ELEMENTAL SUBROUTINE swap(a,b)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (inout) :: a, b
     END SUBROUTINE swap
  END INTERFACE

  n = SIZE(a,1)

  ! 'normalize' equations k and below
  DO i = k, n
     maxv = MAXVAL(ABS(A(i,:)))
     A(i,:) = A(i,:) / maxv
     B(i) = B(i) / maxv
  ENDDO

  ! find position of max element in k column below the diagonal.
  i = k 
  DO j = k+1, n
     IF ( ABS(a(i,k)) < ABS(a(j,k)) ) i = j
  ENDDO
  ! in i column a(i,k) is maximum.

  IF (i /= k) THEN              ! exchange k, i rows of A and B
     CALL swap(a(k,:), a(i,:))
     CALL swap(b(k), b(i))     
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


SUBROUTINE nonlinear(x, f)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: x(:)
  DOUBLE PRECISION, INTENT (out) :: f(:)

  f(1) = x(1) + x(2) + x(3) - 3d0
  f(2) = x(1)**2 * x(2) + x(2)**2 * x(3) + x(3)**2 * x(1) - 4d0
  f(3) = x(1)**2  + x(2)**2  + x(3)**2  - 5d0
END SUBROUTINE nonlinear



SUBROUTINE linear(x, a)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: x(:)
  DOUBLE PRECISION, INTENT (out) :: a(:,:)

  a(1,1) = 1d0
  a(1,2) = 1d0
  a(1,3) = 1d0

  a(2,1) = 2d0 * x(1) * x(2) + x(3)**2
  a(2,2) = 2d0 * x(2) * x(3) + x(1)**2
  a(2,3) = 2d0 * x(3) * x(1) + x(2)**2

  a(3,1) = 2d0 * x(1)
  a(3,2) = 2d0 * x(2)
  a(3,3) = 2d0 * x(3)

END SUBROUTINE linear


PROGRAM solve
  IMPLICIT NONE

  INTERFACE
     SUBROUTINE nonlinear(x, f)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: x(:)
       DOUBLE PRECISION, INTENT (out) :: f(:)
     END SUBROUTINE nonlinear

     SUBROUTINE linear(x, a)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: x(:)
       DOUBLE PRECISION, INTENT (out) :: a(:,:)

     END SUBROUTINE linear

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
  
  DOUBLE PRECISION :: a(3,3), x(3), b(3), y(3)
 
  x(1) = 1.1d0
  x(2) = 2.1d0
  x(3) = 3.1d0

  DO 
     CALL nonlinear(x,b)
     IF (ALL(ABS(b) < 1d-7)) EXIT

     CALL linear(x,a)

     CALL TRIANG(a,b)
     CALL BACKSU(a,b, y)
     x = x - y
  ENDDO

  PRINT *, x, b

END PROGRAM solve
