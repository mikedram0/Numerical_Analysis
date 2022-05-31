! R(x_i) = y_i => P(x_i) = y_i Q(x_i) =>
! a_0 + a_1 x_i + a_2 x_i^2 = y_i (1 + b_1 x_i + b_2 x_i^2 + b_3 x_i^3 =>
! 1 a_0 + x_i a_1 + x_i^2 a_2 - y_i x_i b_1 - y_i x_i^2 b_2 - y_i x_i^3 b_3 = y_i
!
!  | 1  x_1   x_1^2   -y_1*x_1   -y_1*x_1^2  -y_1*x_1^3 |  | a_0 |    | y_1 |
!  | 1  x_2   x_2^2   -y_2*x_2   -y_2*x_2^2  -y_2*x_2^3 |  | a_1 |    | y_2 |
!  | 1  x_3   x_3^2   -y_3*x_3   -y_3*x_3^2  -y_3*x_3^3 |  | a_2 |    | y_3 |
!  | 1  x_4   x_4^2   -y_4*x_4   -y_4*x_4^2  -y_4*x_4^3 |* | b_1 | =  | y_4 |
!  | 1  x_5   x_5^2   -y_5*x_5   -y_5*x_5^2  -y_5*x_5^3 |  | b_2 |    | y_5 |
!  | 1  x_6   x_6^2   -y_6*x_6   -y_6*x_6^2  -y_6*x_6^3 |  | b_3 |    | y_6 |
!

PROGRAM rational
  IMPLICIT NONE

  INTEGER, PARAMETER :: n = 6
  DOUBLE PRECISION :: A(n,n), B(n), X(n), Y(n), C(n)

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

     FUNCTION R(C, x)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: C(6)
       DOUBLE PRECISION, INTENT (in) :: x
       DOUBLE PRECISION :: R
     END FUNCTION R
  END INTERFACE

  INTEGER :: i
  
  X = (/ 0.9d0, 1.1d0, 1.5d0, 2.0d0, 2.9d0, 3.5d0 /)
  Y = (/ 5.607d0, 4.576d0, 3.726d0, 3.354d0, 3.14d0, 3.087d0 /)

  A(:,1) = 1d0
  A(:,2) = X
  A(:,3) = X*X
  A(:,4) = -Y*X
  A(:,5) = -Y*X*X
  A(:,6) = -Y*X*X*X

  B = Y
  
  ! solve A*C=B to get C, the coefficients a,b
  
  CALL triang(A,B)
  CALL backsu(A,B,C)

  PRINT *, "a_0  a_1  a_2 = ", C(1:3)
  PRINT *, "b_1  b_2  b_3 = ", C(4:6)


  DO i=1,n
     PRINT *, R(C, X(i)), Y(i)
  ENDDO

END PROGRAM rational



FUNCTION R(C, x)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: C(6)
  DOUBLE PRECISION, INTENT (in) :: x
  DOUBLE PRECISION :: R
  
  DOUBLE PRECISION :: p, q

  p = C(1) + C(2) * x + C(3) *x*x
  q = 1d0 + C(4) * x + C(5) * x*x + C(6) * x*x*x

  R = p/q
  
END FUNCTION R



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


  DO K = 1, N-1
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
     maxv = MAX(maxv, ABS(b(i)))
     
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

