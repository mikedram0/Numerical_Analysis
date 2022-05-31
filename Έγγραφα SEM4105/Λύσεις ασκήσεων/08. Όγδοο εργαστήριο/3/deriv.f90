PROGRAM derivative
  IMPLICIT NONE
  DOUBLE PRECISION, ALLOCATABLE :: xp(:), yp(:)
  DOUBLE PRECISION :: x, der1, der2

  INTERFACE
     SUBROUTINE read_data(fname,  x, y)
       IMPLICIT NONE
       CHARACTER(*), INTENT (in) :: fname
       DOUBLE PRECISION, ALLOCATABLE, INTENT (out) :: x(:), y(:)
     END SUBROUTINE read_data

     FUNCTION deriv(order, xp, yp, x)
       IMPLICIT NONE
       INTEGER,          INTENT (in) :: order
       DOUBLE PRECISION, INTENT (in) :: xp(:), yp(:)
       DOUBLE PRECISION, INTENT (in) :: x
       DOUBLE PRECISION :: deriv
     END FUNCTION deriv

  END INTERFACE

  CALL read_data("points.dat", xp, yp)

  PRINT *, "Δώσε το σημείο υπολογισμού των παραγώγων: "
  READ *, x

  der1 = deriv(1, xp, yp, x)
  der2 = deriv(2, xp, yp, x)

  PRINT *, "Η πρώτη παράγωγος είναι", der1
  PRINT *, "Η διαφορά από την ακριβή τιμή είναι ", der1 - SIN(2d0*x)
  PRINT *
  PRINT *, "Η δεύτερη παράγωγος είναι", der2
  PRINT *, "Η διαφορά από την ακριβή τιμή είναι ", der2 - 2d0*COS(2d0*x)
END PROGRAM derivative


FUNCTION deriv(order, xp, yp, x)

  IMPLICIT NONE
  INTEGER,          INTENT (in) :: order
  DOUBLE PRECISION, INTENT (in) :: xp(:), yp(:)
  DOUBLE PRECISION, INTENT (in) :: x
  DOUBLE PRECISION :: deriv

  DOUBLE PRECISION, ALLOCATABLE :: A(:,:), b(:), solution(:)
  INTEGER :: n, k

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

  n = SIZE(xp)

  ALLOCATE(A(n,n), b(n))

  SELECT CASE (order)
  CASE (1)
     b(1) = 0.0d0   ! πρώτη παράγωγος του 1
     b(2) = 1.0d0   ! πρώτη παράγωγος του x     
     DO k=3,n
        b(k)   = (k-1) * x**(k-2)  ! πρώτη παράγωγος του x**(k-1)
     ENDDO
  CASE (2) 
     b(1) = 0.0d0   ! δεύτερη παράγωγος του 1
     b(2) = 0.0d0   ! δεύτερη παράγωγος του x     
     b(3) = 2.0d0   ! δεύτερη παράγωγος του x^2

     DO k=4,n
        b(k) = (k-1) * (k-2) * x**(k-3)  ! δεύτερη παράγωγος του x**(k-1)
     ENDDO
  CASE default
     STOP 'not implemented'
  END SELECT

  a(1,:) = 1.0d0
  DO k=2,n
     a(k,:) = xp**(k-1)
  ENDDO

  CALL TRIANG(A, B)

  ALLOCATE (solution(n))
  CALL BACKSU(A, B, solution)

  deriv = DOT_PRODUCT(solution, yp)

END FUNCTION deriv


!     Read data file
SUBROUTINE read_data(fname, x, y)
  IMPLICIT NONE
  CHARACTER(*), INTENT (in) :: fname
  DOUBLE PRECISION, ALLOCATABLE, INTENT (out) :: x(:), y(:)

  INTEGER, PARAMETER ::  UNITNO = 54
  INTEGER :: n, i

  OPEN (UNIT = UNITNO, FILE = fname, status = 'old', action = 'read')
  READ (UNITNO, *) n

  ALLOCATE(X(n), Y(n))

  DO I=1,n
     READ (UNITNO, *) X(I), Y(I)
  ENDDO

  CLOSE (UNITNO)
END SUBROUTINE read_data


!!!! Gauss 


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

