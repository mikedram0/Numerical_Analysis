PROGRAM GAUSS
  IMPLICIT NONE
  DOUBLE PRECISION, ALLOCATABLE :: A(:,:), B(:), X(:)

  INTERFACE 
     SUBROUTINE TRIANG(A, B)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (inout) :: A(:,:), B(:)
     END SUBROUTINE TRIANG

     SUBROUTINE GET_MATRICES(FNAME, A, B)
       IMPLICIT NONE
       CHARACTER(*), INTENT (in) :: FNAME
       DOUBLE PRECISION, ALLOCATABLE, INTENT (out) :: A(:,:), B(:)
     END SUBROUTINE GET_MATRICES

     SUBROUTINE PRINT_MATRIX(A)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: A(:,:)
     END SUBROUTINE PRINT_MATRIX

     SUBROUTINE BACKSU(A, B, X)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in)  :: A(:,:), B(:)
       DOUBLE PRECISION, INTENT (out) :: X(:)
     END SUBROUTINE BACKSU
  END INTERFACE

  CALL GET_MATRICES("gauss2.txt", A, B)
  CALL TRIANG(A, B)

  ALLOCATE (X(SIZE(B)))
  CALL BACKSU(A, B, X)

  PRINT *, X

  DEALLOCATE(a,b,x)
END PROGRAM GAUSS



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
     CALL pivot(a, b, k)                  ! <---------- this is new
     IF (ABS(a(K,K)) < 1d-8) CYCLE
     DO I=K+1, N
        G = -A(I,K) / A(K,K)
        A(I,K:N) = A(I,K:N) + A(K,K:N) * G
        B(I) = B(I) + B(K) * G
     ENDDO
  ENDDO

END SUBROUTINE TRIANG



!     READ A,B FROM FILE 
SUBROUTINE GET_MATRICES(FNAME, A, B)
  IMPLICIT NONE
  CHARACTER(*), INTENT (in) :: FNAME
  DOUBLE PRECISION, ALLOCATABLE, INTENT (out) :: A(:,:), B(:)

  INTEGER :: I, N
  INTEGER, PARAMETER :: unitno = 12

  OPEN(unitno, FILE=FNAME)
  READ (unitno,*) N

  ALLOCATE(A(N,N), B(N))  

  DO I=1,N
     READ (unitno,*) A(I,1:N)
  ENDDO

  DO I=1,N
     READ (unitno,*) B(I)
  ENDDO

  CLOSE(unitno)
END SUBROUTINE GET_MATRICES


!     PRINT MATRIX
SUBROUTINE PRINT_MATRIX(A)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: A(:,:)

  INTEGER :: I, J

  DO I = 1,SIZE(A,1)
     DO J = 1, SIZE(A,2)
        WRITE (*,"(f15.7)",  advance = "no") A(I,J)
     ENDDO
     WRITE (*,*)
  ENDDO
  WRITE (*,*)
END SUBROUTINE PRINT_MATRIX


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


! //// This is new ////////////////////

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
     CALL swap(a(k,k:), a(i,k:))
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

