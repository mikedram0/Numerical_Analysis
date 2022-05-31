! a<->b
ELEMENTAL SUBROUTINE swap(a,b)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (inout) :: a, b

  DOUBLE PRECISION :: temp

  temp = a
  a    = b
  b    = temp

END SUBROUTINE swap


!     MAKE LOWER TRIANGULAR
SUBROUTINE TRIANG(A, change)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (inout) :: A(:,:)
  INTEGER, INTENT (out) :: change

  INTEGER :: N, I, K, M
  DOUBLE PRECISION :: G
  DOUBLE PRECISION, PARAMETER :: eps = 1d-7

  INTERFACE 
     ELEMENTAL SUBROUTINE swap(a,b)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (inout) :: a, b
     END SUBROUTINE swap
  END INTERFACE

  N = SIZE(A,1)

  change = 0
  DO K=1, N-1
     M = MAXLOC(ABS(A(K:N,K)),1) + K-1
     
     IF (M /= K) THEN
        CALL swap(A(K,:), A(M,:))
        change = change + 1
     ENDIF
     
     DO I=K+1, N
        G = -A(I,K) / A(K,K)
        A(I,K:N) = A(I,K:N) + A(K,K:N) * G
     ENDDO
  ENDDO

END SUBROUTINE TRIANG

! gives the determinant of a
FUNCTION det(a)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: a(:,:)
  DOUBLE PRECISION :: det

  DOUBLE PRECISION :: b(SIZE(a,1),SIZE(a,2))
  INTEGER :: i, change

  INTERFACE
     SUBROUTINE triang(a, change)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (inout) :: a(:,:)
       INTEGER, INTENT (out) :: change
     END SUBROUTINE triang
  END INTERFACE

  b = a
  CALL triang(b, change) ! now b is lower triangular

  det = 1.0d0
  DO i = 1,SIZE(a,1)
     det = det * b(i,i)
  ENDDO

  IF (MOD(change, 2) == 1) det = -det
END FUNCTION det


PROGRAM deter
  IMPLICIT NONE

  INTEGER, PARAMETER :: n = 4
  DOUBLE PRECISION :: a(n,n)

  INTERFACE
     FUNCTION det(a)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: a(:,:)
       DOUBLE PRECISION :: det
     END FUNCTION det
  END INTERFACE

  A(1,1) = 2.1D0
  A(1,2) = 3.9D0
  A(1,3) = 0.3D0
  A(1,4) = -4.1D0

  A(2,1) = 4.3D0
  A(2,2) = -1.3D0
  A(2,3) = 0.8D0
  A(2,4) = 1.5D0

  A(3,1) = 1.0D0
  A(3,2) = -2.8D0
  A(3,3) = 4.3D0
  A(3,4) = -8.1D0

  A(4,1) = 2.4D0
  A(4,2) = 6.1D0
  A(4,3) = -1.1D0
  A(4,4) = 12.5D0

  PRINT *, det(a)

END PROGRAM deter
