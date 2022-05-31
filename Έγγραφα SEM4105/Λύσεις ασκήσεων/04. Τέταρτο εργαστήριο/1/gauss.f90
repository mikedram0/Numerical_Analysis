PROGRAM GAUSS
  IMPLICIT NONE
  DOUBLE PRECISION, ALLOCATABLE :: A(:,:), B(:), X(:)

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

     SUBROUTINE GET_MATRICES(FNAME, A, B)
       IMPLICIT NONE
       CHARACTER(*), INTENT (in) :: FNAME
       DOUBLE PRECISION, ALLOCATABLE, INTENT (out) :: A(:,:), B(:)
     END SUBROUTINE GET_MATRICES

     SUBROUTINE PRINT_MATRIX(A)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: A(:,:)
     END SUBROUTINE PRINT_MATRIX
  END INTERFACE

  CALL GET_MATRICES("gauss1.txt", A, B)
  CALL TRIANG(A, B)
  CALL PRINT_MATRIX(A)

  ALLOCATE (X(SIZE(b)))
  CALL BACKSU(A, B, X)

  PRINT *, X

  DEALLOCATE(a,b,x)
END PROGRAM GAUSS




! Κάτω τριγωνοποίηση: Στη στήλη Κ (από 1 έως Ν-1) μηδενίζουμε τα στοιχεία 
! του Α κάτω από τη διαγώνιο (δηλ. στη γραμμή Ι από Κ+1 έως το τέλος). 
! Αυτό γίνεται προσθέτοντας τα στοχεία της γραμμής Κ πολλαπλασιασμένα 
! με κατάλληλη σταθερά (G), που αλλάζει σε κάθε γραμμή.
!
! Την ίδια μετατροπή κάνουμε και στο διάνυσμα Β (με την ίδια σταθερά G 
! όπως και στον Α).

SUBROUTINE TRIANG(A, B)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (inout) :: A(:,:), B(:)

  INTEGER :: N, I, K
  DOUBLE PRECISION :: G

  N = SIZE(B)

  DO K=1, N-1
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

  INTEGER :: I, J, N
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
        WRITE (*,"(f15.7)", advance = "no") A(I,J)
     ENDDO
     WRITE (*,*)
  ENDDO
  WRITE (*,*)
END SUBROUTINE PRINT_MATRIX


SUBROUTINE BACKSU(A, B, X)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in)  :: A(:,:), B(:)
  DOUBLE PRECISION, INTENT (out) :: X(:)

  INTEGER :: I, N

  N = SIZE(B)

  ! Οπισθοδρόμηση: Ξεκινώντας από το τέλος προς την αρχή, υπολογίζουμε 
  ! τα στοιχεία του άγνωστου διανύσματος Χ.

  DO I = N,1,-1
     X(I) = (B(I) - DOT_PRODUCT(A(I,I+1:N), X(I+1:N)) ) / A(I,I)
  ENDDO

END SUBROUTINE BACKSU
