      PROGRAM EIGEN

      INTEGER N
      PARAMETER  (N = 4)
      DOUBLE PRECISION  A(N,N)
      DOUBLE PRECISION  X1, X2, X

      EXTERNAL SECANT, F
      DOUBLE PRECISION F

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
      
      X1 = -1.0D0
      X2 =  1.0D0
      
      CALL SECANT(X1, X2, X, 1D-8, F, A, N)

      PRINT *, "AN EIGENVALUE IS ",  X, F(N,A,X)
      END 


      SUBROUTINE SECANT(X1, X2, X, TOLER, FUNC, A, N)
      INTEGER N
      DOUBLE PRECISION X1, X2, X, TOLER, A(N,N)

      EXTERNAL FUNC
      DOUBLE PRECISION  F1, F2, FUNC

 10   F1 = FUNC(N, A, X1)
      F2 = FUNC(N, A, X2)

      X = X2 - F2 * (X2 - X1) / (F2 - F1)

      IF (ABS(FUNC(N, A, X)) .LT. TOLER) RETURN

      X1 = X2
      X2 = X
      GOTO 10

      END



C GIVES THE DETERMINANT OF A
      FUNCTION DET(N,A)
      INTEGER N
      DOUBLE PRECISION A(N,N), DET

      DOUBLE PRECISION B(N,N)
      INTEGER I, J, CHANGE
      
      DO I=1,N
        DO J =1, N
          B(I,J) = A(I,J)
        ENDDO
      ENDDO
      
      CHANGE = 0

      CALL TRIANG(N,B,CHANGE)
C     NOW B IS DIAGONAL
      
      DET = 1.0D0
      DO I= 1,N
        DET = DET * B(I,I)
      ENDDO

      IF (MOD(CHANGE,2).NE.0) DET = -DET

      END

C     MAKE  A  LOWER TRIANGULAR
      SUBROUTINE TRIANG(N, A, CHANGE)
      INTEGER N, CHANGE
      DOUBLE PRECISION A(N,N)
      EXTERNAL PIVOT
      
      DOUBLE PRECISION  G
      INTEGER I,J, K

C     LOWER TRIANGLE
      DO I=1,N-1
        CALL PIVOT(N, A, I, CHANGE)
        DO K=I+1,N
          G = -A(K,I) / A(I,I)
          DO J=I,N
            A(K,J) = A(K,J) + A(I,J) * G
          ENDDO
        ENDDO
      ENDDO
      END 
      
      
C     MERIKI ODIGISI (PARTIAL PIVOTING)
      SUBROUTINE PIVOT(N, A, I, CHANGE)
      INTEGER N, CHANGE
      DOUBLE PRECISION A(N,N)
      INTEGER I
      EXTERNAL SWAP
      
      INTEGER J, K
      
      J = I
      DO K = I+1, N
        IF ( ABS(A(K,I)) .GT. ABS(A(J,I)) ) J = K
      ENDDO
C     IN I COLUMN A(J,I) IS MAXIMUM.
      
      IF (J .GT. I) THEN 
         CHANGE = CHANGE + 1
C     EXCHANGE I, J ROWS
         DO K = 1,N
           CALL SWAP(A(J,K), A(I,K)) 
         ENDDO
      ENDIF
      END 

C     A<->B
      SUBROUTINE SWAP(A,B)
      DOUBLE PRECISION A, B
      
      DOUBLE PRECISION TEMP

      TEMP = A
      A    = B
      B    = TEMP
      
      END



C FUNCTION F(X) WHERE X IS LAMBDA
      FUNCTION F(N,A,LAMBDA)
      INTEGER N
      DOUBLE PRECISION A(N,N), LAMBDA, F
      DOUBLE PRECISION DET
      EXTERNAL DET

      DOUBLE PRECISION  B(N,N)
      INTEGER I, J

C     SET B = A - LAMBDA I
      DO I=1,N
        DO J = 1, N
          B(I,J) = A(I,J)
        ENDDO
        B(I,I) = A(I,I) - LAMBDA
      ENDDO


      F = DET(N,B)

      END 

