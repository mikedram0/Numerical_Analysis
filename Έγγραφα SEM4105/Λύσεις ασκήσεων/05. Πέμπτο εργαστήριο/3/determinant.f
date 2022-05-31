C     A<->B
      SUBROUTINE SWAP(A,B)
      DOUBLE PRECISION A,B
      
      DOUBLE PRECISION TEMP
      
      TEMP = A
      A    = B
      B    = TEMP
      
      END 
      

C     FIND LINE M WHERE A(K+1:N, K) IS NOT ZERO
      SUBROUTINE FIND(A, N, K, M)
      INTEGER N, K, M
      DOUBLE PRECISION A(N,N)
      
      DOUBLE PRECISION EPS
      PARAMETER (EPS = 1D-7)
      
      DO 10 M = K+1, SIZE(A,1)
        IF (ABS(A(M,K)) .GT. EPS) RETURN
 10   CONTINUE

      END 

C     MAKE  A LOWER TRIANGULAR
      SUBROUTINE TRIANG(A, N, CHANGE)
      INTEGER N, CHANGE
      DOUBLE PRECISION A(N,N)
      
      DOUBLE PRECISION  G
      INTEGER  I, K, J, M
      DOUBLE PRECISION EPS
      PARAMETER (EPS = 1D-7)

      CHANGE = 0
      
      DO 10 K=1,N-1
        DO 20 I=K+1,N
          
          IF (ABS(A(K,K)) .LT. EPS) THEN
             CALL FIND(A, N, K, M)
             DO 30 J=1,N
               CALL SWAP(A(K,J), A(M,J))
 30          CONTINUE
             CHANGE = CHANGE + 1
          ENDIF
          
          G = -A(I,K) / A(K,K)
          DO 40 J = K,N
            A(I,J) = A(I,J) + A(K,J) * G
 40       CONTINUE
 20     CONTINUE
 10   CONTINUE
      
      END
      

C     GIVES THE DETERMINANT OF A
      FUNCTION DET(A, N)
      INTEGER N
      DOUBLE PRECISION A(N,N)
      DOUBLE PRECISION DET
      EXTERNAL TRIANG
      
      DOUBLE PRECISION B(N,N)
      INTEGER I, J, CHANGE

C     COPY A TO B
      DO 10 I=1,N
        DO 20 J=1,N
          B(I,J) = A(I,J)
 20     CONTINUE
 10   CONTINUE

      CALL TRIANG(B,N, CHANGE)
C     NOW B IS LOWER TRIANGULAR

      DET = 1.0D0
      DO 30 I= 1,N
        DET = DET * B(I,I)
 30   CONTINUE

      IF (MOD(CHANGE, 2).EQ.1) DET = -DET

      END 


      PROGRAM DETER
      INTEGER N
      PARAMETER(N=4)
      DOUBLE PRECISION A(N,N), DET
      EXTERNAL DET
      
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
      
      PRINT *, DET(A, N)
      
      END
