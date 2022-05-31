      PROGRAM GAUSS
      DOUBLE PRECISION A(100,100), B(100), X(100)
      INTEGER N, I
      EXTERNAL GETMAT, PRIMAT, TRIANG, BACKSU
      
      CALL GETMAT("gauss2.txt", N, A, B)
      CALL TRIANG(N, A, B)
      CALL BACKSU(N, A, B, X)

      DO 10 I=1,N
        PRINT *, X(I)
 10   CONTINUE

      END

C     READ A,B FROM FILE 
      SUBROUTINE GETMAT(FNAME, N, A, B)
      CHARACTER(*) FNAME
      INTEGER N
      DOUBLE PRECISION A(100,100), B(*)

      INTEGER I, J
      INTEGER UNITNO
      PARAMETER (UNITNO=12)

      OPEN(UNITNO,FILE=FNAME)
      READ (UNITNO,*) N

      DO 10 I=1,N
        READ (UNITNO,*) (A(I,J), J=1,N)
 10   CONTINUE

      DO 20 I=1,N
        READ (UNITNO,*) B(I)
 20   CONTINUE

      CLOSE(UNITNO)
      END


C     PRINT MATRIX
      SUBROUTINE PRIMAT(N, A)
      INTEGER N
      DOUBLE PRECISION A(100,100)

      INTEGER I,J

      DO 10 I = 1,N
        WRITE (*, "(100F15.7)") (A(I,J), J=1,N)
 10   CONTINUE
      WRITE (*,*)
      END



C     MAKE LOWER TRIANGULAR
      SUBROUTINE TRIANG(N, A, B)
      INTEGER N
      DOUBLE PRECISION A(100,100), B(100)
      
      INTEGER I, J, K
      DOUBLE PRECISION G
      EXTERNAL PIVOT

      DO 10 K=1,N-1
C     ------------------------------------------  this is new
        CALL PIVOT(N, A, B, K)    
        IF (ABS(A(K,K)) .LT. 1D-8) GOTO 10
        
        DO 20 I=K+1,N
          G = -A(I,K) / A(K,K)
          DO 30 J = K,N
            A(I,J) = A(I,J) + A(K,J) * G
 30       CONTINUE
          B(I) = B(I) + B(K) * G
 20     CONTINUE
 10   CONTINUE

      END


C     BACKSUBSTITUTION
      SUBROUTINE BACKSU(N, A, B, X)
      INTEGER N
      DOUBLE PRECISION A(100,100), B(100), X(100)
      
      DOUBLE PRECISION SUM
      INTEGER K, J
      
      DO 10 K=N,1,-1
        SUM = 0.0D0
        DO 20 J=K+1,N
          SUM = SUM + A(K,J) * X(J)
 20     CONTINUE
        X(K) = (B(K) - SUM) / A(K,K)
 10   CONTINUE

      END

C     //// This is new ////////////////////

C     MERIKI ODIGISI (PARTIAL PIVOTING)

      SUBROUTINE PIVOT(N, A, B, K)
      INTEGER N
      DOUBLE PRECISION A(100,100), B(100)
      INTEGER K
      
      INTEGER I, J
      DOUBLE PRECISION MAXV
      EXTERNAL SWAP

C     'NORMALIZE' EQUATIONS FROM K AND BELOW
      DO 10 I = K, N

        MAXV = ABS(B(I))
        DO 20 J = 1, N
          IF (MAXV .LT. ABS(A(I,J))) MAXV = ABS(A(I,J))
 20     CONTINUE

        DO 30 J = 1, N
          A(I,J) = A(I,J) / MAXV
 30     CONTINUE
        B(I) = B(I) / MAXV

 10   CONTINUE
      

C     FIND POSITION OF MAX ELEMENT IN K COLUMN BELOW THE DIAGONAL.
      I = K
      DO 40 J = K+1, N
        IF ( ABS(A(I,K)) .GT. ABS(A(J,K)) ) I = J
 40   CONTINUE
C     IN I COLUMN A(I,K) IS MAXIMUM.
      
      IF (I .NE. K) THEN              
C     EXCHANGE I, K ROWS OF A AND B
         DO 50 J=K,N
           CALL SWAP(A(I,J), A(K,J))
 50      CONTINUE
         CALL SWAP(B(I), B(K))    
      END IF
      END 

C     A<->B
      SUBROUTINE SWAP(A,B)
      DOUBLE PRECISION A, B
      
      DOUBLE PRECISION TEMP
      
      TEMP = A
      A    = B
      B    = TEMP
      END

