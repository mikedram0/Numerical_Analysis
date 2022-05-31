      PROGRAM GAUSS
      DOUBLE PRECISION A(100,100), B(100), X(100)
      INTEGER N, I
      EXTERNAL GETMAT, PRIMAT, TRIANG, BACKSU
 
      CALL GETMAT("gauss1.txt", N, A, B)
      CALL PRIMAT(N, A)
      CALL TRIANG(N, A, B)
      CALL PRIMAT(N, A)
      CALL BACKSU(N, A, B, X)

      DO 10 I=1,N
        PRINT *, X(I)
 10   CONTINUE

      END


C     READ A,B FROM FILE 
      SUBROUTINE GETMAT(FNAME, N, A, B)
      CHARACTER(*) FNAME
      INTEGER N
      DOUBLE PRECISION A(100,100), B(100)

      INTEGER I, J
      INTEGER UNITNO
      PARAMETER (UNITNO=12)

      OPEN (UNITNO,FILE=FNAME)
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
        WRITE(*, "(100F15.7)") (A(I,J), J=1,N)
 10   CONTINUE
      WRITE (*,*)
      END

C     MAKE LOWER TRIANGULAR
      SUBROUTINE TRIANG(N, A, B)
      INTEGER N
      DOUBLE PRECISION A(100, 100), B(100)

      INTEGER I, J, K
      DOUBLE PRECISION G

      DO 10 K=1,N-1
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
      INTEGER K, J, N
      DOUBLE PRECISION A(100,100), B(100), X(100)
      DOUBLE PRECISION SUM
      
      DO 10 K=N,1,-1
        SUM = 0.0D0
        DO 20 J=K+1,N
          SUM = SUM + A(K,J) * X(J)
 20     CONTINUE
        X(K) = (B(K) - SUM) / A(K,K)
 10   CONTINUE

      END
