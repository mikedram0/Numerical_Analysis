C     MAKE LOWER TRIANGULAR
      SUBROUTINE TRIANG(N, A, B)
      INTEGER N
      DOUBLE PRECISION A(1000, 1000), B(1000)
      INTEGER I, J, K

      DOUBLE PRECISION G
      EXTERNAL PIVOT

      DO K=1,N-1
        CALL PIVOT(N, A, B, K)    

        DO I=K+1,N
          G = -A(I,K) / A(K,K)
          DO J = 1,N
            A(I,J) = A(I,J) + A(K,J) * G
          ENDDO
          B(I) = B(I) + B(K) * G
        ENDDO
      ENDDO

      END


C     BACKSUBSTITUTION
      SUBROUTINE BACKSU(N, A, B, X)
      INTEGER K, J, N
      DOUBLE PRECISION A(1000,1000), B(1000), X(1000)
      DOUBLE PRECISION SUM
      
      DO K=N,1,-1
        SUM = 0.0D0
        DO J=K+1,N
          SUM = SUM + A(K,J) * X(J)
        ENDDO
        X(K) = (B(K) - SUM) / A(K,K)
      ENDDO

      END

C     MERIKI ODIGISI (PARTIAL PIVOTING)
      SUBROUTINE PIVOT(N, A, B, I)
      DOUBLE PRECISION A(1000,1000), B(1000)
      INTEGER I
      INTEGER J, K

C     FIND POSITION OF MAX ELEMENT IN I COLUMN BELOW THE DIAGONAL.
      J = I
      DO K = I+1, N
        IF ( ABS(A(K,I)) > ABS(A(J,I)) ) J = K
      ENDDO
C     IN I COLUMN A(J,I) IS MAXIMUM.
      
      IF (J .GT. I) THEN              
C     EXCHANGE I, J ROWS OF A AND B
         DO K=1,N
           CALL SWAP(A(I,K), A(J,K))
         ENDDO
         CALL SWAP(B(I), B(J))     
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


      SUBROUTINE CREATE(FNAME)
      CHARACTER(*) FNAME

      INTEGER N
      PARAMETER (N = 15)
      DOUBLE PRECISION A, B, STEP
      PARAMETER(A = 2.0D0, B = 4.0D0, STEP = (B-A)/(N-1))
      DOUBLE PRECISION X, Y
      INTEGER I

      OPEN(22, FILE=FNAME)

      WRITE (22, *) N
      DO I= 1, N
        X = A + STEP * (I-1)
        Y = SIN(X)

        WRITE(22,*) X, Y
      ENDDO

      CLOSE(22)
      END

      SUBROUTINE READDA(FNAME, X, Y, N)
      CHARACTER(*) FNAME
      DOUBLE PRECISION X(*), Y(*)
      INTEGER N

      INTEGER I
      
      OPEN(22, FILE=FNAME)
      
      READ (22, *) N
      
      DO I= 1, N
        READ (22,*) X(I), Y(I)
      ENDDO
      
      CLOSE(22)
      END 


      FUNCTION MAXVAL(N,X)
      INTEGER N
      DOUBLE PRECISION X(N), MAXVAL

      DOUBLE PRECISION P
      INTEGER I

      P=X(1)
      DO I=2,N
        IF (P .LT. X(I)) P = X(I)
      ENDDO
      
      MAXVAL = P
      END

      FUNCTION MINVAL(N,X)
      INTEGER N
      DOUBLE PRECISION X(N), MINVAL

      DOUBLE PRECISION P
      INTEGER I

      P=X(1)
      DO I=2,N
        IF (P .GT. X(I)) P = X(I)
      ENDDO
      
      MINVAL = P
      END


      SUBROUTINE SPLINE(N,X,Y,A,B,C,D)
      INTEGER N
      DOUBLE PRECISION X(0:N), Y(0:N)
      DOUBLE PRECISION A(0:N), B(0:N), C(0:N), D(0:N)

      DOUBLE PRECISION  MATA(1000,1000), MATB(1000), MATX(1000)
      
      INTEGER I, J

      DO I=1,4*N
        DO J=1, 4*N
          MATA(I,J) = 0.0D0
        ENDDO
          MATB(I) = 0.0D0
      ENDDO

      DO I=0,N-1
        MATA(I+1,4*I+4) = 1.0D0
        
        MATB(I+1) = Y(I)
      ENDDO

      DO I=0,N-1
        MATA(I+1+N,4*I+1) = (X(I+1)-X(I))**3
        MATA(I+1+N,4*I+2) = (X(I+1)-X(I))**2
        MATA(I+1+N,4*I+3) = (X(I+1)-X(I))
        MATA(I+1+N,4*I+4) = 1.0D0
        
        MATB(I+1+N) = Y(I+1)
      ENDDO


      DO I=1,N-1
        MATA(I+2*N, 4*(I-1)+1) = 3.0D0 * (X(I)-X(I-1))**2
        MATA(I+2*N, 4*(I-1)+2) = 2.0D0 * (X(I)-X(I-1))
        MATA(I+2*N, 4*(I-1)+3) = 1.0D0
        
        MATA(I+2*N, 4*I+3) = -1.0D0
      ENDDO
      
      DO I=1,N-1
        MATA(I-1+3*N, 4*(I-1)+1) = 3.0D0 * (X(I)-X(I-1))
        MATA(I-1+3*N, 4*(I-1)+2) = 1.0D0
        
        MATA(I-1+3*N, 4*I+2) = -1.0D0
      ENDDO
      
      MATA(4*N-1,2) = 1.0D0
      
      MATA(4*N,4*N-3) = 3.0D0 * (X(N)-X(N-1))
      MATA(4*N,4*N-2) = 1.0D0
      

      CALL TRIANG(4*N, MATA, MATB)
      CALL BACKSU(4*N, MATA, MATB, MATX)
  
      DO I=0,N-1
        A(I) = MATX(4*I+1)
        B(I) = MATX(4*I+2)
        C(I) = MATX(4*I+3)
        D(I) = MATX(4*I+4)
      ENDDO

      END 


      FUNCTION POLY(N, X, A, B, C, D, Z)
      INTEGER N
      DOUBLE PRECISION X(0:N), A(0:N), B(0:N), C(0:N), D(0:N), Z
      DOUBLE PRECISION POLY

      INTEGER I
      DOUBLE PRECISION  DELTA

      DO I=0,N-1
        IF (Z .LE. X(I+1)) GOTO 10
      ENDDO
      
 10   DELTA = Z-X(I)
      POLY = A(I) * DELTA**3 + B(I) * DELTA**2 + C(I) * DELTA + D(I) 
      
      END



      PROGRAM SP
      INTEGER K
      PARAMETER (K=1000)
      DOUBLE PRECISION X(0:K), Y(0:K), A(0:K), B(0:K), C(0:K), D(0:K)
      INTEGER M
      PARAMETER (M = 100)
      DOUBLE PRECISION  MIN_X, MAX_X, STEP, XOUT, YOUT
      INTEGER I, N
      DOUBLE PRECISION POLY, MAXVAL, MINVAL
      EXTERNAL  MAXVAL, MINVAL
      
      CALL CREATE("points.dat")
      CALL READDA("points.dat", X, Y, N)

      CALL SPLINE(N-1,X,Y,A,B,C,D)
      
      MIN_X = MINVAL(N,X)
      MAX_X = MAXVAL(N,X)
      STEP = (MAX_X - MIN_X) / (M-1)

      DO I = 1,M
       XOUT = MIN_X + STEP * (I-1)
       YOUT = POLY(N-1,X, A,B,C,D, XOUT)

       PRINT *, XOUT, YOUT
      ENDDO

      END

