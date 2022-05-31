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
      DO 10 I= 1, N
        X = A + STEP * (I-1)
        Y = SIN(X)

        WRITE(22,*) X, Y
 10   CONTINUE

      CLOSE(22)
      END

      SUBROUTINE READDA(FNAME, X, Y, N)
      CHARACTER(*) FNAME
      DOUBLE PRECISION X(0:*), Y(0:*)
      INTEGER N

      INTEGER I
      
      OPEN(22, FILE=FNAME)
      
      READ (22, *) N

      N = N-1
      
      DO 10 I = 0, N
        READ (22,*) X(I), Y(I)
 10   CONTINUE
      
      CLOSE(22)
      END 

      FUNCTION Q(I,X,Z)
      INTEGER I
      DOUBLE PRECISION X(0:*), Z, Q

      INTEGER J
      
      Q = 1D0
      DO 10 J = 0, I-1
        Q = Q * (Z-X(J))
 10   CONTINUE

      END
      
      
      SUBROUTINE COEFS(N, X, Y, A)
      INTEGER N      
      DOUBLE PRECISION X(0:*), Y(0:*), A(0:*)

      DOUBLE PRECISION S, Q
      INTEGER I
      EXTERNAL Q
      
      DO 20 I = 0, N        
        S = 0D0
        DO 10 J = 0, I-1
          S = S + A(J) * Q(J, X, X(I)) 
 10     CONTINUE
        
        A(I) = (Y(I) - S) / Q(I, X, X(I))
 20   CONTINUE
      
      END 


      FUNCTION POLY(N, A, X, Z)
      INTEGER N
      DOUBLE PRECISION A(0:*), X(0:*), Z
      DOUBLE PRECISION POLY

      DOUBLE PRECISION Q
      EXTERNAL Q

      INTEGER I

      POLY = 0D0
      DO 10 I = 0, N
        POLY = POLY + A(I) * Q(I, X, Z)
 10   CONTINUE
      
      END





      FUNCTION MAXVAL(N,X)
      INTEGER N
      DOUBLE PRECISION X(0:*), MAXVAL

      DOUBLE PRECISION P
      INTEGER I

      P=X(0)
      DO 10 I=1,N
        IF (P .LT. X(I)) P = X(I)
 10   CONTINUE
      
      MAXVAL = P
      END

      FUNCTION MINVAL(N,X)
      INTEGER N
      DOUBLE PRECISION X(0:*), MINVAL

      DOUBLE PRECISION P
      INTEGER I

      P=X(0)
      DO 10 I=1,N
        IF (P .GT. X(I)) P = X(I)
 10   CONTINUE

      MINVAL = P
      END



      PROGRAM NEWTON
      DOUBLE PRECISION X(1000), Y(1000), A(1000)
      INTEGER M
      PARAMETER(M = 100)
      DOUBLE PRECISION MINX, MAXX, STEP, XOUT, YOUT
      INTEGER I, N
      DOUBLE PRECISION POLY, MAXVAL, MINVAL
      EXTERNAL MAXVAL, MINVAL, CREATE, READDA, POLY

      CALL CREATE("points.dat")
      CALL READDA("points.dat", X, Y, N)
      
      CALL COEFS(N,X,Y,A)


      MINX = MINVAL(N,X)
      MAXX = MAXVAL(N,X)
      STEP = (MAXX - MINX) / (M-1)

      DO 10 I = 1,M
        XOUT = MINX + STEP * (I-1)
        YOUT = POLY(N, A, X, XOUT)

        PRINT *, XOUT, YOUT
 10   CONTINUE

      END 

