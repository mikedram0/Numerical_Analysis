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
      

      FUNCTION ELL(I, Z, N, X)
      INTEGER I, N
      DOUBLE PRECISION Z, X(N)
      DOUBLE PRECISION ELL
      
      DOUBLE PRECISION P
      INTEGER  J

      P = 1.0D0
      DO J = 1, SIZE(X)
        IF (I/=J) P = P * (Z-X(J)) / (X(I)-X(J))
      ENDDO

      ELL = P 

      END


      FUNCTION POLY(N, X, Y, Z)
      INTEGER N
      DOUBLE PRECISION X(N), Y(N), Z
      DOUBLE PRECISION POLY

      INTEGER I
      DOUBLE PRECISION P, ELL
      EXTERNAL ELL

      P = 0.0D0
      DO I= 1, N
        P = P + ELL(I,Z,N,X) * Y(I)
      END DO

      POLY = P
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

      PROGRAM LAGRAN
      DOUBLE PRECISION X(1000), Y(1000)
      INTEGER M
      PARAMETER(M = 100)
      DOUBLE PRECISION MIN_X, MAX_X, STEP, XOUT, YOUT
      INTEGER I, N
      DOUBLE PRECISION POLY, MAXVAL, MINVAL
      EXTERNAL MAXVAL, MINVAL, CREATE, READDA, POLY

      CALL CREATE("points.dat")
      CALL READDA("points.dat", X, Y, N)
      
      MIN_X = MINVAL(N,X)
      MAX_X = MAXVAL(N,X)
      STEP = (MAX_X - MIN_X) / (M-1)

      DO I = 1,M
        XOUT = MIN_X + STEP * (I-1)
        YOUT = POLY(N, X, Y, XOUT)

        PRINT *, XOUT, YOUT
      ENDDO

      END 
