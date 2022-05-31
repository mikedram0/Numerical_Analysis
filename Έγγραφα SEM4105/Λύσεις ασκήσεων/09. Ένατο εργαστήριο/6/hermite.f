      PROGRAM GAUHER
      INTEGER N, I
      PARAMETER (N=4)
      DOUBLE PRECISION X(N), W, FVAL, SUM
      DOUBLE PRECISION PI
      PARAMETER (PI=3.14159265358979323846D0)
      EXTERNAL F, WEIGHT
      DOUBLE PRECISION F, WEIGHT

      X(1) = SQRT ( (3.0D0-SQRT(6.0D0)) / 2.0D0 )
      X(2) = -X(1)
      X(3) = SQRT ( (3.0D0+SQRT(6.0D0)) / 2.0D0 )
      X(4) = -X(3) 
      
      SUM = 0.0D0
      DO I=1,N
        W = WEIGHT(N, X(I))
        FVAL = F(X(I))
        SUM = SUM + W * FVAL
      ENDDO

      PRINT *, "The integral is", SUM
      PRINT *, "The correct value is", SQRT(PI) / 2.0D0
      PRINT *, "Exactly the same (Why?)"
      END

      FUNCTION F(X)
      DOUBLE PRECISION X,F
      F = X*X
      END

      FUNCTION PARAG(N)
      INTEGER N, PARAG
      INTEGER TEMP, I

      TEMP = 1
      DO I=2,N
         TEMP = TEMP * I
      ENDDO
      PARAG = TEMP
      END

      FUNCTION WEIGHT(N, X)
      INTEGER N
      DOUBLE PRECISION X, WEIGHT
      EXTERNAL HERMIT
      DOUBLE PRECISION HERMIT
      DOUBLE PRECISION PI
      PARAMETER (PI=3.14159265358979323846D0)
      EXTERNAL PARAG
      INTEGER PARAG
  
      WEIGHT = 2.0D0**(N-1) * PARAG(N) * SQRT(PI) /
     &     (N * HERMIT(N-1,X))**2

      END

      FUNCTION HERMIT(N, X)
      INTEGER N
      DOUBLE PRECISION X, HERMIT

      DOUBLE PRECISION RES
      
      IF (N.EQ.0) RES = 1.0D0
      IF (N.EQ.1) RES = 2.0D0 * X
      IF (N.EQ.2) RES = 4.0D0 * X * X - 2.0D0
      IF (N.EQ.3) RES = X * (8.0D0 * X * X - 12.0D0)
      IF (N.EQ.4) RES = X * X * (16.0D0 * X * X - 48.0D0) + 12.0D0
      IF (N.GE.5) STOP "NOT IMPLEMENTED"

      HERMIT = RES

      END
