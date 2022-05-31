      PROGRAM STEFBO
      INTEGER N
      DOUBLE PRECISION A1, A0, X(1000), Y(1000), R2

      DOUBLE PRECISION EMBADO

      CALL READAT("stefbol.dat", N, X, Y)

C     Y = A0 X^A1
      CALL LSQPOW(N, X, Y, A1, A0, R2)

      EMBADO = 0.05D-4

      PRINT *, "Σταθερά Stefan-Boltzmann: ", 
     &        A0 / EMBADO, "J/K^4/m^2/s"
      PRINT *, "Εκθέτης: ",  A1
      PRINT * , "r^2 = ", R2
      END


C     Read data file
      SUBROUTINE READAT(FNAME, N, X, Y)
      CHARACTER(*) FNAME
      INTEGER N
      DOUBLE PRECISION X(*), Y(*)

      INTEGER UNITNO, I
      PARAMETER (UNITNO = 54)

      OPEN (UNIT = UNITNO, FILE = FNAME)
      READ(UNITNO, *) N

      DO 10 I=1,N
        READ (UNITNO, *) X(I), Y(I)
 10   CONTINUE

      CLOSE (UNITNO)
      END


C     Y = A1 * X + A0
      SUBROUTINE LSQLIN(N, X, Y, A1, A0, R2)
      INTEGER N
      DOUBLE PRECISION X(*), Y(*), A1, A0, R2

      DOUBLE PRECISION SX, SY, SX2, SY2, SXY, D
      INTEGER I

      SX  = 0.0D0
      SY  = 0.0D0
      SX2 = 0.0D0
      SY2 = 0.0D0
      SXY = 0.0D0

      DO 10 I = 1, N
        SX  = SX  + X(I)
        SY  = SY  + Y(I)
        SX2 = SX2 + X(I) * X(I)
        SY2 = SY2 + Y(I) * Y(I)
        SXY = SXY + X(I) * Y(I)
 10   CONTINUE

      D = N * SX2 - SX * SX

      A1 = (N * SXY - SX * SY) / D 
      A0 = (SX2 * SY - SX * SXY) / D 

      R2 = (N * SXY - SX * SY)**2 / (N * SX2 - SX**2)/(N * SY2 - SY**2)
      END


C     Y = A0 X^A1 => LN(Y) = LN(A0) + A1 * LN(X)
      SUBROUTINE LSQPOW(N, X, Y, A1, A0, R2)
      INTEGER N     
      DOUBLE PRECISION X(*), Y(*), A1, A0, R2

      DOUBLE PRECISION ZX(N), ZY(N), ZA0
      INTEGER I

      DO 10 I = 1, N
        IF ( (X(I) .LE. 0.D0) .OR. (Y(I) .LE. 0.D0) ) THEN
           PRINT *, "Data not appropriate"
           RETURN
        ENDIF

        ZX(I) = LOG(X(I))
        ZY(I) = LOG(Y(I))
 10   CONTINUE
      
      CALL LSQLIN(N, ZX, ZY, A1, ZA0, R2)

      A0 = EXP(ZA0)

      END

