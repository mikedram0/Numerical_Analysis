      PROGRAM LEASQU
      INTEGER N
      DOUBLE PRECISION A1, A0, X(1000), Y(1000), R2

      INTEGER METHOD

      CALL READAT("points.dat", N, X, Y)

      PRINT *, "Μέθοδος: "
      PRINT *, "1) y = a * x + b"
      PRINT *, "2) y = a * x^b"
      PRINT *, "3) y = a + b * exp(x)"
      PRINT *, "4) y = a + b * ln(x)"
      PRINT *, "Δώσε μέθοδο."

      READ (*,*) METHOD
      
      IF (METHOD .EQ. 1) THEN 
         CALL LSQLIN(N, X, Y, A1, A0, R2)
         PRINT *, "y = ", A1, " x + ", A0
      ENDIF

      IF (METHOD .EQ. 2) THEN 
         CALL LSQPOW(N, X, Y, A1, A0, R2)
         PRINT *, "y = ", A0, " * x^", A1
      ENDIF

      IF (METHOD .EQ. 3) THEN 
         CALL LSQEXP(N, X, Y, A1, A0, R2)
         PRINT *, "y = ", A0, " + ", A1, " * exp(x)"
      ENDIF

      IF (METHOD .EQ. 4) THEN 
         CALL LSQLOG(N, X, Y, A1, A0, R2)
         PRINT *, "y = ", A0, " + ", A1, " * ln(x)"
      ENDIF
         
      PRINT *, "r^2 = ", R2

      END

C     Read data file
      SUBROUTINE READAT(FNAME, N, X, Y)
      CHARACTER(*) FNAME
      INTEGER N
      DOUBLE PRECISION X(*), Y(*)

      INTEGER UNITNO, I
      PARAMETER (UNITNO = 54)

      OPEN (UNIT = UNITNO, FILE = FNAME, STATUS = 'OLD')
      READ(UNITNO, *) N

      DO I=1,N
        READ (UNITNO, *) X(I), Y(I)
      ENDDO

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
      A0 = (SY - A1 * SX) / N

      R2 = (N * SXY - SX * SY)**2 / D / (N * SY2 - SY**2)
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


C     Y = A0 + A1 * EXP(X)
      SUBROUTINE LSQEXP(N, X, Y, A1, A0, R2)
      INTEGER N
      DOUBLE PRECISION X(*), Y(*), A1, A0, R2

      DOUBLE PRECISION ZX(N)
      INTEGER I

      DO 10 I = 1, N
        ZX(I) = EXP(X(I))
 10   CONTINUE
      
      CALL LSQLIN(N, ZX, Y, A1, A0, R2)

      END


C     Y = A0 + A1 * LOG(X)
      SUBROUTINE LSQLOG(N, X, Y, A1, A0, R2)
      INTEGER N
      DOUBLE PRECISION X(*), Y(*), A1, A0, R2

      DOUBLE PRECISION ZX(N)
      INTEGER I

      DO 10 I = 1, N
        IF ( (X(I) .LE. 0.0D0) ) THEN
           PRINT *, "Data not appropriate"
           RETURN
        ENDIF

        ZX(I) = LOG(X(I))
 10   CONTINUE
      
      CALL LSQLIN(N, ZX, Y, A1, A0, R2)
      
      END
