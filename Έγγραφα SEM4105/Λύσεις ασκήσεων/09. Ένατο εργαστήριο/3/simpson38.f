      PROGRAM SI1338
      DOUBLE PRECISION X(0:10000), Y(0:10000)

      EXTERNAL RDATA, FINT, SIMINT
      DOUBLE PRECISION FINT, SIMINT
      INTEGER N

      CALL RDATA("points.dat", X, Y, N)

      PRINT *, SIMINT(X, Y, N), FINT(3.0D0) - FINT(0.0D0)

      END

      FUNCTION FINT(X)
      DOUBLE PRECISION X, FINT
      FINT = (SIN(X) - COS(X)) / 2.0D0 * EXP(X)
      END


C     READ DATA FILE
      SUBROUTINE RDATA(FNAME, X, Y, N)
      CHARACTER*100 FNAME
      DOUBLE PRECISION  X(0:*), Y(0:*)
      INTEGER N

      INTEGER I, UNITNO
      PARAMETER (UNITNO = 54)

      OPEN (UNIT = UNITNO, FILE = FNAME, STATUS = 'OLD')
      READ (UNITNO, *) N
      N = N-1
      DO 10 I=0,N
        READ (UNITNO, *) X(I), Y(I)
 10   CONTINUE

      CLOSE (UNITNO)
      END


      FUNCTION SIMINT(X, Y, N)
      INTEGER N
      DOUBLE PRECISION X(0:N), Y(0:N), SIMINT
      EXTERNAL SIMP38, SIMP13
      DOUBLE PRECISION SIMP38, SIMP13

      IF (MOD(N,2) .EQ. 1) THEN 
         SIMINT = SIMP38(X(0), X(3), Y) + SIMP13(X(3), X(N), Y(3), N-3) 
      ELSE
         SIMINT = SIMP13(X(0), X(N), Y, N) 
      ENDIF

      END


      FUNCTION SIMP13(A, B, Y, N)
      INTEGER N
      DOUBLE PRECISION A, B, Y(0:N), SIMP13

      DOUBLE PRECISION SUM
      INTEGER I

      IF (MOD(N, 2) .NE. 0) STOP "ODD SIZE IN SIMP13"

      SIMP13 = Y(0) + Y(N)

      SUM = 0.0D0
      DO 10 I = 1,N-1,2
        SUM = SUM + Y(I)
 10   CONTINUE

      SIMP13 = SIMP13 + 4.0D0 * SUM

      SUM = 0.0D0
      DO 20 I = 2,N-2,2
        SUM = SUM + Y(I)
 20   CONTINUE

      SIMP13 = SIMP13 + 2.0D0 * SUM

      SIMP13 = ((B-A) / N) * SIMP13 / 3.0D0
      END 


      FUNCTION SIMP38(A, B, Y)
      DOUBLE PRECISION A, B, Y(0:3), SIMP38

      DOUBLE PRECISION  H

      H = (B-A) / 3D0

      SIMP38 = 3.0D0 * H * (Y(0) + 3D0 * (Y(1) + Y(2)) + Y(3)) / 8D0

      END



