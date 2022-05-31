      FUNCTION SIMPSO(F, A, B, N)

      DOUBLE PRECISION A, B, SIMPSO
      INTEGER N
      DOUBLE PRECISION F
      EXTERNAL F

      INTEGER I
      DOUBLE PRECISION X, SUM1, SUM2, SUM
      DOUBLE PRECISION H

      IF (MOD(N,2) .NE. 0) STOP

      H = (B-A) / N

      SUM1 = 0.0D0
      DO I= 1,N-1,2
         X = A + I * H

         SUM1 = SUM1 + F(X)
      ENDDO

      SUM2 = 0.0D0
      DO I= 2,N-2,2
         X = A + I * H

         SUM2 = SUM2 + F(X)
      ENDDO


      SUM = (F(A) + F(B)) + 4.0D0 * SUM1 + 2.0D0 * SUM2

      SIMPSO = SUM * H / 3.0D0
      END


      FUNCTION F(X)
      DOUBLE PRECISION X, F
      F = SIN(X)
      END

      PROGRAM SIMPS

      DOUBLE PRECISION SIMPSO, F
      EXTERNAL SIMPSO, F

      DOUBLE PRECISION PI
      PARAMETER (PI = 3.14159265358979323846D0)
      DOUBLE PRECISION INT
      INTEGER N, K

      PRINT *, "Ν       ΤΙΜΗ                     ΣΦΑΛΜΑ" 

      DO K = 1, 9
        N = 2**K
        INT = SIMPSO(F, 0.0D0, PI, N)
        PRINT *, N, INT, 2.0D0-INT
      ENDDO

      END
