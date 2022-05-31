C     Epilisi tis f(x)=0 me f(x)=3x e^x -1
C     me th methodo Newton-Raphson.

      PROGRAM NEWTRA

      DOUBLE PRECISION X, F, DF
C     declare that F, DF are functions
      EXTERNAL F, DF

C     A small constant number
      DOUBLE PRECISION TOLER
      PARAMETER (TOLER = 1D-8)

C     Initial approximation 
      X = 1.0D0

 10   X = X - F(X) / DF(X)
      IF (ABS(F(X)) .GT. TOLER) GOTO 10

      PRINT *, "The root is approximately ", X
      PRINT *, "The function value is ", F(X)

      END

C     The function
      FUNCTION F(X)
      DOUBLE PRECISION X, F

      F = 3.0D0 * X * EXP(X) - 1.0D0
      END

C     The first derivative 
      FUNCTION DF(X)
      DOUBLE PRECISION X, DF

      DF = 3.0D0 * EXP(X) * (1.0D0 + X)
      END
