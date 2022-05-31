      SUBROUTINE SECANT(X1, X2, X, TOLER, FUNC)

      DOUBLE PRECISION X1, X2
      DOUBLE PRECISION X, TOLER
      DOUBLE PRECISION FUNC
      EXTERNAL FUNC

      DOUBLE PRECISION F1, F2

      F1 = FUNC(X1)
      F2 = FUNC(X2)

 10   X = X2 - F2 * (X2 - X1) / (F2 - F1)

      X1 = X2
      F1 = F2

      X2 = X
      F2 = FUNC(X)
      
      IF (ABS(F2) .GT. TOLER) GOTO 10

      END


C     The function
      FUNCTION F(X)
      DOUBLE PRECISION X, F

      F = 3.0D0 * LOG(X) + 5.0D0
      END

      PROGRAM SEC
C     /////   SECANT   ///////

C     A small constant number
      DOUBLE PRECISION TOLER
      PARAMETER (TOLER = 1D-8)

      DOUBLE PRECISION X, F
      EXTERNAL F
      
      DOUBLE PRECISION A,B
C     Initial approximations
      A = 0.1D0
      B = 0.2D0
      
      CALL SECANT(A, B, X, TOLER, F)

      PRINT *, "A root is approximately ", X
      PRINT *, "The function value is ", F(X)

      END
