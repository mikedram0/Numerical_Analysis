      PROGRAM ASK3

      EXTERNAL F, G, REGFAL, BISECT
      DOUBLE PRECISION X, F, G
      DOUBLE PRECISION TOLER
      PARAMETER(TOLER = 1D-6)

      INTEGER ITER

      CALL REGFAL(0.4D0, 0.6D0, F, TOLER, X, ITER)
      PRINT *, "Question b: X is ", X

      PRINT *, "Question c: "
      CALL REGFAL(0.0D0, 1.4D0, G, TOLER, X, ITER)
      PRINT *, "X is ", X, " by false position method, in ", 
     &        ITER, " iterations"

      CALL BISECT(0.0D0, 1.4D0, G, TOLER, X, ITER)
      PRINT *, "X is ", X, " by bisection method, in ",
     &        ITER, " iterations"

      PRINT *, "False position needs many steps because g(x1) << g(x2)."
      PRINT *, "The correction in x is very small."
      END


      FUNCTION F(X)
      DOUBLE PRECISION F, X

      F = -2.0D0 + 6.2D0 * X - 4.0D0 * X**2 + 0.7D0 * X**3
      END

      FUNCTION G(X)
      DOUBLE PRECISION G, X

      G = X**10 - 0.95D0
      END


C     Regula falsi (false position) method
      SUBROUTINE REGFAL(A, B, FUNC, TOLER, X, ITER)
C     Initial limits
      DOUBLE PRECISION A, B
C     Tolerance
      DOUBLE PRECISION TOLER

C     Final approximation
      DOUBLE PRECISION X
      
      DOUBLE PRECISION FUNC
      EXTERNAL FUNC

C     Number of iterations needed
      INTEGER ITER

      DOUBLE PRECISION X1, X2, F1, F2, FX

      X1 = A
      X2 = B
      ITER = 0

      F1 = FUNC(X1)
      F2 = FUNC(X2)

      IF (F1 * F2 .GT. 0.0D0) THEN 
         PRINT *, "Limits are wrong"
         GOTO 100
      ENDIF

 200  X = X1 - F1 * (X2-X1) / (F2-F1)
      FX = FUNC(X)

C     Check if root is found
      IF (ABS(FX) .LT. TOLER) GOTO 100

C     New limits
      IF (F1 * FX > 0.D0) THEN
         X1 = X
         F1 = FX
      ELSE
         X2 = X
         F2 = FX
      END IF
      
      ITER = ITER + 1
      GOTO 200
      
 100  END



C     Bisection method
      SUBROUTINE BISECT(A, B, FUNC, TOLER, X, ITER)
C     Initial limits
      DOUBLE PRECISION A, B
C     Tolerance
      DOUBLE PRECISION TOLER

C     Final approximation
      DOUBLE PRECISION X

      DOUBLE PRECISION FUNC
      EXTERNAL FUNC

C     Number of iterations needed
      INTEGER ITER

      DOUBLE PRECISION X1, X2, F1, F2, FX

      X1 = A
      X2 = B
      ITER = 0

      F1 = FUNC(X1)
      F2 = FUNC(X2)

      IF (F1 * F2 .GT. 0.0D0) THEN 
         PRINT *, "Limits are wrong"
         GOTO 100
      ENDIF

 200  X = (X1 + X2) / 2.0D0
      FX = FUNC(X)

C     Check if root is found
      IF (ABS(FX) .LT. TOLER) GOTO 100

C     New limits
      IF (F1 * FX > 0.D0) THEN
         X1 = X
         F1 = FX
      ELSE
         X2 = X
         F2 = FX
      END IF

      ITER = ITER + 1
      GOTO 200
      
 100  END
