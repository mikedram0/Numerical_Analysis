C   Rizes tis f(x)=4 cos(x) - exp(-x) me tis methodous 
C   dihotomisis, statheroy shmeioy, Newton-Raphson, temnousas.

      PROGRAM ALLMET

C     A small constant number
      DOUBLE PRECISION TOLER
      PARAMETER (TOLER = 1D-8)

      DOUBLE PRECISION X, A, B

C     declare that F, DF, G are functions
      EXTERNAL F, DF, G

C     declare the subroutines
      EXTERNAL NEWTRA, BISECT, FIXED, SECANT

C /////    NEWTON-RAPHSON    ///////

C     Initial approximation 
      X = 1.0D0

       CALL NEWTRA(X, TOLER, F, DF)

       CALL OUTPUT("NEWTON", X, F)

C /////    BISECTION    ///////

C     Initial boundaries
      A = 1.0D0
      B = 2.0D0

      CALL BISECT(A, B, X, TOLER, F)

      CALL OUTPUT("BISECTION", X, F)

C /////   FIXED POINT   ///////

C     Initial approximation 
      X = 1.0D0

      CALL FIXED(X, TOLER, F, G)

      CALL OUTPUT("FIXED", X, F)

C /////   SECANT   ///////

C     Initial approximations
      A = -1.0D0
      B = 1.0D0

      CALL SECANT(A, B, X, TOLER, F)

      CALL OUTPUT("SECANT", X, F)

      END

C******************************************

      SUBROUTINE OUTPUT(METHOD, X, FUNC)
      CHARACTER(*) METHOD
      DOUBLE PRECISION X, FUNC
      EXTERNAL FUNC

      PRINT *, "The method is ", METHOD
      PRINT *, "A root is approximately ", X
      PRINT *, "The function value is ", FUNC(X)

      END

C     The function
      FUNCTION F(X)
      DOUBLE PRECISION X, F

      F = 4.0D0 * COS(X) - EXP(-X)
      END

C     The first derivative 
      FUNCTION DF(X)
      DOUBLE PRECISION X, DF

      DF = -4.0D0 * SIN(X) + EXP(-X)
      END


      FUNCTION G(X)
      DOUBLE PRECISION X, G

      G = ACOS(EXP(-X) / 4.0D0)
      END

C******************************************

      SUBROUTINE NEWTRA(X, TOLER, FUNC, DERIV)

      DOUBLE PRECISION X
C     Input  : Initial approximation 
C     Output : Final   approximation

C     accuracy
      DOUBLE PRECISION TOLER

C     Function and derivative
      DOUBLE PRECISION FUNC, DERIV
      EXTERNAL FUNC, DERIV

 10   X = X - FUNC(X) / DERIV(X)

      IF (ABS(FUNC(X)) .GT. TOLER) GOTO 10

      END

C******************************************

      SUBROUTINE BISECT(A, B, X, TOLER, FUNC)

C     Initial boundaries
      DOUBLE PRECISION A, B

      DOUBLE PRECISION X
C     Output : Final approximation

C     accuracy
      DOUBLE PRECISION TOLER

C     Function
      DOUBLE PRECISION FUNC
      EXTERNAL FUNC

 100  X=(A+B)/2.D0

C     Check if root is found
      IF (ABS(FUNC(X)) .LT. TOLER) RETURN

      IF (FUNC(A)*FUNC(X) .LT. 0.D0) THEN 
         B = X
      ELSE
         A = X
      ENDIF
      GOTO 100

      END

C******************************************

      SUBROUTINE FIXED(X, TOLER, F, G)

      DOUBLE PRECISION X, TOLER
      DOUBLE PRECISION F, G
      EXTERNAL F, G

C     a very large number
      DOUBLE PRECISION VFAR
      PARAMETER (VFAR = 1D5)

 100  IF (ABS(X) .GT. VFAR) THEN
         PRINT *, "X is very large; probably g(x) is not appropriate"
         RETURN
      ENDIF

      IF (ABS(F(X)) .LT. TOLER) RETURN

      X = G(X)
      GOTO 100
      
      END

C******************************************

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
