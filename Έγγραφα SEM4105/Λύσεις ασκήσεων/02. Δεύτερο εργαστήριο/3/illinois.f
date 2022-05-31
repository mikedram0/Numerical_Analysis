      PROGRAM METHOD
      DOUBLE PRECISION  X
      DOUBLE PRECISION TOLER
      PARAMETER (TOLER = 1D-9)
      INTEGER ITER, EVAL

      DOUBLE PRECISION F
      EXTERNAL F, BISECT, REGFAL, ILLINO
      
      CALL BISECT(0.0D0, 1.0D0, F, TOLER, X, ITER, EVAL)
      PRINT *, "X IS ", X, " BY BISECTION METHOD, IN ",
     &        ITER, " ITERATIONS AND ", EVAL, " FUNCTION EVALUATIONS"
      
      CALL REGFAL(0.0D0, 1.0D0, F, TOLER, X, ITER, EVAL)
      PRINT *, "X IS ", X, " BY FALSE POSITION METHOD, IN ",
     &        ITER, " ITERATIONS AND ", EVAL, " FUNCTION EVALUATIONS"
      
      CALL ILLINO(0.0D0, 1.0D0, F, TOLER, X, ITER, EVAL)
      PRINT *, "X IS ", X, " BY ILLINOIS METHOD, IN ",
     &        ITER, " ITERATIONS AND ", EVAL, " FUNCTION EVALUATIONS"

      END


      FUNCTION F(X)
      DOUBLE PRECISION X
      DOUBLE PRECISION F
      
      F = X*X -(1.0D0-X)**5
      END 


C     REGULA FALSI (FALSE POSITION) METHOD
      SUBROUTINE REGFAL(A, B, FUNC, TOLER, X, ITER, EVAL)
C     INITIAL LIMITS
      DOUBLE PRECISION A, B

C     TOLERANCE
      DOUBLE PRECISION TOLER

C     FINAL APPROXIMATION
      DOUBLE PRECISION X      

C     NUMBER OF ITERATIONS NEEDED
      INTEGER ITER              

C     NUMBER OF FUNCTION EVALUATIONS
      INTEGER EVAL

      DOUBLE PRECISION FUNC
      EXTERNAL FUNC
      
      DOUBLE PRECISION X1, X2, F1, F2, FX

      X1 = A
      X2 = B
      
      F1 = FUNC(X1)
      F2 = FUNC(X2)

      EVAL = 2
      
      IF (F1 * F2 .GE. 0.0D0) THEN 
         PRINT *, "LIMITS ARE WRONG"
         RETURN
      ENDIF

      ITER = 0
      
 10   X = X1 - F1 * (X2-X1) / (F2-F1)
      FX = FUNC(X)
      EVAL = EVAL + 1

      ITER = ITER + 1
           
C     CHECK IF ROOT IS FOUND
      IF (ABS(FX) .LT. TOLER) RETURN 
      
C     NEW LIMITS
      IF (F1 * FX .GT. 0.0D0) THEN
         X1 = X
         F1 = FX
      ELSE
         X2 = X
         F2 = FX
      END IF
      
      GOTO 10
      
      END 



C     BISECTION METHOD
      SUBROUTINE BISECT(A, B, FUNC, TOLER, X, ITER, EVAL)
C     INITIAL LIMITS
      DOUBLE PRECISION A, B

C     TOLERANCE
      DOUBLE PRECISION TOLER

C     FINAL APPROXIMATION
      DOUBLE PRECISION X      

C     NUMBER OF ITERATIONS NEEDED
      INTEGER ITER

C     NUMBER OF FUNCTION EVALUATIONS
      INTEGER EVAL

      DOUBLE PRECISION FUNC
      EXTERNAL FUNC
      
      DOUBLE PRECISION X1, X2, F1, F2, FX

      X1 = A
      X2 = B
      
      F1 = FUNC(X1)
      F2 = FUNC(X2)

      EVAL = 2
      
      IF (F1 * F2 .GE. 0.0D0) THEN 
         PRINT *, "LIMITS ARE WRONG"
         RETURN
      ENDIF

      ITER = 0
      
 10   X = (X1 + X2) / 2.0D0     
      FX = FUNC(X)
      EVAL = EVAL + 1
      
      ITER = ITER + 1

C     CHECK IF ROOT IS FOUND
      IF (ABS(FX) .LT. TOLER) RETURN

C     NEW LIMITS
      IF (F1 * FX .GT. 0.0D0) THEN
         X1 = X
         F1 = FX
      ELSE
         X2 = X
         F2 = FX
      END IF
      
      GOTO 10
      END


      
C     REGULA FALSI (FALSE POSITION) METHOD; ILLINOIS ALGORITHM.
      SUBROUTINE ILLINO(A, B, FUNC, TOLER, X, ITER, EVAL)
C     INITIAL LIMITS
      DOUBLE PRECISION A, B

C     TOLERANCE
      DOUBLE PRECISION TOLER

C     FINAL APPROXIMATION
      DOUBLE PRECISION X      

C     NUMBER OF ITERATIONS NEEDED
      INTEGER ITER              

C     NUMBER OF FUNCTION EVALUATIONS
      INTEGER EVAL

      DOUBLE PRECISION FUNC
      EXTERNAL FUNC
      
      DOUBLE PRECISION X1, X2, F1, F2, FX
      INTEGER CHNGA, CHNGB
  
      X1 = A
      X2 = B
      F1 = FUNC(X1)
      F2 = FUNC(X2)

      EVAL = 2
      
      IF (F1 * F2 .GE. 0.0D0) THEN 
         PRINT *, "LIMITS ARE WRONG"
         RETURN
      ENDIF

      ITER = 0
      CHNGA = 0
      CHNGB = 0
      
 10   IF (CHNGA .LT. 2 .AND. CHNGB .LT. 2) THEN
         X = (X1 * F2 - X2 * F1) / (F2-F1)
      END IF
      
      IF (CHNGA .EQ. 2) THEN
         X = (X1 * F2 - 2D0 * X2 * F1) / (F2-2D0*F1)
         CHNGA = 0
      END IF

      IF (CHNGB .EQ. 2) THEN
         X = (2D0 * X1 * F2 - X2 * F1) / (2D0*F2-F1)
         CHNGB = 0
      END IF
              
      FX = FUNC(X)
      EVAL = EVAL + 1
      
      ITER = ITER + 1

C     CHECK IF ROOT IS FOUND
      IF (ABS(FX) .LT. TOLER) RETURN
      
      IF (F1 * FX .GT. 0.0D0) THEN
         X1 = X
         F1 = FX
         
         CHNGB = 0
         CHNGA = CHNGA + 1
      ELSE
         X2 = X
         F2 = FX
         
         CHNGA = 0
         CHNGB = CHNGB + 1
      ENDIF

      GOTO 10

      END
