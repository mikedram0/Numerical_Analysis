PROGRAM methods
  IMPLICIT NONE
  DOUBLE PRECISION :: X
  DOUBLE PRECISION, PARAMETER :: TOLER = 1D-9
  INTEGER :: ITER, EVAL

  INTERFACE 
     FUNCTION F(X)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: X
       DOUBLE PRECISION :: F
     END FUNCTION F

     SUBROUTINE BISECT(A, B, FUNC, TOLER, X, ITER, EVAL)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in)  :: A, B   !     Initial limits
       DOUBLE PRECISION, INTENT (in)  :: TOLER  !     Tolerance
       DOUBLE PRECISION, INTENT (out) :: X      !     Final approximation
       INTEGER, INTENT (out)          :: ITER   ! Number of iterations needed
       INTEGER, INTENT (out)          :: EVAL   !     Number of function evaluations
       INTERFACE
          FUNCTION func(x)
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT (in) :: x
            DOUBLE PRECISION :: FUNC
          END FUNCTION func
       END INTERFACE
     END SUBROUTINE BISECT

     SUBROUTINE REGFAL(A, B, FUNC, TOLER, X, ITER, EVAL)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in)  :: A, B   !     Initial limits
       DOUBLE PRECISION, INTENT (in)  :: TOLER  !     Tolerance
       DOUBLE PRECISION, INTENT (out) :: X      !     Final approximation
       INTEGER, INTENT (out)          :: ITER   ! Number of iterations needed
       INTEGER, INTENT (out)          :: EVAL   !     Number of function evaluations
       INTERFACE 
          FUNCTION func(x)
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT (in) :: x
            DOUBLE PRECISION :: FUNC
          END FUNCTION func
       END INTERFACE
     END SUBROUTINE REGFAL

     SUBROUTINE ILLINO(A, B, FUNC, TOLER, X, ITER, EVAL)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in)  :: A, B   !     Initial limits
       DOUBLE PRECISION, INTENT (in)  :: TOLER  !     Tolerance
       DOUBLE PRECISION, INTENT (out) :: X      !     Final approximation
       INTEGER, INTENT (out)          :: ITER   !     Number of iterations needed
       INTEGER, INTENT (out)          :: EVAL   !     Number of function evaluations
       INTERFACE 
          FUNCTION func(x)
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT (in) :: X
            DOUBLE PRECISION :: FUNC
          END FUNCTION func
       END INTERFACE
     END SUBROUTINE ILLINO
  END INTERFACE


  CALL BISECT(0.0D0, 1.D0, F, TOLER, X, ITER, EVAL)
  PRINT *, "X is ", X, " by bisection method, in ", ITER, &
       " iterations and ",  EVAL, " function evaluations"

  CALL REGFAL(0.0D0, 1.D0, F, TOLER, X, ITER, EVAL)
  PRINT *, "X is ", X, " by false position method, in ", ITER, &
       " iterations and ", EVAL, " function evaluations"

  CALL ILLINO(0.0D0, 1.D0, F, TOLER, X, ITER, EVAL)
  PRINT *, "X is ", X, " by illinois method, in ", ITER, &
       " iterations and ",  EVAL, " function evaluations"
END PROGRAM methods


FUNCTION F(X)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: X
  DOUBLE PRECISION :: F

  F = X*X -(1.0D0-X)**5
END FUNCTION F


!     Regula falsi (false position) method
SUBROUTINE REGFAL(A, B, FUNC, TOLER, X, ITER, EVAL)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in)  :: A, B   !     Initial limits
  DOUBLE PRECISION, INTENT (in)  :: TOLER  !     Tolerance
  DOUBLE PRECISION, INTENT (out) :: X      !     Final approximation
  INTEGER, INTENT (out)          :: ITER   !     Number of iterations needed
  INTEGER, INTENT (out)          :: EVAL   !     Number of function evaluations

  INTERFACE 
     FUNCTION func(x)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: X
       DOUBLE PRECISION :: FUNC
     END FUNCTION func
  END INTERFACE

  DOUBLE PRECISION :: X1, X2, F1, F2, FX

  X1 = A
  X2 = B

  F1 = FUNC(X1)
  F2 = FUNC(X2)

  EVAL=2
  
  IF (F1 * F2 >= 0.0D0) THEN 
     PRINT *, "Limits are wrong"
     RETURN
  ENDIF

  ITER = 0

  DO
     X = X1 - F1 * (X2-X1) / (F2-F1)
     FX = FUNC(X)
     EVAL = EVAL + 1
     
     ITER = ITER + 1

     ! Check if root is found
     IF (ABS(FX) < TOLER) EXIT

     !     new limits
     IF (F1 * FX > 0.D0) THEN
        X1 = X
        F1 = FX
     ELSE
        X2 = X
        F2 = FX
     END IF

  ENDDO

END SUBROUTINE REGFAL



!     Bisection method
SUBROUTINE BISECT(A, B, FUNC, TOLER, X, ITER, EVAL)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in)  :: A, B   !     Initial limits
  DOUBLE PRECISION, INTENT (in)  :: TOLER  !     Tolerance
  DOUBLE PRECISION, INTENT (out) :: X      !     Final approximation
  INTEGER, INTENT (out)          :: ITER   !     Number of iterations needed
  INTEGER, INTENT (out)          :: EVAL   !     Number of function evaluations
  
  INTERFACE
     FUNCTION func(x)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: x
       DOUBLE PRECISION :: FUNC
     END FUNCTION func
  END INTERFACE


  DOUBLE PRECISION :: X1, X2, F1, F2, FX

  X1 = A
  X2 = B

  F1 = FUNC(X1)
  F2 = FUNC(X2)

  EVAL = 2
  
  IF (F1 * F2 > 0.0D0) THEN 
     PRINT *, "Limits are wrong"
     RETURN
  ENDIF

  ITER = 0

  DO
     X = (X1 + X2) / 2.0D0     
     FX = FUNC(X)
     EVAL = EVAL + 1

     ITER = ITER + 1
     
     !     Check if root is found
     IF (ABS(FX) < TOLER) EXIT

     !     new limits
     IF (F1 * FX > 0.D0) THEN
        X1 = X
        F1 = FX
     ELSE
        X2 = X
        F2 = FX
     END IF
  ENDDO
END SUBROUTINE BISECT

!     Regula falsi (false position) method; illinois algorithm.
SUBROUTINE ILLINO(A, B, FUNC, TOLER, X, ITER, EVAL)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in)  :: A, B   !     Initial limits
  DOUBLE PRECISION, INTENT (in)  :: TOLER  !     Tolerance
  DOUBLE PRECISION, INTENT (out) :: X      !     Final approximation
  INTEGER, INTENT (out)          :: ITER   !     Number of iterations needed
  INTEGER, INTENT (out)          :: EVAL   !     Number of function evaluations
  
  INTERFACE 
     FUNCTION func(x)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: X
       DOUBLE PRECISION :: FUNC
     END FUNCTION func
  END INTERFACE

  DOUBLE PRECISION :: X1, X2, F1, F2, FX
  INTEGER :: chnga, chngb
  
  X1 = A
  X2 = B
  F1 = FUNC(X1)
  F2 = FUNC(X2)

  EVAL = 2
  
  IF (F1 * F2 > 0.0D0) THEN 
     PRINT *, "Limits are wrong"
     RETURN
  ENDIF

  ITER = 0
  chnga = 0
  chngb = 0

  DO
     IF (chnga < 2 .AND. chngb < 2) THEN
        X = (X1 * F2 - X2 * F1) / (F2-F1)
     END IF

     IF (chnga == 2) THEN
        X = (X1 * F2 - 2 * X2 * F1) / (F2-2*F1)
        chnga = 0
     END IF

     IF (chngb == 2) THEN
        X = (2 * X1 * F2 - X2 * F1) / (2*F2-F1)
        chngb = 0
     END IF
     

     FX = FUNC(X)
     EVAL = EVAL + 1

     ITER = ITER + 1
     
     !     Check if root is found
     IF (ABS(FX) < TOLER) EXIT

     IF (F1 * FX > 0.D0) THEN
        X1 = X
        F1 = FX

        chngb = 0
        chnga = chnga + 1
     ELSE
        X2 = X
        F2 = FX

        chnga = 0
        chngb = chngb + 1
     ENDIF


  ENDDO

END SUBROUTINE ILLINO
