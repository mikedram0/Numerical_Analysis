PROGRAM regula_falsi
  IMPLICIT NONE
  DOUBLE PRECISION :: X
  DOUBLE PRECISION, PARAMETER :: TOLER = 1D-6
  INTEGER :: ITER

  INTERFACE 
     FUNCTION F(X)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: X
       DOUBLE PRECISION :: F
     END FUNCTION F

     FUNCTION G(X)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: X
       DOUBLE PRECISION :: G
     END FUNCTION G

     SUBROUTINE BISECT(A, B, FUNC, TOLER, X, ITER)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in)  :: A, B   !     Initial limits
       DOUBLE PRECISION, INTENT (in)  :: TOLER  !     Tolerance
       DOUBLE PRECISION, INTENT (out) :: X      !     Final approximation
       INTEGER, INTENT (out)          :: ITER   ! Number of iterations needed
       INTERFACE
          FUNCTION func(x)
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT (in) :: x
            DOUBLE PRECISION :: FUNC
          END FUNCTION func
       END INTERFACE
     END SUBROUTINE BISECT

     SUBROUTINE REGFAL(A, B, FUNC, TOLER, X, ITER)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in)  :: A, B   !     Initial limits
       DOUBLE PRECISION, INTENT (in)  :: TOLER  !     Tolerance
       DOUBLE PRECISION, INTENT (out) :: X      !     Final approximation
       INTEGER, INTENT (out)          :: ITER   ! Number of iterations needed

       INTERFACE 
          FUNCTION func(x)
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT (in) :: x
            DOUBLE PRECISION :: FUNC
          END FUNCTION func
       END INTERFACE
     END SUBROUTINE REGFAL

  END INTERFACE

  CALL REGFAL(0.4D0, 0.6D0, F, TOLER, X, ITER)
  PRINT *, "Question b: X is ", X

  PRINT *, "Question c: "
  CALL REGFAL(0.0D0, 1.4D0, G, TOLER, X, ITER)
  PRINT *, "X is ", X, " by false position method, in ", ITER, " iterations"

  CALL BISECT(0.0D0, 1.4D0, G, TOLER, X, ITER)
  PRINT *, "X is ", X, " by bisection method, in ", ITER, " iterations"

  PRINT *, "False position needs many steps because g(x1) << g(x2)."
  PRINT *, "The correction in x is very small."
END PROGRAM regula_falsi


FUNCTION F(X)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: X
  DOUBLE PRECISION :: F

  F = -2.0D0 + 6.2D0 * X - 4.0D0 * X**2 + 0.7D0 * X**3
END FUNCTION F

FUNCTION G(X)
  DOUBLE PRECISION, INTENT (in) :: X
  DOUBLE PRECISION :: G

  G = X**10 - 0.95D0
END FUNCTION G


!     Regula falsi (false position) method
SUBROUTINE REGFAL(A, B, FUNC, TOLER, X, ITER)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in)  :: A, B   !     Initial limits
  DOUBLE PRECISION, INTENT (in)  :: TOLER  !     Tolerance
  DOUBLE PRECISION, INTENT (out) :: X      !     Final approximation
  INTEGER, INTENT (out)          :: ITER   !     Number of iterations needed

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
  
  IF (F1 * F2 >= 0.0D0) THEN 
     PRINT *, "Limits are wrong"
     RETURN
  ENDIF
  
  ITER = 0
  
  DO
     X = X1 - F1 * (X2-X1) / (F2-F1)
     FX = FUNC(X)

     !     Check IF root is found
     IF (ABS(FX) < TOLER) EXIT

     !     NEW limits
     IF (F1 * FX > 0.D0) THEN
        X1 = X
        F1 = FX
     ELSE
        X2 = X
        F2 = FX
     END IF

     ITER = ITER + 1
  ENDDO

END SUBROUTINE REGFAL



!     Bisection method
SUBROUTINE BISECT(A, B, FUNC, TOLER, X, ITER)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in)  :: A, B   !     Initial limits
  DOUBLE PRECISION, INTENT (in)  :: TOLER  !     Tolerance
  DOUBLE PRECISION, INTENT (out) :: X      !     Final approximation
  INTEGER, INTENT (out)          :: ITER   !     Number of iterations needed

  INTERFACE
     FUNCTION func(x)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: x
       DOUBLE PRECISION :: func
     END FUNCTION func
  END INTERFACE


  DOUBLE PRECISION :: X1, X2, F1, F2, FX

  X1 = A
  X2 = B

  F1 = FUNC(X1)
  F2 = FUNC(X2)
  
  IF (F1 * F2 > 0.0D0) THEN 
     PRINT *, "Limits are wrong"
     RETURN
  ENDIF
  
  ITER = 0
  
  DO
     X = (X1 + X2) / 2.0D0     
     FX = FUNC(X)

     !     Check if root is found
     IF (ABS(FX) < TOLER) EXIT

     !     NEW limits
     IF (F1 * FX > 0.D0) THEN
        X1 = X
        F1 = FX
     ELSE
        X2 = X
        F2 = FX
     END IF

     ITER = ITER + 1
  ENDDO
END SUBROUTINE BISECT
