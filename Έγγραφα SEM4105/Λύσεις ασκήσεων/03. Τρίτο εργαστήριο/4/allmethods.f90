!   Rizes tis f(x)=4 cos(x) - exp(-x) me tis methodous 
!   dihotomisis, statheroy shmeioy, Newton-Raphson, temnousas.

PROGRAM all_methods
  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: TOLER = 1D-8 ! A small constant number
  DOUBLE PRECISION :: X, A, B
  
  INTERFACE 
     FUNCTION F(X)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: X
       DOUBLE PRECISION :: F
     END FUNCTION F

     FUNCTION DF(X)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: X
       DOUBLE PRECISION :: DF
     END FUNCTION DF

     FUNCTION G(X)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: X
       DOUBLE PRECISION :: G
     END FUNCTION G

     SUBROUTINE NEWTRA(X, TOLER, FUNC, DERIV)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (inout) :: X
       DOUBLE PRECISION, INTENT (in) :: TOLER 
       INTERFACE 
          FUNCTION FUNC(X)
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT (in) :: X
            DOUBLE PRECISION :: FUNC
          END FUNCTION FUNC

          FUNCTION DERIV(X)
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT (in) :: X
            DOUBLE PRECISION :: DERIV
          END FUNCTION DERIV
       END INTERFACE
     END SUBROUTINE NEWTRA

     SUBROUTINE BISECT(A, B, X, TOLER, FUNC)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (inout) :: A, B 
       DOUBLE PRECISION, INTENT (out)   :: X     
       DOUBLE PRECISION, INTENT (in) :: TOLER 

       INTERFACE 
          FUNCTION FUNC(X)
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT (in) :: X
            DOUBLE PRECISION :: FUNC
          END FUNCTION FUNC
       END INTERFACE
     end SUBROUTINE BISECT

     SUBROUTINE FIXED(X, TOLER, F, G)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (inout) :: X
       DOUBLE PRECISION, INTENT (in)   :: TOLER

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
       END INTERFACE
     END SUBROUTINE FIXED

     SUBROUTINE SECANT(X1, X2, X, TOLER, FUNC)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (inout) :: X1, X2
       DOUBLE PRECISION, INTENT (out)   :: X
       DOUBLE PRECISION, INTENT (in)    :: TOLER
       
       INTERFACE
          FUNCTION FUNC(X)
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT (in) :: X
            DOUBLE PRECISION :: FUNC
          END FUNCTION FUNC
       END INTERFACE
     END SUBROUTINE SECANT
  END INTERFACE

  ! /////    NEWTON-RAPHSON    ///////

  X = 1.0D0    !     Initial approximation 

  CALL NEWTRA(X, TOLER, F, DF)

  CALL OUTPUT("NEWTON", X, F)

  ! /////    BISECTION    ///////

  !     Initial boundaries
  A = 1.0D0
  B = 2.0D0

  CALL BISECT(A, B, X, TOLER, F)

  CALL OUTPUT("BISECTION", X, F)

  ! /////   FIXED POINT   ///////

  X = 1.0D0 !     Initial approximation 

  CALL FIXED(X, TOLER, F, G)

  CALL OUTPUT("FIXED", X, F)

  ! /////   SECANT   ///////

  !     Initial approximations
  A = 1.0D0
  B = 2.0D0

  CALL SECANT(A, B, X, TOLER, F)

  CALL OUTPUT("SECANT", X, F)

END PROGRAM all_methods

!******************************************

SUBROUTINE OUTPUT(METHOD, X, FUNC)
  IMPLICIT NONE
  CHARACTER(*), INTENT (in) :: METHOD
  DOUBLE PRECISION, INTENT (in) ::  X

  INTERFACE 
     FUNCTION FUNC(X)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: X
       DOUBLE PRECISION :: FUNC
     END FUNCTION FUNC
  END INTERFACE

  PRINT *, "The method is ", METHOD
  PRINT *, "A root is approximately ", X
  PRINT *, "The function value is ", FUNC(X)
END SUBROUTINE OUTPUT

!     The function
FUNCTION F(X)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: X
  DOUBLE PRECISION :: F

  F = 4.0D0 * COS(X) - EXP(-X)
END FUNCTION F

!     The first derivative 
FUNCTION DF(X)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: X
  DOUBLE PRECISION :: DF

  DF = -4.0D0 * SIN(X) + EXP(-X)
END FUNCTION DF


FUNCTION G(X)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: X
  DOUBLE PRECISION :: G

  G = ACOS(EXP(-X) / 4.0D0)
END FUNCTION G

!******************************************

SUBROUTINE NEWTRA(X, TOLER, FUNC, DERIV)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (inout) :: X
  !     Input  : Initial approximation 
  !     Output : Final   approximation
  DOUBLE PRECISION, INTENT (in) :: TOLER   !     accuracy

  !     Function and derivative
  INTERFACE 
     FUNCTION FUNC(X)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: X
       DOUBLE PRECISION :: FUNC
     END FUNCTION FUNC

     FUNCTION DERIV(X)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: X
       DOUBLE PRECISION :: DERIV
     END FUNCTION DERIV
  END INTERFACE

  DO 
     X = X - FUNC(X) / DERIV(X)
     IF (ABS(FUNC(X)) < TOLER) EXIT
  ENDDO

END SUBROUTINE NEWTRA

!******************************************

SUBROUTINE BISECT(A, B, X, TOLER, FUNC)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (inout) :: A, B  ! Initial boundaries
  DOUBLE PRECISION, INTENT (out)   :: X     ! Output : Final approximation
  DOUBLE PRECISION, INTENT (in) :: TOLER    ! accuracy

  INTERFACE 
     FUNCTION FUNC(X)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: X
       DOUBLE PRECISION :: FUNC
     END FUNCTION FUNC
  END INTERFACE

  DO
     X = (A+B)/2.D0

     !     Check if root is found
     IF (ABS(FUNC(X)) < TOLER) EXIT

     IF (FUNC(A)*FUNC(X) < 0.D0) THEN 
        B = X
     ELSE
        A = X
     ENDIF
  ENDDO

END SUBROUTINE BISECT

!******************************************

SUBROUTINE FIXED(X, TOLER, F, G)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (inout) :: X
  DOUBLE PRECISION, INTENT (in)   :: TOLER

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
  END INTERFACE

  DOUBLE PRECISION, PARAMETER :: VFAR = 1D5 ! a very large number

  DO
     IF (ABS(X) > VFAR) THEN
        PRINT *, "X is very large; probably g(x) is not appropriate"
        EXIT
     ENDIF

     IF (ABS(F(X)) < TOLER) EXIT

     X = G(X)
  ENDDO

END SUBROUTINE FIXED

!******************************************

SUBROUTINE SECANT(X1, X2, X, TOLER, FUNC)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (inout) :: X1, X2
  DOUBLE PRECISION, INTENT (out)   :: X
  DOUBLE PRECISION, INTENT (in)    :: TOLER
  
  INTERFACE
     FUNCTION FUNC(X)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: X
       DOUBLE PRECISION :: FUNC
     END FUNCTION FUNC
  END INTERFACE

  DOUBLE PRECISION :: F1, F2

  F1 = FUNC(X1)
  F2 = FUNC(X2)
  
  DO 
     X = X2 - F2 * (X2 - X1) / (F2 - F1)

     X1 = X2
     F1 = F2

     X2 = X
     F2 = FUNC(X)

     IF (ABS(F2) < TOLER) EXIT
  END DO

END SUBROUTINE SECANT
