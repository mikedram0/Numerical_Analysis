C  Bisection method to locate roots of 
C   
C  f(x)=sqrt(x)-cos(x) in [0,1].
C
C  Exercise 1b, chapter 2.
C
      PROGRAM BISECT
C     left and right limits bracketing the root 
      DOUBLE PRECISION A,B

C     approximation for the root
      DOUBLE PRECISION X

C     function value
      DOUBLE PRECISION FNVAL, FA, FB

C     function name type declarations
      DOUBLE PRECISION F
      LOGICAL SMSGN

C     declare that F, SMSGN are functions
      EXTERNAL F,SMSGN

C     a small constant number
      DOUBLE PRECISION TOLER
      PARAMETER (TOLER = 1D-8)

C     initial limits
      A=0.D0
      B=1.D0

      FA = F(A)
      FB = F(B)
      
      IF (SMSGN(FA,FB)) STOP
      
 100  X=(A+B)/2.D0

      FNVAL = F(X)

      PRINT *, "The current approximation is ", X 
      PRINT *, "The current value of F(x) is ", FNVAL
      PRINT *

C     Check if root is found
      IF (ABS(FNVAL) .LT. TOLER) STOP

      IF (SMSGN(FNVAL,FA)) THEN 
         B = X
         FB - FNVAL
      ELSE
         A = X
         FA = FNVAL
      ENDIF
      GOTO 100

      END


      FUNCTION F(X)
      DOUBLE PRECISION X,F

      F = SQRT(X)-COS(X)
      END


C     same sign
      FUNCTION SMSGN(A,B)
      DOUBLE PRECISION A,B
      LOGICAL SMSGN

C     SMSGN = A*B > 0.0D0
      SMSGN = SIGN(A,B) .EQ. A
      END
