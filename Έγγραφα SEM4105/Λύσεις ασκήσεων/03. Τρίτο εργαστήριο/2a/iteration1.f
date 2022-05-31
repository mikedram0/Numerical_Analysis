C  Fixed point method to locate a root of 
C   
C   f(x)=x^2-6x+5
C
     PROGRAM FIXPOI

      DOUBLE PRECISION X

C     declare that F,G are functions
      EXTERNAL F,G

C     function names type declaration
      DOUBLE PRECISION F,G

C     function value
      DOUBLE PRECISION FNVAL

C     a small constant number
      DOUBLE PRECISION TOLER
      PARAMETER (TOLER = 1D-8)

C     a very large number
      DOUBLE PRECISION VFAR
      PARAMETER (VFAR = 1D5)

C     initial guess for x
      X = 1.3D0

 100  FNVAL = F(X)

      PRINT *, "The current approximation is ", X 
      PRINT *, "The current value of F(x) is ", FNVAL
      PRINT *

      IF (ABS(X) .GT. VFAR) THEN
         PRINT *, "X is very large; probably g(x) is not appropriate"
         GOTO 200
      ENDIF

      IF (ABS(FNVAL) .LT. TOLER) GOTO 200

      X = G(X)
      GOTO 100

 200  END


      FUNCTION F(X)
      DOUBLE PRECISION X,F

      F = X**2 - 6.0D0 * X + 5.0D0
      END


C     if f(x) above then 
C     f(x) = 0 => x = g(x) 
C     with g as follows:

      FUNCTION G(X)
      DOUBLE PRECISION X,G

      G = (X**2 + 5.0D0) / 6.0D0
      END
