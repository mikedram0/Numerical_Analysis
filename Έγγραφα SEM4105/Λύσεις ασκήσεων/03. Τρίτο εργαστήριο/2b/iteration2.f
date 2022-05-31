C  Fixed point method to locate roots of 
C   
C  f(x)=x-cos^3 x  near 0.6.
C
C  Exercise 4b, chapter 2.

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
      X = 0.7D0

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

      F = X-COS(X)**3
      END


C     if f(x) above then 
C     f(x) = 0 => x = g(x) 
C     with g as follows:
C
      FUNCTION G(X)
      DOUBLE PRECISION X,G

      G = ACOS(X**(1.D0/3.D0))

C     g(x)=cos^3(x) is not appropriate.
      END
