C     FORTRAN PROGRAM TO COMPUTE EXP(X), USING THE SERIES
C         1 + X/1! + X^2/2! + X^3/3! + ...
C 
C     Exercise 4 of chapter 1.
C
      PROGRAM EXPSER

      DOUBLE PRECISION SUM, X, A
      INTEGER N

      PRINT *, "Give x: "
      READ *, X

      SUM = 0.D0

C     initial value of A is the term for n=0.
      N = 0
      A = 1.D0 

 10   SUM = SUM + A
C     next term is ...
      N = N+1
      A = A * X / N
      IF (SUM + A .NE. SUM) GOTO 10

      PRINT *, "The sum of Taylor series for EXP at ", X, " is ", SUM
      PRINT *, "The correct value is ", EXP(X)

      END
