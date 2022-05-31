!     FORTRAN PROGRAM TO COMPUTE EXP(X), USING THE SERIES
!         1 + X/1! + X^2/2! + X^3/3! + ...
! 
!     Exercise 4 of chapter 1.
!
PROGRAM expseries
  IMPLICIT NONE
  DOUBLE PRECISION :: sum, x, a
  INTEGER :: n
  
  PRINT *, "give x: "
  READ *, x
  
  sum = 0.d0   
  n = 0
  a = 1.d0 ! initial value of a is the term for n=0.
  DO 
     sum = sum + a
     !     next term is ...
     n = n+1
     a = a * x / n
     IF (sum + a == sum) EXIT
  ENDDO

  PRINT *, "The sum of Taylor series for EXP at ", x, " is ", sum
  PRINT *, "The correct value is ", EXP(x)
  
END PROGRAM expseries
