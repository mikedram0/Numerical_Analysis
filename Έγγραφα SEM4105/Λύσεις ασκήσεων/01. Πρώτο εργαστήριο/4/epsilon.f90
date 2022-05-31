!     Compute epsilon for real and double precision numbers.
!
PROGRAM epsil 
  IMPLICIT NONE
  REAL :: epsre
  DOUBLE PRECISION :: epsdp
  
  epsre = 1.0
  
  DO 
     epsre = 0.5 * epsre
     IF (1.0 + epsre == 1.0) EXIT
  ENDDO
  
  !     here epsre is half the epsilon.
  PRINT *, "real epsilon is ", 2.0 * epsre
  PRINT *, "from epsilon()", EPSILON(1.0)
  
  epsdp = 1.d0
  
  DO
     epsdp = 0.5d0 * epsdp
     IF (1.0d0 + epsdp == 1.d0) EXIT
  ENDDO
  !     here epsdp is half the epsilon.
  PRINT *, "double precision epsilon is ", 2.0d0 * epsdp
  PRINT *, "from epsilon()", EPSILON(1.0d0)

END PROGRAM epsil
