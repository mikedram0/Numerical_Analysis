C     Compute epsilon for real and double precision numbers.
C
      PROGRAM EPSIL 

      REAL EPSRE
      DOUBLE PRECISION EPSDP

      EPSRE = 1.0

 10   EPSRE = 0.5 * EPSRE
      IF ((1.0 + EPSRE) .NE. 1.0) GOTO 10

C     here epsre is half the epsilon.
      PRINT *, "real epsilon is ", 2.0 * EPSRE
      
      EPSDP = 1.D0

 20   EPSDP = 0.5D0 * EPSDP
      IF ((1.0D0 + EPSDP) .NE. 1.D0) GOTO 20

C     here epsdp is half the epsilon.
      PRINT *, "double precision epsilon is ", 2.0D0 * EPSDP
      
      END
