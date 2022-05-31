      PROGRAM TAYLOR

      DOUBLE PRECISION A, B, H
      PARAMETER (A=-1.0D0, B=1.0D0, H=0.01D0)
      DOUBLE PRECISION YINIT
      PARAMETER (YINIT = 3.0D0)

      DOUBLE PRECISION XOLD, YOLD, XNEW, YNEW, DY(5)

      XOLD = A
      YOLD = YINIT

 10   XNEW = XOLD + H

      CALL DIFEQ(XOLD, YOLD, DY)

      YNEW = YOLD + (DY(1) + (DY(2) + (DY(3)
     &  + (DY(4) + DY(5) * H / 5) * H/4) * H/3) * H/2) * H
      
      PRINT *, XNEW, YNEW

      XOLD = XNEW
      YOLD = YNEW

      IF (XOLD .LT. B) GOTO 10

      END

      SUBROUTINE DIFEQ(X,Y,DY)
      DOUBLE PRECISION X,Y, DY(5)

      DY(1) = COS(X)-SIN(Y) + X*X

      DY(2) = 2.0D0 * X - SIN(X) - COS(Y) * DY(1)

      DY(3) = 2.0D0 - COS(X) - COS(Y) * DY(2) + SIN(Y) * DY(1)**2 

      DY(4) = SIN(X) + COS(Y) * (DY(1)**3 - DY(3))
     &        + 3 * DY(1) * DY(2) * SIN(Y)

      DY(5) = COS(X) + COS(Y) * (6 * DY(1)**2 * DY(2) - DY(4)) 
     &        + (3 * DY(2)**2 + DY(1) * (4 * DY(3) - DY(1)**3)) * SIN(Y)

      END

