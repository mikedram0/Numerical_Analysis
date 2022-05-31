      SUBROUTINE DERIV(X,Y,Z, DY, DZ)
      DOUBLE PRECISION X, Y, Z, DY(4), DZ(4)

      DY(1) = Y + Z*Z - X**3
      DZ(1) = Y**3 + Z + COS(X)

      DY(2) =-3.0D0 * X*X + DY(1) + 2.0D0 * Z * DZ(1)
      DZ(2) = 3.0D0 * Y*Y * DY(1) + DZ(1) - SIN(X)
 
      DY(3) = -6.0D0 * X + DY(2) + 2.0D0 * (DZ(1)**2 + Z * DZ(2))
      DZ(3) = 3.0D0 * Y * (2.0D0 * DY(1)**2 + Y * DY(2)) 
     &     + DZ(2) - COS(X)

      DY(4) =-6.0D0 + DY(3) + 2.0D0 * (3.0D0*DZ(1) * DZ(2)+ Z*DZ(3)) 
      DZ(4) = DZ(3) + SIN(X) + 6.0D0 * DY(1)**3 
     &     + 3.0D0 * Y * (6.0D0 * DY(1) * DY(2) + Y * DY(3))
      END


      PROGRAM TAYLOR

      DOUBLE PRECISION A, B, H
      PARAMETER (A=0.0D0, B=1.0D0, H=0.1D0)
      DOUBLE PRECISION YINIT, ZINIT
      PARAMETER (YINIT = 0.3D0, ZINIT = 0.1D0)

      DOUBLE PRECISION XOLD, YOLD, ZOLD, XNEW, YNEW, ZNEW
      DOUBLE PRECISION DY(4), DZ(4)
      INTEGER K, NSTEPS
      
      XOLD = A
      YOLD = YINIT
      ZOLD = ZINIT

      NSTEPS = NINT((B-A)/H)

      DO K=1,NSTEPS
        XNEW = XOLD + H
        
        CALL DERIV(XOLD, YOLD, ZOLD, DY, DZ)
        
        YNEW = YOLD + H * (DY(1) + H / 2D0 * (DY(2) + H / 3D0 * (DY(3)
     &       + H / 4D0 *  DY(4))))
        
        ZNEW = ZOLD + H * (DZ(1) + H / 2D0 * (DZ(2) + H / 3D0 * (DZ(3)
     &       + H / 4D0 *  DZ(4))))
        
        PRINT *, XNEW, YNEW, ZNEW
        
        XOLD = XNEW
        YOLD = YNEW
        ZOLD = ZNEW
      ENDDO
      
      END


      
