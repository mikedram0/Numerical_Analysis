      PROGRAM GAUSIN
      DOUBLE PRECISION GAUSS2, GAUSS3, FINT
      EXTERNAL GAUSS2, GAUSS3, F, FINT
      DOUBLE PRECISION A, B
      A = 2.1D0
      B = 5.2D0
      
      PRINT *, "Gauss 2 points:", GAUSS2(A, B, F)
      PRINT *, "Gauss 3 points:", GAUSS3(A, B, F)
      PRINT *, "correct: ", FINT(B) - FINT(A)
      END

      FUNCTION F(X)
      DOUBLE PRECISION X, F
      F = X**3 * EXP(-X)
      END

C     D FINT(X) / D X = F
      FUNCTION FINT(X)
      DOUBLE PRECISION X, FINT
      FINT = -(6.0D0 + 6.0D0 * X + 3.0D0 * X**2 + X**3) * EXP(-X)
      END


      FUNCTION GAUSS2(A, B, F)
      DOUBLE PRECISION A, B, GAUSS2, F
      EXTERNAL F

      DOUBLE PRECISION X(2), C(2)
      DOUBLE PRECISION SUM, Y
      INTEGER I
 
      C(1) = 1.0D0
      C(2) = 1.0D0

      X(1) = -1.0D0 / SQRT(3.0D0)
      X(2) =  1.0D0 / SQRT(3.0D0)

      SUM = 0.0D0
      DO I=1,2
        Y = ( (B-A) * X(I) + (B+A) ) / 2.0D0
        SUM = SUM + C(I) * F(Y)
      ENDDO

      GAUSS2 = SUM * (B-A) / 2.0D0

      END

      FUNCTION GAUSS3(A, B, F)
      DOUBLE PRECISION A, B, GAUSS3, F
      EXTERNAL F

      DOUBLE PRECISION X(3), C(3)
      DOUBLE PRECISION SUM, Y
      INTEGER I

      C(1) = 5.0D0 / 9.0D0
      C(2) = 8.0D0 / 9.0D0
      C(3) = 5.0D0 / 9.0D0

      X(1) = -SQRT(0.6D0)
      X(2) = 0.0D0
      X(3) =  SQRT(0.6D0)

      SUM = 0.0D0
      DO I=1,3
        Y = ( (B-A) * X(I) + (B+A) ) / 2.0D0
        SUM = SUM + C(I) * F(Y)
      ENDDO

      GAUSS3 = SUM * (B-A) / 2.0D0

      END
