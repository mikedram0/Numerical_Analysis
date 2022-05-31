      PROGRAM TRAP1
      
      DOUBLE PRECISION TRAPEZ, G
      EXTERNAL TRAPEZ, G
      
      DOUBLE PRECISION PI
      PARAMETER (PI = 3.14159265358979323846D0)
      DOUBLE PRECISION INTEGR
      INTEGER N, K

      DO K = 1, 9
        N = 2**K
        INTEGR = TRAPEZ(G, 0.0D0, PI, N)
        PRINT *, N, INTEGR, 2.0D0-INTEGR
      ENDDO

      END

      FUNCTION G(X)
      DOUBLE PRECISION X, G
      G = SIN(X)
      END


      FUNCTION TRAPEZ(F, A, B, N)

      DOUBLE PRECISION A, B, TRAPEZ
      INTEGER N
      DOUBLE PRECISION F
      EXTERNAL F

      INTEGER I
      DOUBLE PRECISION X, SUM, H

      H = (B-A) / N

      SUM = 0.5D0 * (F(A) + F(B)) 

      DO I= 1,N-1
        X = A + I * H

        SUM = SUM + F(X)
      ENDDO

      TRAPEZ = SUM * H
      END
