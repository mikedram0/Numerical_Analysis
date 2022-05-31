      PROGRAM MAINPR

      DOUBLE PRECISION RICHAR, G, TRAPEZ
      EXTERNAL RICHAR, G, TRAPEZ
      DOUBLE PRECISION PI
      PARAMETER (PI = 3.14159265358979323846D0)
      DOUBLE PRECISION A, B
      PARAMETER (A=0.0D0, B = PI)
      INTEGER N
      PARAMETER (N=16)
      DOUBLE PRECISION CORR
      PARAMETER (CORR=2.0D0)

      PRINT *, "     Spaces      Error with TRAPEZ"
      PRINT *,   N, CORR - TRAPEZ(G, A, B,   N)
      PRINT *, 2*N, CORR - TRAPEZ(G, A, B, 2*N)
      PRINT *, 4*N, CORR - TRAPEZ(G, A, B, 4*N)

      PRINT *
      PRINT *, "Error with Richardson: ", CORR - RICHAR(G, A, B, N)

      END

      FUNCTION G(X)
      DOUBLE PRECISION X, G
      G = SIN(X)
      END

C     RICHARDSON INTEGRATION: THREE STEPS H, H/2, H/4
C     
C     I = (I_H - 20 I_{H/2} + 64 I_{H/4}) / 45
C
      FUNCTION RICHAR(F, A, B, N)

      DOUBLE PRECISION A, B, RICHAR
      INTEGER N
      DOUBLE PRECISION F
      EXTERNAL F

      DOUBLE PRECISION TRAPEZ
      EXTERNAL TRAPEZ
      DOUBLE PRECISION IH, IH2, IH4

      IH  = TRAPEZ(F, A, B,   N)
      IH2 = TRAPEZ(F, A, B, 2*N)
      IH4 = TRAPEZ(F, A, B, 4*N)

      RICHAR = (IH - 20.0D0 * IH2 + 64.0D0 * IH4) / 45.0D0

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

      DO 10 I= 1,N-1
        X = A + I * H

        SUM = SUM + F(X)
 10   CONTINUE

      TRAPEZ = SUM * H
      END
