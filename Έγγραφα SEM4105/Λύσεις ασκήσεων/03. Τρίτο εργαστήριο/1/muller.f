C Find the roots of f(x)=0 with f(x)=sin(x) - x^2
C using the Muller's method.
C use complex arithmetic

      PROGRAM MULLER
      DOUBLE PRECISION EPSI, DELTA
      PARAMETER (EPSI = 1D-6, DELTA = 1D-6)
      
      DOUBLE COMPLEX X,X0,X1,X2, W1, W0,A,B,C, D, D1, D2, P
      DOUBLE COMPLEX F0,F1,F2, T
      DOUBLE COMPLEX F

      X0 = (0.1D0, 0.0D0)
      X1 = (0.2D0, 0.0D0)
      X2 = (0.5D0, 0.0D0)

      F0 = F(X0)
      F1 = F(X1)
      F2 = F(X2)

 10   W0 = (F2-F0)/(X2-X0)
      W1 = (F2-F1)/(X2-X1)      
      
      A = (W1-W0)/(X1-X0)
      
      B = W0+A*(X2-X0)
      C = F2
      
      P = B*B-4*A*C

      D1 = B + SQRT(P)
      D2 = B - SQRT(P)
      
      IF (ABS(D1) > ABS(D2)) THEN
         D = D1
      ELSE
         D = D2
      ENDIF
      
      X = X2 - 2*C / D

      T = F(X)
C     termination condition
      IF((ABS(X - X2) .LT. DELTA) .OR. (ABS(T) .LT. EPSI)) GOTO 20

      X0 = X1
      X1 = X2
      X2 = X
      F0 = F1
      F1 = F2
      F2 = T
      
      GOTO 10

 20   WRITE (*,*) "Η ρίζα είναι", X 
      WRITE (*,*) "Η τιμή της συνάρτησης", F(X)

      END 

C     Define f(x)
      FUNCTION F(X)
      DOUBLE COMPLEX  X, F
      F = SIN(X) - X**2
      END 

