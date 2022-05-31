      SUBROUTINE RK4(X, Y, H, F, XNEW, YNEW, N)
      INTEGER N
      DOUBLE PRECISION X, Y(N), H, XNEW, YNEW(N), F
      
      DOUBLE PRECISION K(4,N), TEMP(N)
      INTEGER I, J
      
      DO 10 J = 1,N
         TEMP(J) = Y(J)
 10   CONTINUE
      
      DO 20 I = 1,N
        K(1,I) = H * F(I,X,TEMP,N)
 20   CONTINUE

      DO 30 J = 1,N
         TEMP(J) = Y(J) + K(1,J)/2.0D0
 30   CONTINUE
      
      DO 40 I = 1,N
        K(2,I) = H * F(I,X + H/2.0D0, TEMP,N)
 40   CONTINUE

      DO 50 J = 1,N
         TEMP(J) = Y(J) + K(2,J)/2.0D0
 50   CONTINUE

      DO 60 I = 1,N
         K(3,I) = H * F(I,X + H/2.0D0, TEMP,N)
 60   CONTINUE

      DO 70 J = 1,N
         TEMP(J) = Y(J) + K(3,J)
 70   CONTINUE

      DO 80 I = 1,N
        K(4,I) = H * F(I,X + H, TEMP,N)
 80   CONTINUE
      
      DO 90 I = 1,N
       YNEW(I) = Y(I) 
     &          + (K(1,I) + 2.0D0 * (K(2,I) + K(3,I)) + K(4,I)) / 6.0D0
 90   CONTINUE
      
      XNEW = X + H
      END



C     THE SYSTEM IS 
C     D THETA / D T  = Z
C     D Z / D T      = -SIN(THETA)
C     

C     THETA -> Y(1), Z-> Y(2).  D THETA / D T -> F(1)
C     D Z / D T -> F(2). T -> X
      FUNCTION F(I, X, Y, N)
      
      INTEGER I, N
      DOUBLE PRECISION X, Y(N)
      DOUBLE PRECISION  F
      
      IF (I .EQ. 1) F = Y(2)
      IF (I .EQ. 2) F = -SIN(Y(1)) 
      END 


C     SOLUTION OF THETA'' = - THETA:   THETA = A COS(T) + B SIN(T)
C     THETA(0) = 45 DEG. , THETA'(0) = 0  
      SUBROUTINE SOLUT(T,Y,N)
      INTEGER N
      DOUBLE PRECISION T, Y(N)

      DOUBLE PRECISION PI, A
      PARAMETER (PI = 3.14159265358979323846D0)
      PARAMETER (A = 45.0D0 / 180.0D0 * PI)

      Y(1) =  A * COS(T)
      Y(2) = -A * SIN(T)
      END 


      PROGRAM PNDLM

      DOUBLE PRECISION TINIT, TFIN, H
      DOUBLE PRECISION YINIT(2), YOLD(2), YNEW(2), YAPPRO(2)
      DOUBLE PRECISION TOLD, TNEW, F
      EXTERNAL F

      DOUBLE PRECISION PI
      PARAMETER (PI = 3.14159265358979323846D0)

      TINIT = 0.0D0
      TFIN  = 10.0D0
      
      H= 1D-5
      YINIT(1) = 45.0D0 / 180.0D0 * PI
      YINIT(2) = 0.0D0
      
      TOLD    = TINIT
      YOLD(1) = YINIT(1)
      YOLD(2) = YINIT(2)

 10   CALL RK4(TOLD,YOLD,H,F,TNEW,YNEW,2)
      TOLD    = TNEW
      YOLD(1) = YNEW(1)
      YOLD(2) = YNEW(2)
      
      CALL SOLUT(TNEW, YAPPRO, 2)
      
      PRINT "(5F16.8)", TNEW, YNEW(1), YNEW(2), YAPPRO(1), YAPPRO(2)
      
      IF (TNEW .LT. TFIN) GOTO 10
      
      END 
