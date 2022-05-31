C R(x_i) = y_i => P(x_i) = y_i Q(x_i) =>
C a_0 + a_1 x_i + a_2 x_i^2 = y_i (1 + b_1 x_i + b_2 x_i^2 + b_3 x_i^3 =>
C 1 a_0 + x_i a_1 + x_i^2 a_2 - y_i x_i b_1 - y_i x_i^2 b_2 - y_i x_i^3 b_3 = y_i
C
C  | 1  x_1   x_1^2   -y_1*x_1   -y_1*x_1^2  -y_1*x_1^3 |  | a_0 |    | y_1 |
C  | 1  x_2   x_2^2   -y_2*x_2   -y_2*x_2^2  -y_2*x_2^3 |  | a_1 |    | y_2 |
C  | 1  x_3   x_3^2   -y_3*x_3   -y_3*x_3^2  -y_3*x_3^3 |  | a_2 |    | y_3 |
C  | 1  x_4   x_4^2   -y_4*x_4   -y_4*x_4^2  -y_4*x_4^3 |* | b_1 | =  | y_4 |
C  | 1  x_5   x_5^2   -y_5*x_5   -y_5*x_5^2  -y_5*x_5^3 |  | b_2 |    | y_5 |
C  | 1  x_6   x_6^2   -y_6*x_6   -y_6*x_6^2  -y_6*x_6^3 |  | b_3 |    | y_6 |
C
      PROGRAM RATION
      
      INTEGER N
      PARAMETER (N = 6)
      DOUBLE PRECISION  A(N,N), B(N), X(N), Y(N), C(N)
      
      DATA X / 0.9D0, 1.1D0, 1.5D0, 2.0D0, 2.9D0, 3.5D0 /
      DATA Y / 5.607D0, 4.576D0, 3.726D0, 3.354D0, 3.14D0, 3.087D0 /

      INTEGER I
      DOUBLE PRECISION R
      EXTERNAL R, TRIANG,BACKSU

      DO 10 I = 1,N
        A(I,1) = 1D0
        A(I,2) = X(I)
        A(I,3) = X(I)**2
        A(I,4) = -Y(I)*X(I)
        A(I,5) = -Y(I)*X(I)**2
        A(I,6) = -Y(I)*X(I)**3
        
        B(I) = Y(I)
 10   CONTINUE


C SOLVE A*C=B TO GET C, THE COEFFICIENTS a,b
  
      CALL TRIANG(N,A,B)
      CALL BACKSU(N,A,B,C)

      PRINT *, "A_0  A_1  A_2 = ", C(1), C(2), C(3)
      PRINT *, "B_1  B_2  B_3 = ", C(4), C(5), C(6)

      
      DO 20 I=1,N
        PRINT *, R(C, X(I)), Y(I)
 20   CONTINUE
      

      END


      
      FUNCTION R(C, X)
      DOUBLE PRECISION C(6)
      DOUBLE PRECISION X
      DOUBLE PRECISION R  
      DOUBLE PRECISION P, Q
      
      P = C(1) + C(2) * X + C(3) *X*X
      Q = 1D0 + C(4) * X + C(5) * X*X + C(6) * X*X*X
      
      R = P/Q
  
      END


C     MAKE LOWER TRIANGULAR
      SUBROUTINE TRIANG(N, A, B)
      INTEGER N
      DOUBLE PRECISION A(N,N), B(N)

      INTEGER I, J, K
      DOUBLE PRECISION G
      EXTERNAL PIVOT

      DO 10 K=1,N-1
        CALL PIVOT(N, A, B, K)    

        DO 20 I=K+1,N
          G = -A(I,K) / A(K,K)

          DO 30 J = 1,N
            A(I,J) = A(I,J) + A(K,J) * G
 30       CONTINUE
          B(I) = B(I) + B(K) * G
 20     CONTINUE
 10   CONTINUE

      END


C     MERIKI ODIGISI (PARTIAL PIVOTING)

      SUBROUTINE PIVOT(N, A, B, K)
      INTEGER N
      DOUBLE PRECISION A(N,N), B(N)
      INTEGER K
      
      INTEGER I, J
      DOUBLE PRECISION MAXV
      EXTERNAL SWAP

C     'NORMALIZE' EQUATIONS FROM K AND BELOW
      DO 10 I = K, N

        MAXV = ABS(B(I))
        DO 20 J = 1, N
          IF (MAXV .LT. ABS(A(I,J))) MAXV = ABS(A(I,J))
 20     CONTINUE

        DO 30 J = 1, N
          A(I,J) = A(I,J) / MAXV
 30     CONTINUE
        B(I) = B(I) / MAXV

 10   CONTINUE
      

C     FIND POSITION OF MAX ELEMENT IN K COLUMN BELOW THE DIAGONAL.
      I = K
      DO 40 J = K+1, N
        IF ( ABS(A(I,K)) .GT. ABS(A(J,K)) ) I = J
 40   CONTINUE
C     IN I COLUMN A(I,K) IS MAXIMUM.
      
      IF (I .NE. K) THEN              
C     EXCHANGE I, K ROWS OF A AND B
         DO 50 J=1,N
           CALL SWAP(A(I,J), A(K,J))
 50      CONTINUE
         CALL SWAP(B(I), B(K))    
      END IF
      END 

C     A<->B
      SUBROUTINE SWAP(A,B)
      DOUBLE PRECISION A, B
      DOUBLE PRECISION TEMP
      
      TEMP = A
      A    = B
      B    = TEMP
      END



C     BACKSUBSTITUTION
      SUBROUTINE BACKSU(N, A, B, X)
      INTEGER N
      DOUBLE PRECISION A(N,N), B(N), X(N)
      
      DOUBLE PRECISION SUM
      INTEGER K, J
      
      DO 10 K=N,1,-1
        SUM = 0.0D0
        DO 20 J=K+1,N
          SUM = SUM + A(K,J) * X(J)
 20     CONTINUE
        X(K) = (B(K) - SUM) / A(K,K)
 10   CONTINUE

      END
