      PROGRAM DRVTV
      DOUBLE PRECISION XP(100), YP(100)
      DOUBLE PRECISION X, DER1, DER2, P1, P2
      DOUBLE PRECISION DERIV
      EXTERNAL DERIV, REDATA

      CALL REDATA("points.dat", XP, YP, N)

      PRINT *, "Dose shmeio ypologismoy paragwgwn"
      READ *, X
      
      DER1 = DERIV(1, XP, YP, X, N)
      DER2 = DERIV(2, XP, YP, X, N)
      P1 = DER1 - SIN(2D0*X)
      P2 = DER2 - 2D0*COS(2D0*X)
      
      PRINT *, "H proti paragogos einai ", DER1
      PRINT *, "H diafora apo tin akribi timi einai ", P1
      PRINT *
      PRINT *, "H deyteri paragogos einai ", DER2
      PRINT *, "H diafora apo tin akribi timi einai ", P2
      END 

      FUNCTION DERIV(ORDER, XP, YP, X, N)      
      INTEGER          ORDER, N
      DOUBLE PRECISION XP(100), YP(100)
      DOUBLE PRECISION X
      DOUBLE PRECISION DERIV

      INTEGER K,J
      DOUBLE PRECISION A(100,100), B(100), SOL(100)
      EXTERNAL TRIANG, BACKSU

      IF (ORDER .EQ. 1) THEN
         B(1) = 0.0d0           ! πρώτη παράγωγος του 1
         B(2) = 1.0d0           ! πρώτη παράγωγος του x     
         
         DO 10 K=3,N
            B(K) = (K-1) * X**(K-2) ! πρώτη παράγωγος του x**(k-1)
 10      CONTINUE
      ENDIF

      IF (ORDER.EQ.2) THEN
         B(1) = 0.0d0           ! δεύτερη παράγωγος του 1
         B(2) = 0.0d0           ! δεύτερη παράγωγος του x     
         B(3) = 2.0d0           ! δεύτερη παράγωγος του x^2
         
         DO 20 K=4,N
            B(K) = (K-1) * (K-2) * X**(K-3) ! δεύτερη παράγωγος του x**(k-1)
 20      CONTINUE
      ENDIF
      
      DO 30 J=1,N
         A(1,J) = 1.0D0
 30   CONTINUE
      
      DO 40 K=2,N
         DO 80 J=1,N
            A(K,J) = XP(J)**(K-1)
 80      CONTINUE
 40   CONTINUE
      
      CALL TRIANG(N, A, B)
      CALL BACKSU(N, A, B, SOL)

       DERIV = 0.0D0
       DO 50 J=1,N
          DERIV = DERIV + SOL(J) * YP(J)
 50    CONTINUE
       END 


C     Read data file
      SUBROUTINE REDATA(FNAME, X, Y, N)
      CHARACTER(*) FNAME
      DOUBLE PRECISION  X(100), Y(100)

      INTEGER UNITNO, I, N
      PARAMETER(UNITNO=55)

      OPEN (UNIT = UNITNO, FILE = FNAME, STATUS = 'OLD')
      READ (UNITNO, *) N
      
      DO 10 I=1,N
         READ (UNITNO, *) X(I), Y(I)
 10   CONTINUE
      
      CLOSE (UNITNO)
      END 




C     MAKE LOWER TRIANGULAR
      SUBROUTINE TRIANG(N, A, B)
      INTEGER N
      DOUBLE PRECISION A(100,100), B(100)
      
      INTEGER I, J, K
      DOUBLE PRECISION G
      EXTERNAL PIVOT

      DO 10 K=1,N-1
        CALL PIVOT(N, A, B, K)    
        IF (ABS(A(K,K)) .LT. 1D-8) GOTO 10
        
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
      DOUBLE PRECISION A(100,100), B(100)
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
      DOUBLE PRECISION A(100,100), B(100), X(100)
      
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
