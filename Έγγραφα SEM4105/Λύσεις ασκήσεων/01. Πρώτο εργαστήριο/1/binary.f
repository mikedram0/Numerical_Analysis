      PROGRAM BINARY
      INTEGER K, I, J, B(32)
      
      WRITE (*,*) "DOSE AKERAIO"
      READ (*,*) K
      
      DO 100 I = 1, 32
         B(I) = MOD(K,2)
      
         K = K / 2
         
         IF (K .EQ. 0) GOTO 10
 100  CONTINUE
      
 10   WRITE(*, "(32I1)") (B(J), J=I,1,-1)
      END

