      SUBROUTINE JACOBI(N, A, B, X)
      INTEGER N
      DOUBLE PRECISION A(N,N), B(N), X(N)

      DOUBLE PRECISION Y(N)
      INTEGER I, J
      LOGICAL FOUND
      DOUBLE PRECISION TOL
      PARAMETER (TOL=1D-7)

 10   FOUND = .TRUE.

      DO I=1,N
        Y(I) = B(I)
        DO J = 1,N
          Y(I) = Y(I) - A(I,J) * X(J)
        ENDDO
      ENDDO

      DO I=1,N
        FOUND = FOUND .AND. (ABS(Y(I)) < TOL)

        X(I) = X(I) + Y(I) / A(I,I)
      ENDDO

      IF (FOUND) GOTO 20
      GOTO 10

 20   END 


      SUBROUTINE SEIDEL(N, A, B, X)
      INTEGER N
      DOUBLE PRECISION A(N,N), B(N), X(N)

      DOUBLE PRECISION Y
      INTEGER I, J
      LOGICAL FOUND
      DOUBLE PRECISION TOL
      PARAMETER (TOL=1D-7)

 10   FOUND = .TRUE.

      DO I=1,N
        Y = B(I)
        DO J = 1,N
          Y = Y - A(I,J) * X(J)
        ENDDO

        FOUND = FOUND .AND. (ABS(Y) < TOL)

        X(I) = X(I) + Y / A(I,I)
      ENDDO

      IF (FOUND) GOTO 20
      GOTO 10

 20   END 


      PROGRAM GJS
      INTEGER N
      PARAMETER (N=4)

      INTEGER I
      DOUBLE PRECISION A(N,N), B(N), X(N)

      A(1,1) = 12.1D0
      A(1,2) = 3.9D0
      A(1,3) = 0.3D0
      A(1,4) = -4.1D0
      
      A(2,1) = 4.3D0
      A(2,2) = -11.3D0
      A(2,3) = 0.8D0
      A(2,4) = 1.5D0
      
      A(3,1) = 1.0D0
      A(3,2) = -2.8D0
      A(3,3) = 14.3D0
      A(3,4) = -8.1D0
      
      A(4,1) = 2.4D0
      A(4,2) = 6.1D0
      A(4,3) = -1.1D0
      A(4,4) = 12.5D0
      
      B(1) = 1.2D0
      B(2) = 2.3D0
      B(3) = 3.4D0
      B(4) = 4.5D0
      

      DO I=1,N
        X(I) = 0.0D0
      ENDDO

      CALL JACOBI(N,A,B,X)
      PRINT *, "solution with jacobi:"
      PRINT *, (X(I), I=1,N)

      DO I=1,N
        X(I) = 0.0D0
      ENDDO

      CALL SEIDEL(N,A,B,X)
      PRINT *, "solution with seidel:"
      PRINT *, (X(I), I=1,N)

      END
