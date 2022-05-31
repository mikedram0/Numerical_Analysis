SUBROUTINE jacobi(a, b, x)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in)    :: a(:,:)
  DOUBLE PRECISION, INTENT (in)    :: b(:)
  DOUBLE PRECISION, INTENT (inout) :: x(:)

  DOUBLE PRECISION :: y(SIZE(x))
  INTEGER :: i
  DOUBLE PRECISION, PARAMETER :: TOL = 1d-7
  
  DO
     y = b - MATMUL(a,x)

     DO i=1,SIZE(x)
        x(i) = x(i) + y(i) / a(i,i)
     ENDDO

     IF (ALL(ABS(y) < TOL)) EXIT
  ENDDO

END SUBROUTINE jacobi


SUBROUTINE seidel(a, b, x)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in)    :: a(:,:)
  DOUBLE PRECISION, INTENT (in)    :: b(:)
  DOUBLE PRECISION, INTENT (inout) :: x(:)

  DOUBLE PRECISION :: y(SIZE(x))
  INTEGER :: i
  DOUBLE PRECISION, PARAMETER :: TOL = 1d-7
  
  DO
     DO i=1,SIZE(x)
        y(i) = b(i) - DOT_PRODUCT(a(i, :), x) 
        
        x(i) = x(i) + y(i) / a(i,i)
     ENDDO

     IF (ALL(ABS(y) < TOL)) EXIT
  ENDDO

END SUBROUTINE seidel


PROGRAM gjs
  IMPLICIT NONE

  INTEGER, PARAMETER :: n = 4
  DOUBLE PRECISION :: a(n,n), b(n), x(n)

  INTERFACE 
     SUBROUTINE jacobi(a, b, x)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in)    :: a(:,:)
       DOUBLE PRECISION, INTENT (in)    :: b(:)
       DOUBLE PRECISION, INTENT (inout) :: x(:)
     END SUBROUTINE jacobi

     SUBROUTINE seidel(a, b, x)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in)    :: a(:,:)
       DOUBLE PRECISION, INTENT (in)    :: b(:)
       DOUBLE PRECISION, INTENT (inout) :: x(:)
     END SUBROUTINE seidel
  END INTERFACE


  a(1,1) = 12.1d0
  a(1,2) = 3.9d0
  a(1,3) = 0.3d0
  a(1,4) = -4.1d0

  a(2,1) = 4.3d0
  a(2,2) = -11.3d0
  a(2,3) = 0.8d0
  a(2,4) = 1.5d0

  a(3,1) = 1.0d0
  a(3,2) = -2.8d0
  a(3,3) = 14.3d0
  a(3,4) = -8.1d0

  a(4,1) = 2.4d0
  a(4,2) = 6.1d0
  a(4,3) = -1.1d0
  a(4,4) = 12.5d0

  b(1) = 1.2d0
  b(2) = 2.3d0
  b(3) = 3.4d0
  b(4) = 4.5d0

  x = 0.0d0
  CALL jacobi(a,b,x)
  PRINT *, "solution with jacobi:"
  PRINT *, x

  x = 0.0d0
  CALL seidel(a,b,x)
  PRINT *, "solution with seidel:"
  PRINT *, x

END PROGRAM gjs
