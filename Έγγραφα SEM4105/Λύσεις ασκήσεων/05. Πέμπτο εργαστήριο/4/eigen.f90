PROGRAM eigen
  IMPLICIT NONE
  INTEGER, PARAMETER  :: n = 4
  DOUBLE PRECISION :: a(n,n)
  DOUBLE PRECISION :: x1, x2, x

  INTERFACE
     SUBROUTINE secant(x1, x2, x, toler, func, a)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (inout) :: x1, x2
       DOUBLE PRECISION, INTENT (out)   :: x
       DOUBLE PRECISION, INTENT (in)    :: toler
       DOUBLE PRECISION, INTENT (in)    :: a(:,:)

       INTERFACE
          FUNCTION func(a, x)
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT (in) :: a(:,:), x
            DOUBLE PRECISION :: func
          END FUNCTION func
       END INTERFACE
     END SUBROUTINE secant

     FUNCTION f(a,lambda)
       IMPLICIT NONE

       DOUBLE PRECISION, INTENT (in) :: a(:,:), lambda
       DOUBLE PRECISION :: f
     END FUNCTION f
  END INTERFACE

  A(1,1) = 2.1D0
  A(1,2) = 3.9D0
  A(1,3) = 0.3D0
  A(1,4) = -4.1D0
  A(2,1) = 4.3D0
  A(2,2) = -1.3D0
  A(2,3) = 0.8D0
  A(2,4) = 1.5D0
  A(3,1) = 1.0D0
  A(3,2) = -2.8D0
  A(3,3) = 4.3D0
  A(3,4) = -8.1D0
  A(4,1) = 2.4D0
  A(4,2) = 6.1D0
  A(4,3) = -1.1D0
  A(4,4) = 12.5D0

  x1 = -1.0d0
  x2 =  1.0d0

  CALL secant(x1, x2, x, 1d-8, f, a)

  PRINT *, "An eigenvalue is ",  x, "with determinant value ", f(a,x)
END PROGRAM eigen


SUBROUTINE secant(x1, x2, x, toler, func, a)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (inout) :: x1, x2
  DOUBLE PRECISION, INTENT (out)   :: x
  DOUBLE PRECISION, INTENT (in)    :: toler
  DOUBLE PRECISION, INTENT (in)    :: a(:,:)

  INTERFACE
     FUNCTION func(a, x)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: a(:,:), x
       DOUBLE PRECISION :: func
     END FUNCTION func
  END INTERFACE

  DOUBLE PRECISION :: f1, f2

  DO 
     f1 = func(a, x1)
     f2 = func(a, x2)

     x = x2 - f2 * (x2 - x1) / (f2 - f1)

     IF (ABS(func(a, x)) < toler) EXIT

     x1 = x2
     x2 = x
  END DO

END SUBROUTINE secant



! gives the determinant of a
FUNCTION det(a)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: a(:,:)
  DOUBLE PRECISION :: det

  DOUBLE PRECISION :: b(SIZE(a,1),SIZE(a,2))
  INTEGER :: i, changes

  INTERFACE
     SUBROUTINE matrix_low_triang(a, changes)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (inout) :: a(:,:)
       INTEGER, INTENT (out) :: changes
     END SUBROUTINE matrix_low_triang
  END INTERFACE

  b = a
  CALL matrix_low_triang(b, changes)
  ! now b is diagonal

  det = 1.0d0
  DO i= 1,SIZE(a,1)
     det = det * b(i,i)
  ENDDO

  IF (MOD(changes,2) /= 0) det = -det

END FUNCTION det

! make  A  lower triangular
SUBROUTINE matrix_low_triang(a, changes)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (inout) :: a(:,:)
  INTEGER, INTENT (out) :: changes

  INTERFACE
     SUBROUTINE pivot(a, i, changes)
       IMPLICIT NONE

       DOUBLE PRECISION, INTENT (inout) :: A(:,:)
       INTEGER, INTENT (in) :: i
       INTEGER, INTENT (inout) :: changes
     END SUBROUTINE pivot
  END INTERFACE

  DOUBLE PRECISION :: g
  INTEGER :: n, i, k

  n = SIZE(a,1)

  changes = 0

  !     lower triangle
  DO i=1,n-1
     CALL pivot(a, i, changes)
     DO k=i+1,n
        g = -a(k,i) / a(i,i)
        a(k,i:n) = a(k,i:n) + a(i,i:n) * g
     ENDDO
  ENDDO
END SUBROUTINE matrix_low_triang


!meriki odigisi (partial pivoting)
SUBROUTINE pivot(a, i, changes)
  IMPLICIT NONE

  DOUBLE PRECISION, INTENT (inout) :: A(:,:)
  INTEGER, INTENT (in) :: i
  INTEGER, INTENT (inout) :: changes

  INTERFACE
     ELEMENTAL SUBROUTINE swap(a,b)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (inout) :: a, b
     END SUBROUTINE swap
  END INTERFACE

  INTEGER :: j, k

  j = i
  DO k = i+1, SIZE(a,1)
     IF ( ABS(a(k,i)) > ABS(a(j,i)) ) j = k
  ENDDO
  ! in i column a(j,i) is maximum.

  IF (j > i) THEN 
     CALL swap(a(j,:), a(i,:))  ! exchange i, j rows
     changes = changes + 1
  END IF

END SUBROUTINE pivot

! a<->b
ELEMENTAL SUBROUTINE swap(a,b)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (inout) :: a, b

  DOUBLE PRECISION :: temp

  temp = a
  a    = b
  b    = temp

END SUBROUTINE swap



! function f(x) where x is lambda
FUNCTION f(a,lambda)
  IMPLICIT NONE

  DOUBLE PRECISION, INTENT (in) :: a(:,:), lambda
  DOUBLE PRECISION :: f

  INTERFACE
     FUNCTION det(a)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: a(:,:)
       DOUBLE PRECISION :: det
     END FUNCTION det
  END INTERFACE

  DOUBLE PRECISION :: b(SIZE(a,1),SIZE(a,2))
  INTEGER :: i

  ! set b = a - lambda i
  b = a
  DO i = 1, SIZE(a,1)
     b(i,i) = a(i,i) - lambda
  ENDDO

  f = det(b)

END FUNCTION f

