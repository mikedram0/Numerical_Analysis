FUNCTION simpson(f, a, b, n)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) ::  a, b
  INTEGER, INTENT (in) :: n
  DOUBLE PRECISION :: simpson

  INTERFACE
     FUNCTION f(x)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: x
       DOUBLE PRECISION :: f
     END FUNCTION f
  END INTERFACE

  INTEGER :: i
  DOUBLE PRECISION :: x, sum1, sum2, h, sum


  IF (MOD(n,2)/=0) STOP

  h = (b-a) / n

  sum1 = 0.0d0
  DO i= 1,n-1,2
     x = a + i * h

     sum1 = sum1 + f(x)
  ENDDO

  sum2 = 0.0d0
  DO i= 2,n-2,2
     x = a + i * h

     sum2 = sum2 + f(x)
  ENDDO

  sum = f(a) + f(b) + 4.0d0 * sum1 + 2.0d0 * sum2

  simpson = sum * h / 3.0d0
END FUNCTION simpson


FUNCTION f(x)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: x
  DOUBLE PRECISION :: f

  f = SIN(x)
END FUNCTION f


PROGRAM simpson_integration
  IMPLICIT NONE

  DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979323846d0
  DOUBLE PRECISION :: int
  INTEGER :: n, k
  INTERFACE
     FUNCTION simpson(f, a, b, n)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) ::  a, b
       INTEGER, INTENT (in) :: n
       DOUBLE PRECISION :: simpson

       INTERFACE
          FUNCTION f(x)
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT (in) :: x
            DOUBLE PRECISION :: f
          END FUNCTION f
       END INTERFACE
     END FUNCTION simpson


     FUNCTION f(x)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: x
       DOUBLE PRECISION :: f
     END FUNCTION f

  END INTERFACE
  PRINT *, "ν                τιμη                      σφαλμα" 

  DO k = 1, 9
     n = 2**k
     int = simpson(f, 0.0d0, pi, n)
     PRINT *, n, int, 2.0d0-int
  ENDDO

END PROGRAM simpson_integration
