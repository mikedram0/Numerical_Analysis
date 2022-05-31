PROGRAM least_squares
  IMPLICIT NONE
  INTEGER :: method
  DOUBLE PRECISION, ALLOCATABLE :: x(:), y(:)
  DOUBLE PRECISION :: A1, A0, R2

  INTERFACE
     SUBROUTINE read_data(fname,  x, y)
       IMPLICIT NONE
       CHARACTER(*), INTENT (in) :: fname
       DOUBLE PRECISION, ALLOCATABLE, INTENT (out) :: x(:), y(:)
     END SUBROUTINE read_data

     SUBROUTINE LSQLIN(X, Y, A1, A0, R2)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: X(:), Y(:)
       DOUBLE PRECISION, INTENT (out) :: A1, A0, R2
     END SUBROUTINE LSQLIN

     SUBROUTINE LSQPOW(X, Y, A1, A0, R2)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: X(:), Y(:)
       DOUBLE PRECISION, INTENT (out) :: A1, A0, R2
     END SUBROUTINE LSQPOW

     SUBROUTINE LSQEXP(X, Y, A1, A0, R2)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: X(:), Y(:)
       DOUBLE PRECISION, INTENT (out) :: A1, A0, R2
     END SUBROUTINE LSQEXP

     SUBROUTINE LSQLOG(X, Y, A1, A0, R2)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: X(:), Y(:)
       DOUBLE PRECISION, INTENT (out) :: A1, A0, R2
     END SUBROUTINE LSQLOG
  END INTERFACE

  CALL read_data("points.dat", X, Y)

  PRINT *, "Μέθοδος: "
  PRINT *, "1) y = a * x + b"
  PRINT *, "2) y = a * x^b"
  PRINT *, "3) y = a + b * exp(x)"
  PRINT *, "4) y = a + b * ln(x)"
  PRINT *, "Δώσε μέθοδο."

  READ (*,*) METHOD

  SELECT CASE (method)
  CASE (1)
     CALL LSQLIN(X, Y, A1, A0, R2)
     PRINT *, "y = ", A1, " x + ", A0
  CASE (2)
     CALL LSQPOW(X, Y, A1, A0, R2)
     PRINT *, "y = ", A0, " * x^", A1
  CASE (3)
     CALL LSQEXP(X, Y, A1, A0, R2)
     PRINT *, "y = ", A0, " + ", A1, " * exp(x)"
  CASE (4)
     CALL LSQLOG(X, Y, A1, A0, R2)
     PRINT *, "y = ", A0, " + ", A1, " * ln(x)"
  END SELECT

  PRINT *, "r^2 = ", R2

END PROGRAM least_squares

!     Read data file
SUBROUTINE read_data(fname, x, y)
  IMPLICIT NONE
  CHARACTER(*), INTENT (in) :: fname
  DOUBLE PRECISION, ALLOCATABLE, INTENT (out) :: x(:), y(:)

  INTEGER, PARAMETER ::  UNITNO = 54
  INTEGER :: n, i

  OPEN (UNIT = UNITNO, FILE = FNAME, status = 'old', action = 'read')
  READ (UNITNO, *) N

  ALLOCATE(X(n), Y(n))

  DO I=1,N
     READ (UNITNO, *) X(I), Y(I)
  ENDDO

  CLOSE (UNITNO)
END SUBROUTINE read_data


!     Y = A1 * X + A0
SUBROUTINE LSQLIN(X, Y, A1, A0, R2)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: X(:), Y(:)
  DOUBLE PRECISION, INTENT (out) :: A1, A0, R2

  DOUBLE PRECISION :: SX, SY, SX2, SY2, SXY, D
  INTEGER :: n

  n = SIZE(x)

  SX  = SUM(x)
  SY  = SUM(y)
  SX2 = SUM(x*x)
  SY2 = SUM(y*y)
  SXY = SUM(x*y)

  D = N * SX2 - SX * SX

  A1 = (N * SXY - SX * SY) / D 
  A0 = (SY - A1 * SX) / N

  R2 = (N * SXY - SX * SY)**2 / D / (N * SY2 - SY**2)
END SUBROUTINE LSQLIN


!     Y = A0 X^A1 => LN(Y) = LN(A0) + A1 * LN(X)
SUBROUTINE LSQPOW(X, Y, A1, A0, R2)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: X(:), Y(:)
  DOUBLE PRECISION, INTENT (out) :: A1, A0, R2

  DOUBLE PRECISION :: ZX(SIZE(x)), ZY(SIZE(x)), ZA0

  INTERFACE
     SUBROUTINE LSQLIN(X, Y, A1, A0, R2)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: X(:), Y(:)
       DOUBLE PRECISION, INTENT (out) :: A1, A0, R2
     END SUBROUTINE LSQLIN
  END INTERFACE

  IF ( ANY( (X <= 0.D0) .OR. (Y <= 0.D0) ) ) THEN
     PRINT *, "Data not appropriate"
     RETURN
  ENDIF
  
  ZX = LOG(X)
  ZY = LOG(Y)
  
  CALL LSQLIN(ZX, ZY, A1, ZA0, R2)
  
  A0 = EXP(ZA0)

END SUBROUTINE LSQPOW


!     Y = A0 + A1 * EXP(X)
SUBROUTINE LSQEXP(X, Y, A1, A0, R2)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: X(:), Y(:)
  DOUBLE PRECISION, INTENT (out) :: A1, A0, R2

  DOUBLE PRECISION :: ZX(SIZE(x))

  INTERFACE
     SUBROUTINE LSQLIN(X, Y, A1, A0, R2)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: X(:), Y(:)
       DOUBLE PRECISION, INTENT (out) :: A1, A0, R2
     END SUBROUTINE LSQLIN
  END INTERFACE

  ZX = EXP(X)

  CALL LSQLIN(ZX, Y, A1, A0, R2)

END SUBROUTINE LSQEXP


!     Y = A0 + A1 * LOG(X)
SUBROUTINE LSQLOG(X, Y, A1, A0, R2)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: X(:), Y(:)
  DOUBLE PRECISION, INTENT (out) :: A1, A0, R2

  DOUBLE PRECISION :: ZX(SIZE(x))

  INTERFACE
     SUBROUTINE LSQLIN(X, Y, A1, A0, R2)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: X(:), Y(:)
       DOUBLE PRECISION, INTENT (out) :: A1, A0, R2
     END SUBROUTINE LSQLIN
  END INTERFACE

  IF ( ANY(X <= 0.D0) ) THEN
     PRINT *, "Data not appropriate"
     RETURN
  ENDIF

  ZX = LOG(X)

  CALL LSQLIN(ZX, Y, A1, A0, R2)

END SUBROUTINE LSQLOG
