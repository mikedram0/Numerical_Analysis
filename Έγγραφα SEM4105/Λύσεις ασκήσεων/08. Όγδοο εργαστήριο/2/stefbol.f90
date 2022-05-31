PROGRAM stefan_boltzmann
  IMPLICIT NONE

  DOUBLE PRECISION :: A1, A0, EMBADO, R2
  DOUBLE PRECISION, ALLOCATABLE :: x(:), y(:)
  INTERFACE 
     !     Read data file
     SUBROUTINE read_data(fname,  x, y)
       IMPLICIT NONE
       CHARACTER(*), INTENT (in) :: fname
       DOUBLE PRECISION, ALLOCATABLE, INTENT (out) :: x(:), y(:)
     END SUBROUTINE read_data

     SUBROUTINE LSQPOW(X, Y, A1, A0, R2)
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT (in) :: X(:), Y(:)
       DOUBLE PRECISION, INTENT (out) :: A1, A0, R2
     END SUBROUTINE LSQPOW
  END INTERFACE

  CALL read_data("stefbol.dat", X, Y)

  !     Y = A0 X^A1
  CALL LSQPOW(X, Y, A1, A0, R2)

  EMBADO = 0.05D-4

  PRINT *, "Σταθερά Stefan-Boltzmann: ", A0 / EMBADO, "J/K^4/m^2/s"
  PRINT *, "Εκθέτης: ",  A1
  PRINT * , "r^2 = ", R2
END PROGRAM stefan_boltzmann


!     Read data file
SUBROUTINE read_data(fname,  x, y)
  IMPLICIT NONE
  CHARACTER(*), INTENT (in) :: fname
  DOUBLE PRECISION, ALLOCATABLE, INTENT (out) :: x(:), y(:)

  INTEGER, PARAMETER ::  UNITNO = 54
  INTEGER :: n, i

  OPEN (UNIT = UNITNO, FILE = FNAME, status = 'old', action = 'read')
  READ (UNITNO, *) N

  ALLOCATE(x(n), y(n))

  DO I=1,N
     READ (UNITNO, *) x(I), y(I)
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
  A0 = (SX2 * SY - SX * SXY) / D 

  R2 = (N * SXY - SX * SY)**2 / (N * SX2 - SX**2) / (N * SY2 - SY**2)
END SUBROUTINE LSQLIN


!     Y = A0 X^A1 => LN(Y) = LN(A0) + A1 * LN(X)
SUBROUTINE LSQPOW(X, Y, A1, A0, R2)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (in) :: X(:), Y(:)
  DOUBLE PRECISION, INTENT (out) :: A1, A0, R2

  INTEGER :: I
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
