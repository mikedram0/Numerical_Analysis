ELEMENTAL FUNCTION isPowerOf2(n)
  IMPLICIT NONE
  INTEGER, INTENT (in) :: n
  LOGICAL :: isPowerOf2

  isPowerof2 = (n /= 0) .AND. (IAND(n,n-1) /= 0) 
END FUNCTION isPowerOf2

! the m-th fourier coefficient 
RECURSIVE SUBROUTINE fft(m, f, C)
  IMPLICIT NONE
  INTEGER, INTENT (in) :: m
  DOUBLE COMPLEX, INTENT (in) :: f(0:)
  DOUBLE COMPLEX, INTENT (out) :: C

  DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979323846d0
  
  INTEGER :: n
  DOUBLE COMPLEX :: Ce, Co
  DOUBLE COMPLEX :: a
  
  n = SIZE(f)
  
  IF (.NOT.isPowerof2(n)) STOP

  IF (n == 1) THEN
     C = f(0)
     RETURN
  END IF

  CALL fft(m, f(0:n-1:2), Ce)
  CALL fft(m, f(1:n-1:2), Co)

  a = EXP(-2.d0 * pi * m * (0d0, 1.0d0) / n)

  C = (Ce + a * Co) / 2d0
END SUBROUTINE fft
  
