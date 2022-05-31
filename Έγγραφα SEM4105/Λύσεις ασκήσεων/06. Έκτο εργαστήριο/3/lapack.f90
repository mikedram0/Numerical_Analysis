PROGRAM eigenvalues
  IMPLICIT NONE

  INTERFACE
     SUBROUTINE DGEEV(JOBVL,JOBVR,N,A,LDA,WR,WI,VL,LDVL,VR,LDVR,&
          WORK, LWORK,INFO)
       IMPLICIT NONE
       CHARACTER,        INTENT (in)    :: JOBVL,JOBVR
       INTEGER,          INTENT (in)    :: N, LDA, LDVL, LDVR, LWORK
       INTEGER,          INTENT (out)   :: INFO       
       DOUBLE PRECISION, INTENT (inout) :: A(LDA,*)
       DOUBLE PRECISION, INTENT (out)   :: WR(*), WI(*)
       DOUBLE PRECISION, INTENT (out)   :: VL(LDVL,*), VR(LDVR,*), WORK(*)
     END SUBROUTINE DGEEV
  END INTERFACE

  integer, parameter :: n = 3
  integer, parameter :: ldvl = 1
  integer, parameter :: ldvr = 1

  DOUBLE PRECISION :: a(n,n), wr(n), wi(n), vl(ldvl,n), vr(ldvr,n), wrk(1)
  DOUBLE PRECISION, ALLOCATABLE :: work(:)
  INTEGER :: info,lwork
  
  INTEGER :: i
  
  a(1,1) = 6.3d0
  a(1,2) = 2.10d0
  a(1,3) = 4.15d0
  a(2,1) = 3.1d0
  a(2,2) = 5.14d0
  a(2,3) = 1.03d0
  a(3,1) = -11d0
  a(3,2) = 12.3d0
  a(3,3) = -8.8d0

  lwork = -1
  info = 0
  call dgeev('n', 'n', n, a, n, wr, wi, vl, ldvl, vr, ldvr, wrk, lwork, info)

  lwork = INT(wrk(1))
  ALLOCATE(work(lwork))
  info = 0
  call dgeev('n', 'n', n, a, n, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
  DEALLOCATE (work)

  IF (info /= 0) THEN
     PRINT *, "error"
     STOP
  END IF
  
  DO i =1, N
     PRINT *, WR(i), WI(i)
  END DO
END PROGRAM eigenvalues
