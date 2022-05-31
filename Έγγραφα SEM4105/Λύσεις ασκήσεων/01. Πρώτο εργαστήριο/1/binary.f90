PROGRAM binary
  IMPLICIT NONE
  INTEGER :: k, b(32), i, j

  WRITE (*,*) "Dose mh arnhtiko akeraio"
  READ (*,*) k

  DO i = 1,SIZE(b)
     b(i) = MOD(k,2)
     k = k/2
     IF (k == 0) EXIT
  ENDDO

  DO j=i,1,-1
     WRITE(*, "(I1)", advance="no") b(j)
  END DO

  WRITE(*,*)
  
END PROGRAM binary

