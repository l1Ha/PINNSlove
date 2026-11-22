SUBROUTINE lineqd (c, a, n, m)
USE Numeric_Kinds_Module
      USE FileUnits_Module
      INTEGER n, m, i, ip1, j, k
      REAL(Kind=WP_Kind) a(n,n), c(n,m), Temp
!-----------------------------------------------------------------------
!  solves a linear system of equations ax=c, the solution
!  is returned in the matrix or vector c.
!-----------------------------------------------------------------------
DO i=1, n
   ip1=i+1
   temp=1.d0/a(i,i)
   !IF(a(i,i)==0.d0)THEN
    !   STOP "HERE"
   !ENDIF
   IF(n == i) GOTO 20
   DO j=ip1, n
      a(i,j)=a(i,j)*temp
   ENDDO
20 DO j=1, m
      c(i,j)=c(i,j)*temp
   ENDDO
   DO 70 k=1, n
      IF(k == i) GOTO 70
      IF(n == i) GOTO 50
      DO j=ip1, n
         a(k,j)=a(k,j)-a(k,i)*a(i,j)
      ENDDO
50    DO j=1, m
         c(k,j)=c(k,j)-a(k,i)*c(i,j)
      ENDDO
70 CONTINUE
ENDDO
RETURN
END
