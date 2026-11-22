SUBROUTINE asym(a,n,err)
USE Numeric_Kinds_Module
USE Numbers_Module
!-----------------------------------------------------------------------
!
IMPLICIT NONE
INTEGER n,j,i
REAL(Kind=WP_Kind) a(n,n), vnum, vden, err
!
vnum = 0.d0
vden = 0.d0
DO j=1,n
   DO i=j+1,n
      vnum = vnum + (a(j,i)-a(i,j))**2
      vden = vden + (abs(a(j,i))+abs(a(i,j)))**2
   ENDDO
ENDDO
IF(vnum==Zero)THEN
   err=Zero
ELSE
   err = SQRT(vnum/vden)
ENDIF
!
RETURN
ENDSUBROUTINE asym
