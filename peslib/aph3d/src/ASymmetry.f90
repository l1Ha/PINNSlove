SUBROUTINE ASymmetry(a,n,err)
USE Numeric_Kinds_Module
!-----------------------------------------------------------------------
! Author: Gregory A. Parker
! Calculates the asymmetry of the matrix A
IMPLICIT NONE
INTEGER n,j,i
CHARACTER(LEN=21), PARAMETER:: ProcName='ASymmetry'
CHARACTER(LEN=6) Print_Flag
REAL(Kind=WP_Kind) a(n,n), vnum, vden, err

CALL PoptAsy(ProcName, Print_Flag)
!
vnum = 0.d0
vden = 0.d0
DO j=1,n
   DO i=j+1,n
      vnum = vnum + (a(j,i)-a(i,j))**2
      vden = vden + (abs(a(j,i))+abs(a(i,j)))**2
   ENDDO
ENDDO
err = SQRT(vnum/vden)
! If err=0 the matrix is exaclty symmetric
RETURN
ENDSUBROUTINE ASymmetry