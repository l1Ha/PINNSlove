SUBROUTINE SmxInv (a,n)
USE Numeric_Kinds_Module
USE FileUnits_Module
! 
!	computes in-place inverse of symmetric matrix
!
!	using only lower half of matrix a
!
!	===Input===
!     n         -- size of matrix to be inverted
!     a(i,j)    -- matrix to be inverted
!
!	===Output==:
!     a(i,j)    -- the inverted matrix of a 
!
IMPLICIT NONE
INTEGER n, info, ipiv(n)
REAL(Kind=WP_Kind) a(n,n), ui(n,n), work(n*n)

CALL unit(n, ui)
CALL dsysv ('l', n, n, a, n, ipiv, ui, n, work, n*n, info)

IF(info /= 0)THEN
   WRITE(Msg_Unit,*)"Error occurred in SmxInv", info
   STOP "Error occurred in SmxInv"
ENDIF

a=ui
RETURN
ENDSUBROUTINE SmxInv
