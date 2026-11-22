SUBROUTINE delchk(naph, ucheck, ndel, umat)
!
!-----------------------------------------------------------------
!     This routine provides a unitarity check on the APH to Delves
!     transformation matrix. If the magnitude of the diagonal elements
!     are almost 1 this indicates that the Delves basis nearly
!     complete (sufficiently large enough to expand the APH functions).
!-----------------------------------------------------------------
!this routine is called by: AphDel_Old or AphDel_New
!this routine calls: popt,mxoutl
!-----------------------------------------------------------------------
USE FileUnits_Module
USE crays_Module
USE Numeric_Kinds_Module
IMPLICIT NONE
INTEGER naph,ndel
LOGICAL little, medium, full
REAL(Kind=WP_Kind)ucheck(naph*naph), umat(ndel, naph)

INTEGER ithcll,ithsub,iaph
DATA ithcll/0/, ithsub/0/
DATA little/.false./, medium/.false./, full/.false./
CALL popt('delchk  ', little, medium, full, ithcll, ithsub)
!
WRITE(Out_Unit,*)'begining of routine delchk'
!-----------------------------------------------------------------
!     ucheck=umat(transpose)*umat
!-----------------------------------------------------------------
 CALL dgemm('T','N',naph,naph,ndel,1.0d0,umat,ndel,umat,ndel,0.0d0,ucheck,naph)
!-----------------------------------------------------------------
!     Print the diagonal elements of ucheck if little=.true.
!-----------------------------------------------------------------
IF(little)THEN
   WRITE(Out_Unit,*)'Diagonal elements of delchk unitary matrix'
   WRITE(Out_Unit,*)'IF the magnitude of the diagonal elements'
   WRITE(Out_Unit,*)' are almost one the Delves basis is OK'
   WRITE(Out_Unit,13)(ucheck(iaph), iaph=1, naph*naph, naph+1)
ENDIF
!-----------------------------------------------------------------
!     Print the entire ucheck matrix if full=.true.
!-----------------------------------------------------------------
IF(full)THEN
   WRITE(Out_Unit,*)'ucheck matrix from delchk'
   WRITE(Out_Unit,*)'IF the magnitude of the diagonal elements'
   WRITE(Out_Unit,*)' are almost 1 the Delves basis is OK'
   CALL MxOutL(ucheck, naph, naph, 0, 'aph  ', 'aph  ')
ENDIF
!
WRITE(Out_Unit,*)'END of routine delchk'
RETURN
 13   FORMAT(1x, 10e12.5)
!---------------------***END-delchk-------------------------------
END
