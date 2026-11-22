SUBROUTINE aphchk(ndel, ucheck, naph, umat)
!
!-----------------------------------------------------------------
!     This routine provides a unitarity check on the APH to Delves
!     transformation matrix. If the magnitude of the diagonal elements
!     are almost 1 this indicates that the Delves basis nearly
!     complete (sufficiently large enough to expand the APH functions).
!-----------------------------------------------------------------
!this routine is called by: AphDel_Old or AphDel_New
!this routine calls: popt, mxoutl
!-----------------------------------------------------------------------
USE FileUnits_Module
USE crays_Module
USE Numeric_Kinds_Module
IMPLICIT NONE
INTEGER NDel, NAPH
INTEGER ithcll,ithsub,idel
LOGICAL little, medium, full
REAL(Kind=WP_Kind)ucheck(ndel*ndel), umat(ndel,naph)


DATA ithcll/0/,ithsub/0/
DATA little/.false./,medium/.false./,full/.false./
CALL popt('aphchk  ', little, medium, full, ithcll, ithsub)
!
WRITE(Out_Unit,*)'begining of routine aphchk'
!-----------------------------------------------------------------
!     ucheck=umat*umat(transpose)
!-----------------------------------------------------------------
 CALL dgemm('N','T',ndel,ndel,naph,1.0d0,umat,ndel,umat,ndel,0.d0,ucheck,ndel)
!-----------------------------------------------------------------
!     Print the diagonal elements of ucheck if little=.true.
!-----------------------------------------------------------------
IF(little)THEN
   WRITE(Out_Unit,*)'Diagonal elements of aphchk unitary matrix'
   WRITE(Out_Unit,*)'IF the magnitude of the diagonal elements'
   WRITE(Out_Unit,*)' are almost one the APH basis is OK'
   WRITE(Out_Unit,13)(ucheck(idel), idel=1, ndel*ndel, ndel+1)
ENDIF
!-----------------------------------------------------------------
!     Print the entire ucheck matrix if full=.true.
!-----------------------------------------------------------------
IF(full)THEN
   WRITE(Out_Unit,*)'ucheck matrix from aphchk'
   WRITE(Out_Unit,*)'IF the magnitude of the diagonal elements'
   WRITE(Out_Unit,*)' are almost 1 the Delves basis is OK'
   CALL MxOutL(ucheck, ndel, ndel, 0, 'space', 'space')
ENDIF
!
  WRITE(Out_Unit,*)'END of routine aphchk'
RETURN
 13   FORMAT(1x, 10e12.5)
!---------------------***end-aphchk-------------------------------
ENDSUBROUTINE aphchk
