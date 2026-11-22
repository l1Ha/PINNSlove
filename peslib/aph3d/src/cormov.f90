SUBROUTINE cormov(nelem, idiv, lcor, iwork1, linkel, iwork2, sinlist, philist, work1, work2)
USE Numeric_Kinds_Module
USE FileUnits_Module
!
!     $RCSfile: cormov.f,v $ $Revision: 1.14 $
!     $Date: 89/10/18 14:17:57 $
!     $State: Stable $
!


!     P U R P O S E   O F    S U B R O U T I N E
!-----------------------------------------------------------------------
!     Rearrange the element's corner and divison arrays so that they corresp
!     with the new element numbers due to having new elements.  Note that th
!     effectively renumbers the elements since idiv and lcor are the only th
!     that actually defines the elements. This routine
!     should be called after all the new elements have been made, not just 4
!     ones from one old one.
!     The arrarys idiv, lcor, iwork1, and iwork2 are changed on RETURN.
!-----------------------------------------------------------------------


IMPLICIT NONE

!     L O G I C A L S

LOGICAL idbug

!     I N T E G E R S

INTEGER idiv, iel2, iwork1, iwork2, last, lcor, linkel, m, mm, nelem

!     R E A L S
REAL(Kind=WP_Kind) philist, sinlist, work1, work2

!     D I M E N S I O N S

DIMENSION idiv(nelem),    lcor(2,nelem),   iwork1(nelem), linkel(nelem),  iwork2(2,nelem), sinlist(nelem), &
           philist(nelem), work1(nelem),    work2(nelem)

DATA idbug/.false./
IF(idbug) WRITE(Out_unit,*)'entering cormov'
!
iwork1(1)=idiv(1)
iwork2(1,1)=lcor(1,1)
iwork2(2,1)=lcor(2,1)
work1(1)=philist(1)
work2(1)=sinlist(1)
last=1
!
DO m=2,nelem
   iel2=linkel(last)
!-----------------------------------------------------------------------
!     Make sure that iel2 not greater than nelem.
!-----------------------------------------------------------------------
   IF(iel2>nelem)THEN
      WRITE(Out_unit,*)'error iel2>nelem: iel2,nelem = ',iel2,nelem
      STOP 'cormov'
   ENDIF
   iwork1(m)=idiv(iel2)
   work1(m)=philist(iel2)
   work2(m)=sinlist(iel2)
   DO mm=1,2
      iwork2(mm,m)=lcor(mm,iel2)
   ENDDO
   last=iel2
ENDDO
DO m=1,nelem
   idiv(m)=iwork1(m)
   philist(m)=work1(m)
   sinlist(m)=work2(m)
   lcor(1,m)=iwork2(1,m)
   lcor(2,m)=iwork2(2,m)
ENDDO
IF(idbug)THEN
   WRITE(Out_unit,*)(idiv(m),m=1,nelem)
   WRITE(Out_unit,*)(philist(m),m=1,nelem)
   WRITE(Out_unit,*)'exiting cormov'
ENDIF
RETURN
END
