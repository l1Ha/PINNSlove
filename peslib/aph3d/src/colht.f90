SUBROUTINE colht (mht, nd, lm, neq, nid, nidm, idi)
Use FileUnits_Module
!
!   $RCSfile: colht.f,v $ $Revision: 1.15 $
!   $Date: 89/11/16 11:02:54 $
!   $State: Stable $
!

!            P U R P O S E   O F    S U B R O U T I N E

!-----------------------------------------------------------------------
! Program
!     to calculate the column heights.
!
!           x x 0 0 0 0 0 0            1  3
!             x x 0 x 0 0 0               2  5    10
!               x 0 x 0 x 0                  4     9    18
!                 x x x x 0                     6  8 13 17
!                   x x x 0                        7 12 16
!                     x x x                          11 15 21
!                       x x                             14 20
!                         x                                19
!
!         Matrix with skyline             Storage Scheme
!
!   The matrix is assumed to be symmetric and to have zero elements
!   above the skyline.
!-----------------------------------------------------------------------


!            I N P U T    A R G U M E N T S

!     mht       column heights.
!     nd
!     lm
!     neq
!     nid(NDISCE) Array with the number of terms in the constraint equations
!     nidm
!     idi

IMPLICIT NONE
INTEGER nd, nidm, nce, ls, ii, j, me, i, neq,  jj
INTEGER lm(1), mht(1), nid(1), idi(nidm,1)
ls=neq+1
DO i=1,18
   IF(lm(i)<0)THEN
      nce=-lm(i)
      DO j=1,nid(nce)
         jj=idi(j,nce)
         IF(jj-ls<0)ls=jj
      ENDDO
   ELSEIF(lm(i)>0)THEN
      IF(lm(i)-ls<0)ls=lm(i)
   ENDIF
ENDDO
DO i=1,18
   ii=lm(i)
   IF(ii<0)THEN
      nce=-ii
      DO j=1,nid(nce)
         ii=idi(j,nce)
         me=ii-ls
         IF(me>mht(ii)) mht(ii)=me
      ENDDO
   ELSEIF(ii>0)THEN
      me=ii-ls
      IF(me>mht(ii)) mht(ii)=me
   ENDIF
ENDDO

RETURN
END
