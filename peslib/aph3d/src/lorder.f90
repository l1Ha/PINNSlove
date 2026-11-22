SUBROUTINE lorder(philist, lowel, nelem, mxnode, mxelmnt)
USE Numeric_Kinds_Module
USE FileUnits_Module
!
!     $RCSfile: lorder.f,v $ $Revision: 1.14 $
!     $Date: 89/10/18 14:18:06 $
!     $State: Stable $
!

!     P U R P O S E   O F    S U B R O U T I N E

!     this routine makes a list of element numbers in lowel such that lowel(
!     contains the element number with the highest value of the phirms


!     I N P U T    A R G U M E N T S
!
!     philist
!     lowel
!     nelem
!     mxnode
!     mxelmnt


IMPLICIT NONE

!     L O G I C A L S

LOGICAL little, medium, full

!     I N T E G E R S
INTEGER ithcall, ithsub

INTEGER lowel, nelem, mxnode, mxelmnt, iel, ibot, itop, jk,  islot
!     R E A L S
REAL(Kind=WP_Kind) philist, phim, phicur
!     D I M E N S I O N S
DIMENSION philist(mxelmnt), lowel(mxelmnt)


DATA ithcall/0/, ithsub/0/
DATA little/.false./, medium/.false./, full/.false./
CALL popt ('lorder  ', little, medium, full, ithcall, ithsub)

lowel(1) = 1
DO 155 iel=1,nelem
   phim = philist(iel)
   IF(full) WRITE(Out_unit,*)'phim=',phim,' iel=',iel

!     philist is indexed by element number and contains the largest value of
!     phirms found at its node points.  lowel is indexed numerically with
!     islot = 1 being the element number which has the largest value of the
!     phirms.
!     test for a phirms either below/above the lowest/highest previous
!     value of the phirms for all elements; take care of these cases separat
   IF(philist(lowel(1))<=phim)THEN
      ibot=0
      IF(medium) WRITE(Out_unit,*)'philist(lowel(1)=',philist(lowel(1))
      GOTO 158
   ELSEIF(philist(lowel(iel-1))>=phim)THEN
      ibot = iel-1
      IF(medium) WRITE(Out_unit,*)'philist(lowel(iel-1)=', philist(lowel(iel-1))
      GOTO 158
   ENDIF
!     find where this value of the phirms fits in the high to low list
!     of elements; we use a halving procedure to try and find where the elem
!     belongs as fast as possible
   itop = iel - 1
   ibot = 1
   DO 157 jk=1,mxelmnt
      islot = (itop-ibot)/2 + ibot
      phicur = philist(lowel(islot))
      IF(medium) WRITE(Out_unit,*)'phicur=',phicur
      IF(phicur<phim)THEN
         itop = islot
      ELSEIF(phicur>phim)THEN
         ibot = islot
      ELSE
         ibot = islot
         itop = islot
      ENDIF
      IF(itop-ibot<=1)GOTO 158
 157     CONTINUE
   WRITE(Out_unit,*)'***error***: could not find slot for iel=',iel
   WRITE(Out_unit,*)'phim=',phim
   STOP 'lorder'
 158     CONTINUE
!     insert a new element into the low/high list; this requires moving all
!     elements with a phirms value lower than the new element up on the list
   DO 159 jk=1,iel-ibot-1
      lowel(iel-jk+1) = lowel(iel-jk)
 159     CONTINUE
   lowel(ibot+1) = iel
   IF(full) WRITE(Out_unit,*)'lowel=',(lowel(jk), jk=1,iel)
 155  CONTINUE
IF(full) WRITE(Out_unit,*)'exiting lorder'
RETURN
END
