SUBROUTINE sector1(rholst,rhonow,olast,okeep,rleft,rright,rmax,isect,nrho)
USE FileUnits_Module
USE region_Module
USE fuzzy_Module
USE rhovalue_Module
IMPLICIT NONE
INTEGER irhoen, ithcll, ithsub, isect, j, i, jj, nrho
REAL(Kind=WP_Kind) blast, okeep, olast, delrho, rholst, bkeep, rbegin
REAL(Kind=WP_Kind) deltamin, deltamax, rhonow, rr, rmax, rleft, rright
!#include <Inc.unitno>
INTEGER rhosave
LOGICAL little,medium,full

DATA irhoen/9999/
DATA ithcll/0/,ithsub/0/
DATA little/.false./,medium/.false./,full/.false./
CALL popt('sector  ',little,medium,full,ithcll,ithsub)
isect=isect+1
IF(iregion=='delves '.or.iregion=='djproj ')THEN
blast=0.0d0
bkeep=sqrt(2.0d0**1.5d0*log(1.0d0/okeep))
IF(olast/=0.0d0)THEN
    IF(olast>=0.9999) olast=.99999
    blast=sqrt(2.0d0**1.5d0*log(1.0d0/olast))
    delrho=delrho*rholst/(rholst-delrho)*bkeep/blast
    ELSE
    delrho=.05d0*rholst*bkeep
ENDIF
IF(delrho<deltamin*rholst)THEN
    delrho=deltamin*rholst
ELSEIF(delrho>deltamax*rholst)THEN
    delrho=deltamax*rholst
ENDIF
rhonow=rholst+delrho
rr=(rholst+rhonow)/2.d0
IF(rr<rmax.AND.ABS(rr-rmax)>fuzz) nrho=nrho+1
ELSE
!-----------------------------------------------------------------------
! Find the next sector in the aph region.
!
! We have just finished a sector, so read in the mid-points of
! the next two sectors
! rholst is the mid-point of the next
! sector and rhonow is the mid-point of the sector after that.
!-----------------------------------------------------------------------
    READ(Ovr_Unit,end=32) nrho,rholst,rhonow
    IF(little)WRITE(Out_Unit,10) nrho,rholst,rhonow
    rhosave=nrho
    nrho=irhoen
ENDIF
!-----------------------------------------------------------------------
! Compute the left and right boundaries of the sector.
! On entry, rleft and rright are the left and right boundaries of the
! sector that was just integrated.
!-----------------------------------------------------------------------
IF(isect/=1) rleft=rright
rright=(rholst+rhonow)/2.d0
!-----------------------------------------------------------------------
! If this is the first time into this routine we have to find the sector
! which brackets the starting distance, rbegin.
!-----------------------------------------------------------------------
IF(isect==1.AND.iregion=='aph    ')THEN
rbegin=rleft
rleft=2.0d0*rholst-rright
IF(little)WRITE(Out_Unit,'(A,3ES23.15)')'rbegin,rleft,rright=', rbegin,rleft,rright
IF(rbegin<rright)THEN
!-----------------------------------------------------------------------
! The beginning distance is in or BEFORE the first sector.
!-----------------------------------------------------------------------
    IF(rbegin>=rleft.or.ABS(rbegin-rleft)<fuzz)THEN
    rleft=rbegin
    ELSE
    WRITE(Out_Unit,*)'Integration starts at ',rleft, ', not rbegin=',rbegin
    STOP 'sector'
    ENDIF
ELSE
!-----------------------------------------------------------------------
! The starting distance is in some sector beyond the first one.
!-----------------------------------------------------------------------
    DO 21 j=1,1000
!-----------------------------------------------------------------------
! Read the rest of the overlap file for this sector.
!-----------------------------------------------------------------------
    DO 20 i=1,100

    READ(Ovr_Unit)nrho
    DO jj=1,nrho
        READ(Ovr_Unit)
    ENDDO
!-----------------------------------------------------------------------
! Read the next sector, find the boundaries, and see if this is the
! desired sector. Note that rbegin>=rleft.
!-----------------------------------------------------------------------
        READ(Ovr_Unit,end=31)nrho,rholst,rhonow
        IF(nrho==rhosave+1)THEN
            rhosave=nrho
            goto 25
        ENDIF
20         CONTINUE
    WRITE(Out_Unit,*)' Error: did not find the next sector'
    WRITE(Out_Unit,*)' Did you use more than 100 lambda in the'
    WRITE(Out_Unit,*)' abm part of the code?'
    STOP 'sector'

25         CONTINUE
    isect=isect+1
    rleft=rright
    rright=0.5*(rholst+rhonow)
    IF(full)THEN
        WRITE(Out_Unit,*)'nrho,rholst,rhonow=',nrho,rholst,rhonow
        WRITE(Out_Unit,*)'isect,rleft,rright=',isect,rleft,rright
    ENDIF
    IF(rbegin<rright)THEN
!             rleft=rbegin
        GOTO 23
    ENDIF
21       CONTINUE
    WRITE(Out_Unit,*)'error: did not find rbegin'
    STOP 'sector'
23       CONTINUE
ENDIF
nrho=irhoen
ENDIF

IF( iregion == 'aph    ' )THEN
!-----------------------------------------------------------------------
! If rmax<=rright, then this is the last sector to propagate and
! we set nrho=isect to tell hybrid not to propagate any farther.
!-----------------------------------------------------------------------
IF(full)WRITE(Out_Unit,'(A,2ES23.15,i5)')'rmax, rright, nrho=',rmax, rright,nrho
IF(full)WRITE(Out_Unit,'(A,2ES23.15)')'rmax-rright=',rmax-rright
IF(rmax<=rright)THEN
    nrho=isect
    rright=rmax
ENDIF
ELSE
IF( isect==nrho ) rright=rhonow
ENDIF
IF(rhonow>rmax.AND.ABS(rhonow-rmax)>fuzz) rhonow=rmax
IF(rright>rmax.AND.ABS(rright-rmax)>fuzz) rright=rmax
IF(little)WRITE(Out_Unit,'(A,ES23.15,A,ES23.15,A,I5,A,I6)')'rleft=',rleft,' rright=', rright,' isect=', isect,' nrho=',nrho
RETURN
!
31   WRITE(Out_Unit,*)'error: could not find rbegin ', 'in the overlap file'
STOP 'sector'
32   WRITE(Out_Unit,*)'error: propagated past end of overlap file'
WRITE(Out_Unit,10)nrho,rholst,rhonow
STOP 'sector'
!
10 FORMAT(//,'In Sector1: nrho=',i5,' rholst=',f10.5,' rhonow=',f10.5)
!-----------------------***end-sector***-------------------------------
END
