SUBROUTINE potvib (n, r, u, up, upp)
USE FileUnits_MODULE
USE Narran_Module
USE chltot_MODULE
USE Arrch_Module
USE region_MODULE
USE Integrat_Module
IMPLICIT NONE
INTEGER nn2, n, i, nstrt, nend
REAL(Kind=WP_Kind) r, ra, rb, rc
!#include <Inc.unitno>
!-----------------------------------------------------------------------
! calculates the matrix elements of the potential.
! written by g. a. parker
!-----------------------------------------------------------------------
!#include <narran.parms>
REAL(Kind=WP_Kind) u(*), up(*), upp(*)

!#include <Inc.chltot>
!#include <Inc.arrch>
!#include <Inc.region>
!#include <Inc.chint>

nn2=n*(n+1)/2
DO i=1,nn2
    u(i)=0.d0
    up(i)=0.d0
    upp(i)=0.d0
ENDDO

IF(integrat(1))THEN
    ra=r
    nstrt=1
    nend=nchanl(1)
    karran=1
    CALL potva (n,ra,u,up,upp,nstrt,nend)
ENDIF

IF(integrat(2))THEN
    rb=r
    nstrt=nchanl(1)+1
    nend=nchanl(1)+nchanl(2)
    karran=2
    CALL potvb (n,rb,u,up,upp,nstrt,nend)
ENDIF

IF(integrat(3))THEN
    rc=r
    nstrt=nchanl(1)+nchanl(2)+1
    nend=nchanl(1)+nchanl(2)+nchanl(3)
    karran=3
    CALL potvc (n,rc,u,up,upp,nstrt,nend)
ENDIF

IF(iregion=='delves ')CALL potth(n,r,u,up,upp)

RETURN
END
