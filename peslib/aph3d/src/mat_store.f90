SUBROUTINE mat_store(wmat, nwmat, wsub, nsub1, nsub2, lam, lamp)
USE FileUnits_MODULE
USE pbasis_Module
IMPLICIT NONE
INTEGER lam, start, preadd, postadd, i, ij, isub, j, lamp, nwmat, nsub1, nsub2, megamax
!     >        lammax, megamin, megacoup
REAL(Kind=WP_Kind) wmat(nwmat*(nwmat+1)/2), wsub(nsub1*nsub2)
!      COMMON / basis / nlam(0:1000), lammin, lammax, megamin, megamax, megacoup
!#include <Inc.unitno>
!-----------------------------------------------------------------------
! revised by X.L.
megamax=maxmega
IF(nsub1>nwmat.or.nsub2>nwmat)THEN
    WRITE(Msg_Unit,*)'Error nsub1>nwmat.or.nsub2>nwmat'
    WRITE(Msg_Unit,*)'nsub1,nsub2,nwmat=',nsub1,nsub2,nwmat
    STOP 'mat_store'
ENDIF

! Determine last element of the previous column.
!-----------------------------------------------------------------------
start=0
DO 2 i=lammin,lamp-1
    start = start+nlam(i)
2     CONTINUE
start = start*(start+1)/2
!
preadd=0
DO 4 i=lammin,lam-1
    preadd = preadd+nlam(i)
4     CONTINUE
!
ij = start
isub = 0
IF(lamp==lam)THEN
    DO 8 i = 1, nsub1
    ij = preadd+ij
    DO 10 j = 1, i
        isub = isub +1
        ij = ij +1
        wmat(ij) = wsub(isub)
10          CONTINUE
8        CONTINUE
ELSE
    postadd=0
    DO 6 i=lam+1,lamp-1
    postadd = postadd+nlam(i)
6        CONTINUE

    DO 20 i = 1, nsub1
    ij = preadd+ij
    DO 25 j = 1, nsub2
        isub = isub +1
        ij = ij +1
        wmat(ij) = wsub(isub)
25          CONTINUE
    ij = postadd+ij+i
20        CONTINUE
ENDIF

RETURN
END
