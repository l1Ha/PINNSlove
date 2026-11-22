SUBROUTINE DetQuad
USE Numeric_Kinds_Module
USE Parms_Module
USE FileUnits_Module
USE QDiv_Module
USE QCase_Module
IMPLICIT NONE
LOGICAL debug
INTEGER ithdiv, ichdiv, iang
REAL(Kind=WP_Kind) cos, xhigh, xlow

debug = .false.

WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Theta integration'
ngth = 0
DO ithdiv = 1, nthdiv
   xlow = cos(2.d0*thdiv(ithdiv))
   xhigh = cos(2.d0*thdiv(ithdiv+1))
   CALL glegen (nqth(ithdiv), xth(ngth+1), wth(ngth+1), xlow, xhigh)
   ngth = ngth + nqth(ithdiv)
   WRITE(Out_Unit,*)'ithdiv=',ithdiv,nqth(ithdiv),ngth
   WRITE(Out_Unit,'(4ES20.12)')thdiv(ithdiv),thdiv(ithdiv+1),xlow,xhigh
ENDDO

DO iang = 1, ngth
   wth(iang) = -wth(iang)
ENDDO

IF(debug)THEN
   WRITE(Out_Unit,*)'Theta quadrature points'
   WRITE(Out_Unit,*)(xth(iang),iang=1,ngth)
   WRITE(Out_Unit,*)'Theta weights'
   WRITE(Out_Unit,*)(wth(iang),iang=1,ngth)
ENDIF

WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Chi integration'
ngch = 0
DO ichdiv = 1, nchdiv
   xlow = chdiv(ichdiv)
   xhigh = chdiv(ichdiv+1)
   CALL glegen (nqch(ichdiv), xch(ngch+1), wch(ngch+1), xlow, xhigh)
   ngch = ngch + nqch(ichdiv)
   WRITE(Out_Unit,*)'ichdiv=',ichdiv,nqch(ichdiv),ngch
   WRITE(Out_Unit,'(4ES20.12)')chdiv(ichdiv),chdiv(ichdiv+1),xlow,xhigh
ENDDO

IF(debug)THEN
   WRITE(Out_Unit,*)'Chi quadrature points'
   WRITE(Out_Unit,*)(xch(iang),iang=1,ngch)
   WRITE(Out_Unit,*)'Chi weights'
   WRITE(Out_Unit,*)(wch(iang),iang=1,ngch)
ENDIF

RETURN
ENDSUBROUTINE detquad
