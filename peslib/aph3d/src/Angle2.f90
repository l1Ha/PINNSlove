SUBROUTINE Angle2(ngch, ngth, xch, xth, ThetAPH, chitau, nang, sthaph, cthaph, wch, wth, weight, iarran)
!

!            P U R P O S E   O F    S U B R O U T I N E
!   This routine deterimines the APH theta and chi angles from the
!   quadrature points.


!            I N P U T    A R G U M E N T S
!   ngch  the number of Gauss_Legendre quadrature in the chi direction.
!   ngth  the number of Gauss_Legendre quadrature in the
!         theta direction.
!   xch   Gauss_Legendre quadrature points in the chi direction.
!   xth   Gauss_Legendre quadrature points in the theta direction.

!            O U T P U T    A R G U M E N T S
!   ThetAPH  list of APH theta angles.
!   chitau   list of prin value APH chi angles meas from quad channel.
!   nang     total number of angles.
!   iarran   arrangement to which quadrature point belongs

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
USE Numeric_Kinds_Module
USE TotJ_Module
USE FileUnits_Module
USE Numbers_Module
IMPLICIT NONE
INTEGER ngch, ngth, nang, iang, ngleg, igch, igth, iarran(nang)
REAL(Kind=WP_Kind) sin, cos, max, sign, fac
REAL(Kind=WP_Kind) xch(nang), xth(nang), ThetAPH(nang), chitau(nang)
REAL(Kind=WP_Kind) sthaph(nang), cthaph(nang), wch(nang), wth(nang)
REAL(Kind=WP_Kind) weight(nang)
INTRINSIC sqrt, sin, cos, atan, atan2, max, sign

!-----------------------------------------------------------------------
! Determine printing options.
!-----------------------------------------------------------------------
LOGICAL little, medium, full
INTEGER ithcall, ithsub
DATA ithcall /0/, ithsub /0/
DATA little /.false./, medium /.false./, full /.false./
CALL popt ('angle   ', little, medium, full, ithcall, ithsub)

!
iang = 0
fac = one
ngleg = ngch
IF(symmetry)THEN
   fac = two
   ngleg = ngch/2
ENDIF
DO igch = 1, ngleg
   DO igth = 1, ngth
      iang = iang + 1
      iarran(iang)=1
      weight(iang) =fac*wch(igch)*wth(igth)
      weight(iang) = sqrt(weight(iang))
      ThetAPH(iang) = acos(xth(igth))/2.d0
      cthaph(iang)= cos(ThetAPH(iang))
      sthaph(iang)= sin(ThetAPH(iang))
      chitau(iang)  = xch(igch)
!----------------------------------------------------------------------
!  this chitau is the principal value chi measured from the quadrature
!  arrangement; i.e., its range is -pi/2, pi/2
!----------------------------------------------------------------------
   ENDDO
ENDDO


IF(little)THEN
   WRITE(Out_Unit,*)'routine angle2:'
   WRITE(Out_Unit,*)'ngch = ',ngch,'  ngth = ',ngth
   WRITE(Out_Unit,*)'symmetry=',symmetry
ENDIF

IF(iang/=nang)THEN
   WRITE(Out_Unit,*)'iang/=nang: ',iang,nang
   STOP 'angle2'
ENDIF


IF(full)THEN
   WRITE(Out_Unit,*)' iang,  ThetAPH,  chitau'
   DO iang = 1, nang
      WRITE(Out_Unit,*)iang,ThetAPH(iang),chitau(iang)
   ENDDO
ENDIF


RETURN
ENDSUBROUTINE Angle2
