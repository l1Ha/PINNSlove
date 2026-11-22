SUBROUTINE Angle (nglegn, nhermt, xptg, xpth, ThetAPH, chitau, nang, &
                  iarran, sthaph, cthaph, wptg, wpth, weight)
!


!            P U R P O S E   O F    S U B R O U T I N E
!   This routine deterimines the APH theta and chi angles from the
!   quadrature points.


!            I N P U T    A R G U M E N T S
!   nglegn   an array containing the number of Gauss_Legendre quadrature
!            point for each arrangement channels.
!   nhermt   an array containing the number of Gauss_Hermite quadrature
!            points for each arrangement channel.
!   xptg     Gauss_Legendre quadrature points.
!   xpth     Gauss_Hermite quadrature points.
!   mxglegn  max no of Gauss_Legendre quad pts.
!   mxhermt  max no of Gauss_Hermite quad pts.

!            O U T P U T    A R G U M E N T S
!   ThetAPH  list of APH theta angles.
!   chitau   list of prin value APH chi angles meas from quad channel.
!   nang     total number of angles.
!   nangch   an array containing the starting position
!   iarran   arrangement to which quadrature point belongs.

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
USE Numeric_Kinds_Module
USE Parms_Module
USE TotJ_Module
USE FileUnits_Module
USE QCase_Module
USE MassFactor2_Module
USE AzBzCz_Module
USE Numbers_Module
IMPLICIT NONE
INTEGER karran, nang, iang, ngleg, narr,iglegn, ihermt
INTEGER nglegn(narran), nhermt(narran)
INTEGER iarran(nang)
REAL(Kind=WP_Kind) st, ct, s2td, c2td, sqr, thlitdel, zmc, tanth
REAL(Kind=WP_Kind) sin, cos, atan2, max, sign, fac
REAL(Kind=WP_Kind) xptg(mxglegn,narran), wptg(mxglegn,narran)
REAL(Kind=WP_Kind) xpth(mxhermt,narran), wpth(mxhermt,narran)
REAL(Kind=WP_Kind) sthaph(nang), cthaph(nang), a(3), x(2)
REAL(Kind=WP_Kind) ThetAPH(nang), chitau(nang)
REAL(Kind=WP_Kind) weight(1)
INTRINSIC sqrt, sin, cos, atan, atan2, max, sign

!-----------------------------------------------------------------------
! Determine printing options.
!-----------------------------------------------------------------------
LOGICAL little, medium, full
INTEGER ithcall, ithsub
DATA ithcall /0/, ithsub /0/
DATA little /.false./, medium /.false./, full /.false./
CALL popt ('angle   ', little, medium, full, ithcall, ithsub)

!  alternate quadrature case
IF(little) WRITE(Out_Unit,*) 'qcase=', qcase

IF(qcase)THEN
   CALL angle2(ngch, ngth, xch, xth, ThetAPH, chitau, nang, sthaph, cthaph, wch, wth, weight, iarran)
   RETURN
ENDIF

!
iang = 0
narr=narran
IF(symmetry) narr=2
DO karran = 1, narr
   a(1)= az(karran)
   a(3)=-bz(karran)
   fac = one
   ngleg = nglegn(karran)
   IF(symmetry.AND.karran==1)THEN
      fac = two
      ngleg = nglegn(1)/2
   ENDIF
   DO iglegn = 1, ngleg
!----------------------------------------------------------------------
!  cosine and sine of big theta of quadrature arrangement.
!----------------------------------------------------------------------
      ct=xptg(iglegn,karran)
      st=sqrt(one-ct**2)
      DO ihermt = 1, nhermt(karran)
         iang = iang + 1
         weight(iang) =fac*wptg(iglegn,karran)*wpth(ihermt,karran)
         weight(iang) = sqrt(weight(iang))
         iarran(iang)=karran
         zmc = xpth(ihermt, karran) - cz(karran)
         a(2)=-zmc
!----------------------------------------------------------------------
!  solve for delves little theta of quadrature arrangement.
!  solve quadratic for tanth in a way which is always accurate.
!  numerical recipes, p. 145
!----------------------------------------------------------------------
         zmc=-(a(2)+sign(one,a(2))*sqrt(a(2)**2-four*a(1)*a(3)))
         zmc=zmc/two
         x(1)=zmc/a(1)
         x(2)=a(3)/zmc
         tanth=max(x(1),x(2))
         thlitdel=atan(tanth)
         s2td = sin(two*thlitdel)
         c2td = cos(two*thlitdel)
         sqr  = sqrt(c2td**2+s2td**2*ct**2)
         sthaph(iang)=sqr
         cthaph(iang)=s2td*st
         ThetAPH(iang) = atan2(sthaph(iang),cthaph(iang))
         chitau(iang)  = atan2(s2td*ct/sqr,c2td/sqr)/two
!----------------------------------------------------------------------
!  this chitau is the principal value chi measured from the quadrature
!  arrangement; i.e., its range is -pi/2, pi/2
!----------------------------------------------------------------------
      ENDDO
   ENDDO
ENDDO
IF(iang/=nang)THEN
   WRITE(Out_Unit,*)'iang ne.nang: ',iang,nang
   STOP 'angle'
ENDIF

IF(little)THEN
   WRITE(Out_Unit,*)'routine angle:'
   WRITE(Out_Unit,*)'nglegn = ',nglegn,'  nhermt = ',nhermt
ENDIF
IF(full)THEN
   WRITE(Out_Unit,*)' iang,  ThetAPH,  chitau'
   DO iang = 1, nang
      WRITE(Out_Unit,*)iang,ThetAPH(iang),chitau(iang)
   ENDDO
ENDIF

RETURN
ENDSUBROUTINE Angle
