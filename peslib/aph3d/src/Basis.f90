SUBROUTINE Basis (xptg, wptg, xpth, wpth, dpndx, doscdz, ThetAPH, sthaph, cthaph, chiaph, pn, oscil,  &
                  nang, nangch, parity, parityfn, parityf2, bfaci, dzdthd, dlnb, dxdchi, dxdtha,      &
                  dthddtha, dthddchi, rho, mega, megamin, ioption, weight, nfac, fnorm, chitau, iarran)
!

!            P U R P O S E   O F    S U B R O U T I N E
!   This routine calculates Legendre and harmonic oscillator basis
!   functions.

!            I N P U T    A R G U M E N T S
!  nglegn   number of Gauss_Legendre quadrature points for each
!           arrangement channel.
!  nlegndre maximum order of the Legendre polynomial for each
!           arrangement channel.
!  nhermt   number of Gauss_Hermite quadrature points for each
!           arrangement channel.
!  noscil   maximum order of the oscillator basis for each
!           arrangement.
!  parity   parity of the surface-functions.
!  ioption  IF 1, calculates basis functions for use at present angles.
!           IF 2, calcs oscil and bfaci of previous rho at pres. angles.

!            O U T P U T    A R G U M E N T S
!  xptg     Gauss_Legendre quadrature points.
!  wptg     Gauss_Legendre weights.
!  xpth     Gauss_Hermite quadrature points.
!  wpth     Gauss_Hermite weights.
!  ThetAPH  APH theta angles.
!  chiaph   APH chi angles.
!  sthaph   array of sin(ThetAPH).
!  cthaph   array of cos(ThetAPH).
!  pn       Legendre polynomial basis.
!  oscil    Harmonic oscillator basis.
!  nang     number of APH theta and chi angles.
!  nangch   number of APH theta and chi angles for each arrangement
!           channel.
!  parityfn parity function. This function is 1 for even parity and
!           cos(chi) for odd parity.
!  weight   quadrature weights.

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
USE Numeric_Kinds_Module
USE Parms_Module
USE FileUnits_Module
USE PES_MODULE
USE Quantb_Module
USE Masses_Module
USE QCase_Module
USE Numbers_Module
USE Gaussb_Module
USE MassFactor2_Module
USE AzBzCz_Module
IMPLICIT NONE

!            I N T E G E R S
INTEGER karran, iglegn, ihermt, iarran(mxang)
INTEGER nang, nangch(narran+1), parity, nptchnl(narran)
INTEGER iang, n, l, mega, megamin, ioption

!            R E A L S
REAL(Kind=WP_Kind) denom, schi, cchi, thlitdel
REAL(Kind=WP_Kind) rho, tanth
REAL(Kind=WP_Kind) chi, xbig, sinths, cosths, bfac, zlit
REAL(Kind=WP_Kind) sinthd, costhd, azsbzc, denom3, btemp, y
REAL(Kind=WP_Kind) top, product, cos4
REAL(Kind=WP_Kind) beta, eps, sin4, cos, sin, asin

SAVE   denom, schi, cchi, thlitdel, thliteq
SAVE   tanth, chi, xbig, sinths, cosths, bfac, zlit
SAVE   sinthd, costhd, azsbzc, denom3
SAVE   se, ao, bo, co, btemp, y, top, product, cos4
SAVE   beta, eps, sin4, rhoprm, rhofix, cparm

!            D I M E N S I O N S
REAL(Kind=WP_Kind) xptg(mxglegn, narran), wptg(mxglegn, narran)
REAL(Kind=WP_Kind) xpth(mxhermt, narran), wpth(mxhermt, narran)
REAL(Kind=WP_Kind) pn(0:mxl, mxang, narran),oscil(0:mxn, mxang, narran)
REAL(Kind=WP_Kind) ThetAPH(mxang), chiaph(mxang),parityfn(mxang,narran)
REAL(Kind=WP_Kind) thliteq(narran), nfac(0:2*mxl), fnorm(0:mxn)
REAL(Kind=WP_Kind) sthaph(nang), doscdz(0:mxn, mxang, narran)
REAL(Kind=WP_Kind) dpndx(0:mxl,mxang,narran), parityf2(mxang,narran)
REAL(Kind=WP_Kind) bfaci(mxang,narran), dzdthd(mxang,narran)
REAL(Kind=WP_Kind) dlnb(mxang,narran), dxdchi(mxang,narran)
REAL(Kind=WP_Kind) dxdtha(mxang,narran), dthddtha(mxang,narran)
REAL(Kind=WP_Kind) dthddchi(mxang,narran), weight(mxang)
REAL(Kind=WP_Kind) se(narran), cthaph(mxang), cparm(narran)
REAL(Kind=WP_Kind) ao(narran), bo(narran), co(narran)
REAL(Kind=WP_Kind) rhoprm(narran), rhofix(narran), chitau(nang)

!            I N T R I N S I C    F U N C T I O N S
INTRINSIC cos, atan, sin, sqrt, tan, asin

!            E X T E R N A L S
EXTERNAL popt, glegen, hermit, angle, legendrd, harmonid

!-----------------------------------------------------------------------
! Determine printing options.
!-----------------------------------------------------------------------
LOGICAL little, medium, full, expfit, tmud2, lifhb, lcase
INTEGER ithcall, ithsub, ithcal2
DATA ithcall /0/, ithcal2 /0/, ithsub /0/, lcase/.True./
DATA little /.false./, medium /.false./, full /.false./, lcase/.false./
CALL popt ('basis   ', little, medium, full, ithcall, ithsub)

IF(Medium)WRITE(Out_Unit,*)'old_way=',old_way
WRITE(Msg_Unit,*)'old_way=',old_way

IF(old_way)THEN
  IF(PES_Name=='tmud2 pairwise add. ')THEN
     tmud2 = .true.
     IF(Medium)WRITE(Out_Unit,*)'tmud2=.true.'
  ELSE
     IF(Medium)WRITE(Out_Unit,*)'tmud2=.false.'
     tmud2 = .false.
  ENDIF

  IF(PES_Name=='LiFH_B BondOrder PES')THEN
     expfit = .true.
     lifhb = .true.
     WRITE(Out_Unit,*)'Using exponential fit to cparm'
  ELSE
     expfit = .false.
     lifhb = .false.
     IF(Medium)WRITE(Out_Unit,*)'Using morse fit to cparm'
  ENDIF
ENDIF

!-----------------------------------------------------------------------
! Obtain the quadrature points for the Gauss_Legendre and Gauss_Hermite
! integrations.  DO only on first CALL at first rho.
!-----------------------------------------------------------------------
IF(ioption==1)THEN
      IF(Medium)WRITE(Out_Unit,*)'cases=',lcase,qcase,ithcal2
!----------------------------------------------------------
!modified by X.L if change to qcase=F then we need to be sure ao,bo,co
!are correctly put. 1-10-2006
!---------------------------------------------
      IF(.NOT.qcase.AND.lcase)THEN
      DO karran=1,narran
         ao(karran)=az(karran)
         bo(karran)=bz(karran)
         co(karran)=cz(karran)
      ENDDO
         WRITE(Out_Unit,*)"ao=",ao
         WRITE(Out_Unit,*)"bo=",bo
         WRITE(Out_Unit,*)"co=",co
      ENDIF
!--------------------------------------------------------------------


      IF(ithcal2==1.AND.lcase.neqv.qcase)ithcal2=0
      lcase=qcase
   IF(ithcal2==0)THEN
      ithcal2=1
      IF(qcase)THEN
         CALL detquad
      ELSE
         DO karran = 1, narran
            CALL glegen (nglegn(karran), xptg(1,karran), wptg(1,karran), -one, one)
            CALL hermit (nhermt(karran), xpth(1,karran), wpth(1,karran), 1)
         ENDDO
      ENDIF

      DO karran = 1, narran
            se(karran)=rx(karran)*re(karran)/d(karran)
            IF(old_way)THEN
               IF(tmud2)THEN !  the following is used for tmu+d2
                  rhofix(karran)=1.025d0*se(karran)
               ELSEIF(lifhb)THEN !  the following is used for LiFH
                  rhofix(karran)=1.1*se(karran)
               ENDIF
            ENDIF
         ENDDO

      IF(old_way)THEN
      WRITE(Out_Unit,*)' rhofix=',(rhofix(karran),karran=1,narran)
      ENDIF
!----------------------------------------------------------------------
!  on first CALL at subsequent rho, save parameters of previous rho.
!-----------------------------------------------------------------------
   ELSEIF(mega==megamin)THEN
      DO karran=1,narran
         ao(karran)=az(karran)
         bo(karran)=bz(karran)
         co(karran)=cz(karran)
      ENDDO
   ENDIF

!  ---------------------------------------------------------------------
!  determine equilibrium value of thlitdel.
!  detn parameters relating hermite quadrature variable z to thlitdel.
!  z=a*tan(theta)-b/tan(theta)+c
!  --------------------------------------------------------------------
   DO karran=1,narran
      IF(old_way)THEN
         rhoprm(karran)=rho
      ENDIF
      cparm(karran)=calpha(karran)
      IF(old_way)THEN
         IF(rho<rhofix(karran))THEN
            IF(tmud2)THEN !  the following parameters are for tmu+d2 only
               rhoprm(karran)=rhofix(karran)
               cparm(1)=6.825d0 -2.4533d0*rho
            ELSEIF(lifhb)THEN
               rhoprm(karran)=rhofix(karran)
               y=exp(-balpha(3)*rhofix(3))
               cparm(3)=calpha(3)+dalpha(3)*y
            ELSE !  the following is exp fit to cparm. 
               y=exp(-balpha(karran)*rho)
               cparm(karran)=calpha(karran)+dalpha(karran)*y
            ENDIF      
         ELSE !  they will need to be changed IF we use a basis for another system in which rho<rhofix.
            IF(expfit)THEN  !  the following is exp fit to cparm.
               y=exp(-balpha(karran)*rho)
               cparm(karran)=calpha(karran)+dalpha(karran)*y
            ELSE !  the following is morse fit to cparm.
               y=exp(-balpha(karran)*(rhoprm(karran)-ralpha(karran)))
               cparm(karran)=calpha(karran)+dalpha(karran)*y*(y-two)
            ENDIF
         ENDIF
      ENDIF

      IF(.NOT.old_way)THEN
      CALL equilth (karran,rho,weau(karran),wexeau(karran),usys, zeta(karran),se(karran), thliteq(karran), &
                          delta(karran),beta, eps)

      beta=sqrt(beta)
      eps=eps/(six*beta)
      ELSE
      IF(se(karran)>rhoprm(karran))THEN
         WRITE(Out_Unit,*)'Error: Basis'
         WRITE(Out_Unit,*)'se(karran) must be < rhoprm(karran)'
         WRITE(Out_Unit,*)'se=', se(karran)
         WRITE(Out_Unit,*)'rhoprm', rhoprm(karran)
         
         WRITE(Msg_Unit,*)'Error: Basis'
         WRITE(Msg_Unit,*)'se(karran) must be < rhoprm(karran)'
         WRITE(Msg_Unit,*)'se=', se(karran)
         WRITE(Msg_Unit,*)'rhoprm', rhoprm(karran)
         !STOP 'Basis73' !tmpmod gregparker
      ENDIF

      thliteq(karran)=asin(se(karran)/rhoprm(karran))
      ENDIF


      tanth=tan(thliteq(karran))
      costhd=cos(thliteq(karran))
      sinthd=sin(thliteq(karran))
      sinths=sinthd**2
      sin4=sinths**2
      cos4=costhd**4
      IF(old_way)THEN
         beta=rhoprm(karran)*sqrt(usys*(weau(karran)-wexeau(karran)))*cparm(karran)
         eps=beta*anharm(karran)*sqrt(two*usys*wexeau(karran))*rhoprm(karran)/six
      ELSE
         beta=cparm(karran)*beta
         eps=cparm(karran)*anharm(karran)*eps
      ENDIF
      IF(old_way)THEN
         az(karran)=cos4*costhd*(beta-eps*sinthd)
         bz(karran)=sin4*costhd*(beta+eps*costhd/tanth)
      ELSE
!  parameters of z
         az(karran)=cos4*(beta+eps*tanth)
         bz(karran)=sin4*(beta-eps/tanth)
!  keep a positive.  IF necessary, only fit 2nd deriv, not third.
      IF(az(karran)<1.d-8)THEN
         az(karran)=1.d-8
      ENDIF
      ENDIF
      cz(karran)=bz(karran)/tanth - az(karran)*tanth
      IF(.NOT.old_way)THEN
      IF(az(karran)<=0.d0.or.bz(karran)<=0.d0)THEN
         WRITE(Out_Unit,*)'arrangement channel = ',karran
         WRITE(Out_Unit,*)'beta,eps,tan=',beta, eps, tanth
         WRITE(Out_Unit,*)'a,b,c= ', az(karran),bz(karran),cz(karran)
         STOP 'basis74'
      ENDIF
      ENDIF
   ENDDO

!-----------------------------------------------------------------------
! Determine the APH theta and chi angles.
!  This is the principal value chi measured from quadrature channel.
!-----------------------------------------------------------------------
   CALL angle (nglegn, nhermt, xptg, xpth, ThetAPH, chitau, nang,  &
                iarran, sthaph, cthaph, wptg, wpth, weight)

!-----------------------------------------------------------------------
! Determine functions of the Delves angles.
!  xbig is the cos(thbigdel) of the function arrangement.
!  zlit is dimensionless function of thlitdel of fcn arrangement.
!  karran runs over the function channels.
!  iarran(iang) is the quadrature channel.
!  chitau is  prin value chi of quadrature channel.
!  chi is chi of the function channel.
!  chiaph is usual chi measured from channel 1.
!  bfac is a factor that goes with the Jacobian and sho functions.
! Calculate the normalized Legendre polynomial basis and derivative.
! Calculate the normalized harmonic oscillator basis and derivative.
! Calculate the parity function.
!-----------------------------------------------------------------------
DO karran = 1, narran
   nptchnl(karran)=0
   DO iang = 1, nang
      chi=chitau(iang)-chij(karran,iarran(iang))
      IF(karran==1)chiaph(iang)=chi
      schi   = sin(two*chi)
      cchi   = cos(two*chi)
      top=sqrt(cthaph(iang)**2*cchi**2 + schi**2)
      denom=one/top
      xbig =  sthaph(iang)*schi*denom
      CALL legendrd (nlegndre(karran), mega, pn(0, iang, karran), dpndx(0,iang,karran), xbig, nfac)
      product=sthaph(iang)*cchi

      IF(product<zero)THEN
         tanth = (one-product)*denom
      ELSE
         tanth = top/(one+product)
      ENDIF
      thlitdel=atan(tanth)
      sinthd=sin(thlitdel)
      costhd=cos(thlitdel)
      sinths=sinthd**2
      cosths=costhd**2
      azsbzc=one/(az(karran)*sinths+bz(karran)*cosths)
      bfac=four*sinths*cosths*sqrt(azsbzc)
      IF(.NOT.qcase)THEN
         IF(karran==iarran(iang))THEN
            weight(iang)=weight(iang)*bfac
         ENDIF
      ELSE
         IF(karran==iarran(iang))THEN
            weight(iang)=weight(iang)
         ENDIF
      ENDIF
      zlit=az(karran)*tanth-bz(karran)/tanth+cz(karran)
      CALL harmonid (noscil(karran), oscil(0, iang, karran), doscdz(0,iang,karran), zlit, fnorm)
      IF(ABS(zlit)<4)THEN
         nptchnl(karran) = nptchnl(karran) + 1
      ENDIF
      parityfn(iang,karran)=cos(chi)**parity
!---------------------------------------------------------------------
!  now construct extra functions needed for kinetic energy.
!---------------------------------------------------------------------
      parityf2(iang,karran)=parity*sin(chi)
      bfaci(iang,karran)=one/bfac
      dzdthd(iang,karran)=sixteen*sinths*cosths*bfaci(iang, karran)**2
      dlnb(iang,karran)=four/tan(two*thlitdel)-(az(karran)-bz(karran))*sinthd*costhd*azsbzc
      denom3=denom**3
      dxdchi(iang,karran)=two*sthaph(iang)*cchi*denom3*cthaph(iang)**2
      dxdtha(iang,karran)=cthaph(iang)*schi*denom3
      dthddtha(iang,karran)=-half*cthaph(iang)*cchi*denom
      dthddchi(iang,karran)=xbig
   ENDDO
ENDDO
!  --------------------------------------------------------------------
!  when ioption=2, rerun the  parts of the calc necessary to get
!  oscil and bfaci for the previous rho.  This allows calc of basis of
!  previous rho at quadrature points of present rho.
!  --------------------------------------------------------------------
ELSE

DO karran = 1, narran
   DO iang = 1, nang
      chi=chitau(iang)-chij(karran,iarran(iang))
      schi   = sin(two*chi)
      cchi   = cos(two*chi)
      top=sqrt(cthaph(iang)**2*cchi**2 + schi**2)
      denom=one/top
      product=sthaph(iang)*cchi

      IF(product<zero)THEN
         tanth = (one-product)*denom
      ELSE
         tanth = top/(one+product)
      ENDIF
      thlitdel=atan(tanth)
      sinthd=sin(thlitdel)
      costhd=cos(thlitdel)
      sinths=sinthd**2
      cosths=costhd**2
      azsbzc=one/(ao(karran)*sinths+bo(karran)*cosths)    !TmpMod is this correct?
      !azsbzc=one/(az(karran)*sinths+bz(karran)*cosths)
      btemp=four*sinths*cosths*sqrt(azsbzc)
      zlit=ao(karran)*tanth-bo(karran)/tanth+co(karran)     !TmpMod is this correct?
      !zlit=az(karran)*tanth-bz(karran)/tanth+cz(karran)
      CALL harmonid (noscil(karran), oscil(0, iang, karran), doscdz(0,iang,karran), zlit, fnorm)
      bfaci(iang,karran)=one/btemp
   ENDDO
ENDDO

ENDIF
!-----------------------------------------------------------------------
! Print basis IF desired.
!-----------------------------------------------------------------------
IF(little)THEN
   WRITE(Out_Unit,*)'routine basis'
   WRITE(Out_Unit,*)'narran = ',narran
   WRITE(Out_Unit,*)'nglegn = ',nglegn
   WRITE(Out_Unit,*)'nlegndre = ',nlegndre
   WRITE(Out_Unit,*)'nhermt = ',nhermt
   WRITE(Out_Unit,*)'noscil = ',noscil
   DO karran=1,narran
      WRITE(Out_Unit,*)'arrangement channel = ',karran
      WRITE(Out_Unit,*)'a,b,c= ',az(karran),bz(karran),cz(karran)
   ENDDO
ENDIF
IF(medium)THEN
   DO karran = 1, 3
      WRITE(Out_Unit,*)'arrangement channel = ',karran
      WRITE(Out_Unit,*)'xptg = ',(xptg(iglegn,karran), iglegn=1,nglegn(karran))
      WRITE(Out_Unit,*)'wptg = ',(wptg(iglegn,karran), iglegn=1,nglegn(karran))
      WRITE(Out_Unit,*)'xpth = ',(xpth(ihermt,karran), ihermt=1,nhermt(karran))
      WRITE(Out_Unit,*)'wpth = ',(wpth(ihermt,karran), ihermt=1,nhermt(karran))

   ENDDO
      WRITE(Out_Unit,*)'nptchnl',nptchnl

ENDIF
IF(full)THEN
   DO karran=1, narran
      WRITE(Out_Unit,*)'arrangement channel = ',karran
      WRITE(Out_Unit,*)'a,b,c= ', az(karran),bz(karran),cz(karran)
      WRITE(Out_Unit,*)'pn = ',((pn(l, iang, karran), iang=1,nang), l = jmin(0,karran,mega), nlegndre(karran),jskip(karran))
      WRITE(Out_Unit,*)'oscil = ',((oscil(n, iang, karran),iang=1,nang), n = 0, noscil(karran))
      WRITE(Out_Unit,*)'parityfn= ',(parityfn(iang,karran), iang=1,nang)
      WRITE(Out_Unit,*)'dpndx = ',((dpndx(l, iang, karran),  iang=1,nang), l = jmin(0,karran,mega), nlegndre(karran), jskip(karran))
      WRITE(Out_Unit,*)'doscdz = ',((doscdz(n,iang,karran),iang=1,nang), n = 0, noscil(karran))
      WRITE(Out_Unit,*)'bfaci= ',(bfaci(iang,karran), iang=1,nang)
      WRITE(Out_Unit,*)'parityf2= ',(parityf2(iang,karran), iang=1,nang)
      WRITE(Out_Unit,*)'dzdthd= ',(dzdthd(iang,karran),iang=1,nang)
      WRITE(Out_Unit,*)'dlnb= ',(dlnb(iang,karran), iang=1,nang)
      WRITE(Out_Unit,*)'dxdchi= ',(dxdchi(iang,karran),iang=1,nang)
      WRITE(Out_Unit,*)'dxdtha= ',(dxdtha(iang,karran), iang=1,nang)
      WRITE(Out_Unit,*)'dthddtha= ',(dthddtha(iang,karran),  iang=1,nang)
      WRITE(Out_Unit,*)'dthddchi= ',(dthddchi(iang,karran), iang=1,nang)
   ENDDO
ENDIF

RETURN
ENDSUBROUTINE Basis
