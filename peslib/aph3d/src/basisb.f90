SUBROUTINE basisb (xptg, thlitdel, pn, ngleg, nherm, oscil,parity, parityfn, bfaci, rho, mega)

!            P U R P O S E   O F    S U B R O U T I N E
!   This routine calculates Legendre and harmonic oscillator basis
!   functions for use in AphDel_Old or AphDel_New transformation.
!   Assumes channel, not alternate, quadrature.

!            I N P U T    A R G U M E N T S
!  narran   number of arrangement channels.
!  ngleg   number of Gauss_Legendre quadrature points for each
!           arrangement channel.
!  nlegndre maximum order of the Legendre polynomial for each
!           arrangement channel.
!  nherm   number of Gauss_Hermite quadrature points for each
!           arrangement channel.
!  noscil   maximum order of the oscillator basis for each
!           arrangement.
!  parity   parity of the surface-functions.
!  xptg     Gauss_Legendre quadrature points, cos(bigtheta).
!  thlitdel array of delves theta angles.

!            O U T P U T    A R G U M E N T S
!  pn       Legendre polynomial basis.
!  oscil    Harmonic oscillator basis.
!  parityfn parity function. This function is 1 for even parity and
!           cos(chi) for odd parity.
!  bfaci    coord trans factor that goes with sho fcns of zlit.
!
!this routine is called by:
!         sfunbasb
!this routine calls
!         popt,equilth,legendre,harmonic
!-----------------------------------------------------------------------

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>






!                I N C L U D E   F I L E S
USE fileunits_Module
USE Numeric_Kinds_Module
USE Narran_Module
USE PES_MODULE
USE gaussb_Module
USE Masses_Module
USE MassFactor2_Module
USE Integrat_Module

IMPLICIT NONE
!            I N T E G E R S
INTEGER karran, iglegn, ihermt,&
       parity, ngleg, nherm,&
       n, l, mega, inos

!            R E A L S
REAL(Kind=WP_Kind) xptg, thlitdel, ct, pn, oscil,&
     parityfn,&
     rho,  tanth,&
     chi, sinths, cosths, bfac, bfaci, zlit,&
     sinthd, costhd, azsbzc,&
     y, s2td, c2td, sqr,&
     cos4,half,zero,one,two,four,six,&
     beta, eps, sin4, cos, atan, sin, sqrt, tan, asin


!            D I M E N S I O N S
DIMENSION xptg(mxglegn, narran), thlitdel(maxhermt, narran), pn(0:mxl,mxglegn,narran),oscil(0:mxn,maxhermt,narran),&
         parityfn(mxglegn, maxhermt,narran), bfaci(maxhermt,narran), ngleg(narran), nherm(narran)
    REAL(Kind=WP_Kind), ALLOCATABLE :: thliteq(:),az(:),bz(:),cz(:),&
                              se(:),cparm(:),rhoprm(:),rhofix(:)

!            I N T R I N S I C    F U N C T I O N S
INTRINSIC cos, atan, sin, sqrt, tan, asin

!            E X T E R N A L S
EXTERNAL popt, legendre, harmonic

!      P A R A M E T E R S
PARAMETER (half=0.5d0, zero=0.0d0, one=1.0d0)
PARAMETER (two=2.0d0, four=4.0d0, six=6.0d0)

!-----------------------------------------------------------------------
! Determine printing options.
!-----------------------------------------------------------------------
LOGICAL little, medium, full, expfit, tmud2, lifhb, fh2t5
INTEGER ithcall, ithsub
DATA ithcall /0/, ithsub /0/
DATA little /.false./, medium /.false./, full /.false./
CALL popt ('basisb  ', little, medium, full, ithcall, ithsub)

!-----------------------------------------------------------------------
!  old_way=.true. gives the parametrization that was used from 89/11 to
!  92/5.  old_way=.false. gives the parametrization used since 92/5.
!  Note.  The following parametrization must match that in abm/basis.F.
!-----------------------------------------------------------------------
WRITE(Out_Unit,*)'begining of routine basisb'
IF(Medium)WRITE(Out_Unit,*)'old_way=',old_way
WRITE(Msg_Unit,*)'old_way=',old_way
IF(old_way)THEN
IF(PES_Name=='tmud2 pairwise add. ')THEN
   tmud2 = .true.
   WRITE(Out_Unit,*)'tmud2=.true.'
ELSE
   WRITE(Out_Unit,*)'tmud2=.false.'
   tmud2 = .false.
ENDIF

IF(PES_Name=='LiFH_B BondOrder PES')THEN
   expfit = .true.
   lifhb = .true.
   !WRITE(Out_Unit,*)'Using exponential fit to cparm'
ELSE
   expfit = .false.
   lifhb = .false.
   !WRITE(Out_Unit,*)'Using morse fit to cparm'
ENDIF
ENDIF
IF(PES_Name=='FH2  T5A         PES')THEN
    fh2t5 = .true.
ELSE
    fh2t5 = .false.
ENDIF

!------------------------------------------------------------------------
!ALLOCATE momory
!------------------------------------------------------------------------
ALLOCATE( thliteq(narran), az(narran), bz(narran), cz(narran),&
         se(narran), cparm(narran),rhoprm(narran),&
         rhofix(narran))
!------------------------------------------------------------------------

!  determine scaled equilibrium internuclear distances.

DO karran = 1, narran
   se(karran)=rx(karran)*re(karran)/d(karran)
!  the following is used for tmu+d2 and lifh
   IF(old_way)THEN
      IF(tmud2)THEN
         rhofix(karran)=1.025d0*se(karran)
      ELSEIF(lifhb)THEN
         rhofix(karran)=1.1*se(karran)
      ENDIF
   ENDIF
ENDDO

!  ---------------------------------------------------------------------
!  determine equilibrium value of thlitdel.
!  detn parameters relating hermite quadrature variable z to thlitdel.
!  z=a*tan(theta)-b*cot(theta)+c
!  --------------------------------------------------------------------
DO  karran=1,narran
      IF(old_way)THEN
      rhoprm(karran)=rho
      ENDIF
      cparm(karran)=calpha(karran)
      IF(old_way)THEN
      IF(rho<rhofix(karran))THEN
         IF(tmud2)THEN
            rhoprm(karran)=rhofix(karran)
            cparm(1)=6.825d0 -2.4533d0*rho
!  the above parameters are for tmu+d2
         ELSEIF(lifhb)THEN
            rhoprm(karran)=rhofix(karran)
            y=exp(-balpha(3)*rhofix(3))
            cparm(3)=calpha(3)+dalpha(3)*y
         ELSE  !  the following is exp fit to cparm.
            y=exp(-balpha(karran)*rho)
            cparm(karran)=calpha(karran)+dalpha(karran)*y
         ENDIF

!  they will need to be changed IF we use basis for another system
!  in which rho<rhofix.
      ELSE
!  the following is exp fit to cparm.
         IF(expfit)THEN
            y=exp(-balpha(karran)*rho)
            cparm(karran)=calpha(karran)+dalpha(karran)*y
         ELSE
!  the following is morse fit to cparm.
           y=exp(-balpha(karran)*(rhoprm(karran)-ralpha(karran)))
           cparm(karran)=calpha(karran)+dalpha(karran)*y*(y-two)
         ENDIF
      ENDIF
      ELSE
         IF(fh2t5)THEN
!                  y=exp(-balpha(karran)*(rho-ralpha(karran))**2)
!                  zeta(karran)=dalpha(karran)*y
            WRITE(Out_Unit,*)' zeta(karran)=',zeta(karran)
         ENDIF
      ENDIF

      IF(.NOT.old_way)THEN
      CALL equilth (karran,rho,weau(karran),wexeau(karran),usys, zeta(karran),se(karran), thliteq(karran), delta(karran),beta, eps)
      beta=sqrt(beta)
      eps=eps/(six*beta)
      ELSE
         IF(se(karran)>rhoprm(karran))THEN
            WRITE(Out_Unit,*)'Error: Basisb'
            WRITE(Out_Unit,*)'se(karran) must be <', ' rhoprm(karran)'
            WRITE(Out_Unit,*)'se=', se(karran)
            WRITE(Out_Unit,*)'rhoprm', rhoprm(karran)
            STOP 'Basisb'
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
     beta=rhoprm(karran)*sqrt(usys*(weau(karran)-wexeau(karran))) *cparm(karran)
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
!  keep a positive.  If necessary, only fit 2nd deriv, not third.
        IF(az(karran)<1.d-8) az(karran)=1.d-8
      ENDIF
      cz(karran)=bz(karran)/tanth - az(karran)*tanth
      IF(az(karran)<=0.d0.or.bz(karran)<=0.d0)THEN
         WRITE(Out_Unit,*)'arrangement channel = ',karran
         WRITE(Out_Unit,*)'beta,eps,tan=',beta, eps, tanth
         WRITE(Out_Unit,*)'a,b,c= ',&
                         az(karran),bz(karran),cz(karran)
         STOP 'basisb'
       ENDIF

   ENDDO
!-----------------------------------------------------------------------
!  Determine functions of the  Delves angles.
!  In this version, functions and quadrature are over same arrangement.
!  zlit is dimensionless function of thlitdel used in sho's.
!  chi is  prin value chi.
!  bfac is a factor that goes with the Jacobian and sho functions.
!  Calculate the normalized Legendre polynomial basis.
!  Calculate the normalized harmonic oscillator basis.
!  Calculate the parity function.
!-----------------------------------------------------------------------
WRITE(*,*) "Calling Legendre from BasisB"
DO  karran = 1, narran
   DO  iglegn = 1, ngleg(karran)
      CALL legendre (nlegndre(karran), mega, pn(0,iglegn,karran), xptg(iglegn,karran))
   ENDDO
   DO  ihermt = 1, nherm(karran)
     IF(thlitdel(ihermt,karran)/=0.0)THEN
      sinthd=sin(thlitdel(ihermt,karran))
      costhd=cos(thlitdel(ihermt,karran))
      tanth=sinthd/costhd
      sinths=sinthd**2
      cosths=costhd**2
      azsbzc=one/(az(karran)*sinths+bz(karran)*cosths)
      bfac=four*sinths*cosths*sqrt(azsbzc)
      bfaci(ihermt,karran)=one/bfac
      zlit=az(karran)*tanth-bz(karran)/tanth+cz(karran)
      CALL harmonic (noscil(karran), oscil(0, ihermt, karran), zlit)
         s2td = two*sinthd*costhd
         c2td = cosths-sinths
     ELSE
         bfaci(ihermt,karran)=0d0
         DO  inos=0,noscil(karran)
           oscil(inos,ihermt,karran)=0d0
         ENDDO
      ENDIF
      DO  iglegn = 1, ngleg(karran)
         ct=xptg(iglegn,karran)
         sqr  = sqrt(c2td**2+s2td**2*ct**2)
         IF(sqr==0)THEN
           WRITE(Out_Unit,*)' Error: sqr eq 0'
           WRITE(Out_Unit,*)' Try to decrease nhermt points'
           WRITE(Out_Unit,*)'c2td,s2td,ct=',c2td,s2td,ct
           STOP 'basisb'
         ENDIF
         chi  = atan2(s2td*ct/sqr,c2td/sqr)/two
         parityfn(iglegn,ihermt,karran)=cos(chi)**parity
      ENDDO
   ENDDO
ENDDO
!-----------------------------------------------------------------------
! Print basis if desired.
!-----------------------------------------------------------------------
IF(little)THEN
   WRITE(Out_Unit,*)'routine basisb'
   WRITE(Out_Unit,*)'narran = ',narran
   WRITE(Out_Unit,*)'ngleg = ',ngleg
   WRITE(Out_Unit,*)'nlegndre = ',nlegndre
   WRITE(Out_Unit,*)'nherm = ',nherm
   WRITE(Out_Unit,*)'noscil = ',noscil
ENDIF
IF(medium)THEN
   DO  karran = 1, narran
      IF(Integrat(karran))THEN
         WRITE(Out_Unit,*)'arrangement channel = ',karran
         WRITE(Out_Unit,*)'xptg = ',(xptg(iglegn,karran), iglegn=1,ngleg(karran))
         WRITE(Out_Unit,*)'thlitdel = ',(thlitdel(ihermt,karran),ihermt=1,nherm(karran))
      ENDIF
  ENDDO
ENDIF
IF(full)THEN
   DO  karran=1, narran
      IF(Integrat(karran))THEN
         WRITE(Out_Unit,*)'arrangement channel = ',karran
         WRITE(Out_Unit,*)'a,b,c= ',az(karran),bz(karran),cz(karran)
         WRITE(Out_Unit,*)'pn = ',((pn(l, iglegn, karran), iglegn=1, ngleg(karran)), l = 0, nlegndre(karran))
         WRITE(Out_Unit,*)'oscil = ',((oscil(n, ihermt, karran), ihermt=1, nherm(karran)), n = 0, noscil(karran))
         WRITE(Out_Unit,*)'bfaci= ',(bfaci(ihermt,karran), ihermt=1, nherm(karran))
      ENDIF
   ENDDO
ENDIF
!------------------------------------------------------------------------
!DEALLOCATE momory
!------------------------------------------------------------------------
DEALLOCATE( thliteq, az, bz, cz, se, cparm, rhoprm, rhofix)
!------------------------------------------------------------------------
WRITE(Out_Unit,*)'END of routine basisb'
RETURN
ENDSUBROUTINE BasisB
