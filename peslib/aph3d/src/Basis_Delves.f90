SUBROUTINE Basis_Delves(z,thetd,az,bz,cz,karran)
USE fileunits_Module
USE numeric_kinds_Module
USE Narran_Module
USE PES_MODULE
USE GaussQuady_Module
USE gaussb_Module
USE Masses_Module
USE MassFactor2_Module
USE Integrat_Module
USE rhonum_module
USE thetas_MODULE
USE Numbers_Module
IMPLICIT NONE
INTEGER iglegn, ihermt, n, l,inos,i
INTEGER ic,iy,in,karran

REAL(Kind=WP_Kind)  thlitdel, ct,&
      parityfn,&
      rho,  tanth,&
      chi, sinths, cosths,z(*),  zlit(maxhermt),&
      sinthd, costhd, azsbzc,&
       y, s2td, c2td, sqr,&
      cos4, beta, eps, sin4, cos, sin, tan, asin,pi2,thetd(maxhermt,Narran)
REAL(Kind=WP_Kind) zmc,a(3),x(3),az(3),bz(3),cz(3)

dimension thlitdel(maxhermt, narran)

REAL(Kind=WP_Kind), allocatable :: thliteq(:), se(:), cparm(:), rhoprm(:), rhofix(:)
REAL(Kind=WP_Kind),allocatable :: sqci(:)
INTRINSIC cos, atan, sin, sqrt, tan, asin
REAL(Kind=WP_Kind) atan2, max, sign
INTRINSIC atan2, max, sign
!-----------------------------------------------------------------------
! Determine printing options.
!-----------------------------------------------------------------------
LOGICAL little, medium, full, expfit, tmud2, lifhb, fh2t5
INTEGER ithcall, ithsub, ithcal2
DATA ithcall /0/, ithcal2 /0/, ithsub /0/
IF(old_way)THEN
   IF(PES_Name=='tmud2 pairwise add. ')THEN
      tmud2 = .true.
      WRITE(Out_Unit,*)'tmud2=.true.'
   ELSE
      !WRITE(Out_Unit,*)'tmud2=.false.'
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
IF(.NOT.allocated(thliteq))THEN
   ALLOCATE(thliteq(narran))
   ALLOCATE(se(narran))
   ALLOCATE(cparm(narran))
   ALLOCATE(rhoprm(narran))
   ALLOCATE(rhofix(narran))
   ALLOCATE(sqci(narran))
ENDIF
pi2=2.d0*atan(1.d0)
rho=rhonum
!Step1 try to get thlitdel
ic=karran


!Step2 Get thliteq

!  determine scaled equilibrium internuclear distances.
   se(karran)=rx(karran)*re(karran)
   IF(old_way)THEN
      IF (tmud2)THEN  !  the following is used for tmu+d2
         rhofix(karran)=1.025d0*se(karran)
      ELSE IF (lifhb)THEN !  the following is used for LiFH
         rhofix(karran)=1.1*se(karran)
      ENDIF
   ENDIF

!  ---------------------------------------------------------------------
!  determine equilibrium value of thlitdel.
!  detn parameters relating hermite quadrature variable z to thlitdel.
!  z=a*tan(theta)-b*cot(theta)+c
!  --------------------------------------------------------------------
      IF(old_way)THEN
      rhoprm(karran)=rho
      ENDIF
      cparm(karran)=calpha(karran)
      IF(old_way)THEN
      IF(rho<rhofix(karran))THEN
         IF (tmud2)THEN
            rhoprm(karran)=rhofix(karran)
            cparm(1)=6.825d0 -2.4533d0*rho
!  the above parameters are for tmu+d2
         ELSE IF (lifhb)THEN
            rhoprm(karran)=rhofix(karran)
            y=exp(-balpha(3)*rhofix(3))
            cparm(3)=calpha(3)+dalpha(3)*y
         ELSE
!  the following is exp fit to cparm.
            y=exp(-balpha(karran)*rho)
            cparm(karran)=calpha(karran)+dalpha(karran)*y
         ENDIF

!  they will need to be changed IF we use basis for another system
!  in which rho<rhofix.
      ELSE
!  the following is exp fit to cparm.
         IF (expfit)THEN
            y=exp(-balpha(karran)*rho)
            cparm(karran)=calpha(karran)+dalpha(karran)*y
         ELSE
!  the following is morse fit to cparm.
           y=exp(-balpha(karran)*(rhoprm(karran)-ralpha(karran)))
           cparm(karran)=calpha(karran)+dalpha(karran)*y*(y-two)
         ENDIF
      ENDIF
      ELSE
         IF (fh2t5)THEN
                  y=exp(-balpha(karran)*(rho-ralpha(karran))**2)
                  zeta(karran)=dalpha(karran)*y
            WRITE(Out_Unit,*)' zeta(karran)=',zeta(karran)
         ENDIF
      ENDIF
      IF(.NOT.old_way)THEN
      CALL equilth (karran,rho,weau(karran),wexeau(karran),usys, zeta(karran),se(karran), thliteq(karran), delta(karran),beta, eps)
      beta=sqrt(beta)
      eps=eps/(six*beta)
      ELSE
         IF (se(karran)>rhoprm(karran))THEN
            WRITE(Out_Unit,*)'Error: Basis_Delves'
            WRITE(Out_Unit,*)'se(karran) must be <', ' rhoprm(karran)'
            WRITE(Out_Unit,*)'se=', se(karran)
            WRITE(Out_Unit,*)'rhoprm', rhoprm(karran)
            STOP 'Basis_Delves'
         ENDIF

      thliteq(karran)=asin(se(karran)/rhoprm(karran))
      ENDIF
!Step3 Get zlit
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
!  keep a positive.  IF necessary, only fit 2nd deriv, not third.
        IF(az(karran)<1.d-8) az(karran)=1.d-8
      ENDIF
      cz(karran)=bz(karran)/tanth - az(karran)*tanth
      IF(az(karran)<=0.d0.or.bz(karran)<=0.d0)THEN
         WRITE(Out_Unit,*)'arrangement channel = ',karran
         WRITE(Out_Unit,*)'beta,eps,tan=',beta, eps, tanth
         WRITE(Out_Unit,*)'a,b,c= ',&
                         az(karran),bz(karran),cz(karran)
         STOP 'basis_delves'
       ENDIF
   a(1)=az(karran)
   a(3)=-bz(karran)
      DO ihermt = 1, nhermt(karran)
         zmc = xpth(ihermt, karran) - cz(karran)
!         STOP 'warning'
!module of gaussquad has been changed to be different from gaussb
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
         thetd(ihermt,karran)=atan(tanth)
      ENDDO

RETURN
ENDSUBROUTINE Basis_Delves
