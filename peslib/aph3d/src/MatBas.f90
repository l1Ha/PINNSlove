SUBROUTINE MatBas (rhocent, nbasis, jtot, parity, megamin, megamax, nang, chanl, nangch,   &
                   phi1, v, s, ThetAPH, chiaph, eigen, deltarho, vect1, amatrx, intwt,     &
                   iunit, phi2, vect2, sthaph, cthaph, neigmin, nbasiss, jrot)
!

!            P U R P O S E   O F    S U B R O U T I N E
!   This routine calculates coupling matrix elements for use in the
!   APH3D propagation program.


!            I N P U T    A R G U M E N T S
!  rhocent  current hyperradius.
!  rholast  previous hyperradius.
!  rhonext  next hyperradius.
!  nrhosec  number of rho values to calculate potentials matrix
!           elements in the sector.
!  nbasis   number of surface functions calculated.
!  jtot     total angular momentum.
!  parity   parity.
!  megamin  minimum projection of the total angular momentum on the
!           body-fixed z-axis.
!  megamax  maximum projection of the total angular momentum on the
!           body-fixed z-axis.
!  nang     total number of angles.
!  narran   number of arrangement channels.
!  chanl    array containing the channel number for each basis
!           function.
!  nangch   array containing the number of quadrature points used
!           in each of the arrangement channel integrations.
!  phi      phi(nang,nbasis) primitive basis at each of angles.
!  v        array for storing the potential matrix elements.
!  ThetAPH  list of APH theta angles.
!  chiaph   list of APH chi angles.

!            O U T P U T    A R G U M E N T S
!  v        array for storing the potential matrix elements.
!  s        overlap matrix.

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

USE Numeric_Kinds_Module
USE Parms_Module
USE FileUnits_Module
USE Masses_Module
USE SfnDis_Module
USE CSBasis_Module
USE Boundary_Module
USE Numbers_Module
IMPLICIT NONE
INTEGER jtot, parity, megamin, megamax
INTEGER nangch(Narran+1), nang, krho, iang, ibasis, jbasis
INTEGER ipfac, mega, megap, megapmax, type, iunit
INTEGER unit, case
INTEGER case_min, case_max
INTEGER chanl(mxbasis, 0:mxmega), nbasis(0:mxmega)
INTEGER neigmin(0:mxmega), nbasiss(0:mxmega)
INTEGER jrot(mxbasis,0:mxmega), intwt(narran)
LOGICAL mid_zero
CHARACTER(LEN=11) ident
REAL(Kind=WP_Kind) rho, rhocent, pot, potcent
REAL(Kind=WP_Kind) rholast, rhonext, epslon
REAL(Kind=WP_Kind) deltarho, sqdi, sqti, xlamp, xlampp
REAL(Kind=WP_Kind) tmursi, xtra, xtrafac,temp
REAL(Kind=WP_Kind) ScaleFac, aterm, bterm, cterm, fac
REAL(Kind=WP_Kind) phi1(nang, mxbasis), v(nang), s(mxbasiss,mxbasiss)
REAL(Kind=WP_Kind) ThetAPH(nang), chiaph(nang),phi2(nang,mxbasis)
REAL(Kind=WP_Kind) vect1(mxbasiss,mxbasiss), amatrx(mxbasiss,mxbasiss)
REAL(Kind=WP_Kind) eigen(mxbasiss,0:mxmega),rhos(NRhoSec), cthaph(nang)
REAL(Kind=WP_Kind) vect2(mxbasiss, mxbasiss), sthaph(nang)
REAL(Kind=WP_Kind) stemp(mxbasiss,mxbasiss),vectemp(mxbasiss,mxbasiss)
REAL(Kind=WP_Kind) subd(mxbasiss),eigentemp(mxbasiss)
EXTERNAL potelmnt, potaph, epslon, rdwrsurb, matelmb

!-----------------------------------------------------------------------
! Determine printing options.
!-----------------------------------------------------------------------
LOGICAL little, medium, full
INTEGER ithcall, ithsub
DATA ithcall /0/, ithsub /0/
DATA little /.false./, medium /.false./, full /.false./
DATA mid_zero/.true./
CALL popt ('matbas  ', little, medium, full, ithcall, ithsub)

!  --------------------------------------------------------------
!  set some factors needed.
!  --------------------------------------------------------------
IF(jtot/=0)mid_zero=.false.
ipfac=(-1)**(jtot+parity)
sqti=one/sqrt(two)
!
!  --------------------------------------------------------------
!  evaluate matrix elements at NRhoSec rho values on the sector.
!  This basis code spaces sector centers logarithmically.
!  The next few statements space the sector boundaries linearly.
!  They can be commented out at some later time.
!  --------------------------------------------------------------
rholast=rhocent/(one+deltarho)
rhonext=rhocent*(one+deltarho)
!CALL rhocalc(rhomin, deltarho1, deltarho2, ithrho, rho, rholast, rhonext)
!rhos(1)=(rholast+rhocent)/two
!rhos(2)=rhocent
!rhos(3)=(rhocent+rhonext)/two

!---
IF(IthRho==1)THEN
   RhoLast=rhocent/2.d0
   rhos(1)=RhoLast
ELSE
   RhoLast=Basis_Dist(IthRho-1)
   rhos(1)=Sector_Boundaries(IthRho)
ENDIF
RhoCent=Basis_Dist(IthRho)
RhoNext=Basis_Dist(IthRho+1)
CurrentSector=IthRho
rhos(2)=Basis_Dist(IthRho)
rhos(3)=Sector_Boundaries(IthRho+1)
IF(rhocent<rhos(1).or.rhocent>rhos(3))THEN
   WRITE(Out_Unit,'(1x,A,4E14.7)')'Error in MatBas, Rhos=',rhos,rhocent
   WRITE(Msg_Unit,'(1x,A,4E14.7)')'Error in MatBas, Rhos=',rhos,rhocent
   !!!!STOP 'Error in MATBAS'
ENDIF
IF(Medium)WRITE(Out_Unit,*)RhoLast,Rhos,RhoNext
!--

IF(NRhoSec/=3) STOP 'matbas1'
!  --------------------------------------------------------------
!  the next few statements space the sector boundaries logarithmically,
!  and IF NRhoSec is odd, assure that the central rho value is rhocent.
!  they can be unblocked when the propagator is ready to work the
!  same way.
!  ---------------------------------------------------------------
!     sqdi=one/sqrt(one+deltarho)
!     DO  krho = 1, NRhoSec
!        rho = rhocent*sqdi*(one+deltarho)**((krho-one)/(NRhoSec-one))
!        rhos(krho)=rho
!     ENDDO

!  --------------------------------------------------------------
!  start loop on mega
!  --------------------------------------------------------------
DO mega=megamin, megamax
   REWIND phip_unit
   REWIND phid_unit
   REWIND iunit

!  --------------------------------------------------------------
!  WRITE preliminary info on unit pmat_unit (Pmatrx)
!  --------------------------------------------------------------
   WRITE(pmat_unit) rhocent, neigmin(mega), NRhoSec, (rhos(krho),krho=1,NRhoSec), mid_zero
   WRITE(pmat_unit) (eigen(ibasis,mega), ibasis=1, neigmin(mega))

!  --------------------------------------------------------------
!  angular momentum raising factors
!  --------------------------------------------------------------
   xlamp=0.
   xlampp=0.
   IF(mega<jtot)THEN
      xlamp=sqrt(jtot*(jtot+one)-mega*(mega+one))
   ENDIF
   IF(mega<(jtot-1))THEN
      xlampp=sqrt(jtot*(jtot+one)-(mega+1)*(mega+two))
   ENDIF

!  --------------------------------------------------------------
!  start loop on megap
!  --------------------------------------------------------------
   megapmax=min(jtot,megamax,mega+2)
DO megap=mega, megapmax

!  --------------------------------------------------------------
!  the following is for megap>=mega only.
!  IF megap=mega, calculate potential matrix elements and (for cstest)
!  diagonal part of asymmetric top centrifugal term.
!  --------------------------------------------------------------
   IF(megap==mega)THEN
      type=2
      fac=one
      CALL rdwrsurb(1, iunit, mega, megamin, nang, nbasis(mega),  neigmin(mega), phi1, phi2, vect1, type, &
                    mxbasis, nbasiss(mega), mxbasiss)

!  --------------------------------------------------------------
!  IF megap=mega+1, calculate Coriolis coupling matrix elements.
!  Here it reads functions which are the derivative w.o. chi of
!  the primitive basis.
!  --------------------------------------------------------------
   ELSEIF(megap==mega+1)THEN
      type=3
      fac=xlamp
      IF(mega==0) fac=fac*(1+ipfac)*sqti
      CALL rdwrsurb(1, iunit, megap, megamin, nang, nbasis(megap), neigmin(megap), phi1, phi2, &
                    vect2, type, mxbasis, nbasiss(megap), mxbasiss)

!  --------------------------------------------------------------
!  IF megap=mega+2, calculate asymmetric top coupling matrix elements.
!  --------------------------------------------------------------
   ELSEIF(megap==mega+2)THEN
      type=4
      fac=xlamp*xlampp/four
      IF(mega==0) fac=fac*(1+ipfac)*sqti
      CALL rdwrsurb(1, iunit, megap, megamin, nang, nbasis(megap), neigmin(megap), phi1, phi2, &
                    vect2, type, mxbasis, nbasiss(megap), mxbasiss)

!  --------------------------------------------------------------
!  IF NONE of the above are satisfied, there is no coupling.
!  don't waste time or disk space calculating or writing anything.
!  --------------------------------------------------------------
   ELSE
      type=5
      fac=0.
      WRITE(Out_Unit,*) 'mega, megap=',mega, megap
      STOP 'matbas'
   ENDIF

!  --------------------------------------------------------------
!  start loop on rho values on sector
!  --------------------------------------------------------------
DO krho = 1, NRhoSec
   rho=rhos(krho)
   tmursi=one/(usys2*rho*rho)
IF(type==2)THEN
   case_min = 1
   IF(krho==1)THEN
      case_max = 3
   ELSE
      case_max = 1
   ENDIF
ELSE
   case_min = 1
   case_max = 1
ENDIF

IF(type/=2.AND.krho/=1)GOTO 8
DO case = case_min, case_max
!  --------------------------------------------------------------
!  potential matrix elements. type=2
!  --------------------------------------------------------------
IF(type==2)THEN
   IF(ABS(rho-rhocent)>3.0*epslon(rho).or.(mega==1 .AND..NOT.cstest).or.case/=1)THEN
   IF(case==1)THEN
      DO iang = 1, nang
         CALL potaph (rho, ThetAPH(iang), chiaph(iang), pot)
         CALL potaph (rhocent, ThetAPH(iang), chiaph(iang), potcent)
         v(iang) = pot - (rhocent/rho)**2 *potcent
!  next card is special case for testing.
!              v(iang) = one
      ENDDO
   ELSEIF(case==2)THEN
      DO iang = 1, nang
         aterm = two/(one+sthaph(iang))
         bterm = one/sthaph(iang)**2
         v(iang) = tmursi*(aterm+bterm)/two
      ENDDO
   ELSEIF(case==3)THEN
      DO iang = 1, nang
         aterm = two/(one+sthaph(iang))
         bterm = one/sthaph(iang)**2
         IF(sthaph(iang)/=one)THEN
            cterm = two/(one-sthaph(iang))
         ELSE
            cterm=zero
         ENDIF
         v(iang) = tmursi*(cterm-(aterm+bterm)/two)
      ENDDO
   ENDIF
!  add diagonal part of asymmetric top term, IF necessary.
!  xtra is -(-1)**(j+p)*j*(j+1)*(a-b)/4
      IF(mega==1.AND..NOT.cstest)THEN
         xtrafac=-tmursi*jtot*(jtot+one)*ipfac/four
         DO iang = 1, nang
            xtra=two/(one+sthaph(iang))-one/sthaph(iang)**2
            xtra=xtrafac*xtra
            v(iang)=v(iang)+xtra
         ENDDO
      ENDIF

!  -------------------------------------------------------------------
!  form matrix elements over the primitive basis.
!  ------------------------------------------------------------------
      WRITE(Out_Unit,*)"Primitive Basis Matrix Elements"
      CALL potelmnt (nang, nbasis(mega), chanl(1,mega), nangch, phi1, v, s, rho, intwt, &
                     mxbasis, nbasiss(mega), mxbasiss, jrot(1,mega))

!  -------------------------------------------------------------------
!  DO similarity transform to get matrix elements over surface fcns.
!  first form a=s*vect1. THEN, form s=vect1(trans)*a
!  the new s is THEN the matrix elements over surface fcns.
!  --------------------------------------------------------------------
!GAPFIX 2010
      CALL dgemm('N','N',nbasiss(mega),nbasiss(mega),neigmin(mega),1.d0,s, mxbasiss,vect1,mxbasiss, 0.d0,amatrx,mxbasiss)
      CALL dgemm('T','N',neigmin(mega),nbasiss(mega),neigmin(mega), 1.d0, vect1,mxbasiss,amatrx, mxbasiss,0.d0, s,mxbasiss)
      
      
      !CALL dgemm('N','N',nbasiss(mega),neigmin(mega),nbasiss(mega),1.d0,s, mxbasiss,vect1,mxbasiss, 0.d0,amatrx,mxbasiss)
      !CALL dgemm('T','N',neigmin(mega),neigmin(mega),nbasiss(mega), 1.d0, vect1,mxbasiss,amatrx, mxbasiss,0.d0, s,mxbasiss)
   ELSE
!     bypass calculation at rhocent where v is zero.
      s=zero
   ENDIF

!  --------------------------------------------------------------
!  Coriolis matrix elements. type=3.
!  for symmetric systems, this uses jeven to keep the jrot's even or
!  odd in the incident arrangement, and performs matrix elements
!  coupling functions of different reflection parity in chi via
!  the derivative w.o. chi. rtp 92/8.
!  --------------------------------------------------------------
ELSEIF(type==3)THEN
   IF(krho==1)THEN
   DO iang=1, nang
      v(iang)=-fac*tmursi*cthaph(iang)/sthaph(iang)**2
   ENDDO
!  form matrix elements over primitive basis.
   CALL matelmb(nang, nbasis(mega), nbasis(megap), nangch, phi1, v, phi2, s, rho, intwt,                        &
                chanl(1,mega), chanl(1,megap), nbasiss(mega), nbasiss(megap), jrot(1,mega), jrot(1,megap))
!  similarity transform.  snew=vect1(trans)*s*vect2 is for surface fcns.
      CALL dgemm('N','N',nbasiss(mega),nbasiss(megap),neigmin(megap),1.d0,s, mxbasiss,vect2,mxbasiss, 0.d0, amatrx,mxbasiss)
      CALL dgemm('T','N',neigmin(mega),nbasiss(mega),neigmin(megap), 1.d0, vect1,mxbasiss,amatrx, mxbasiss, 0.d0, s,mxbasiss)
   ELSE
!  after the first rho on the sector, get others by scaling as these
!  matrix elements just ScaleFac as rho**-2
      ScaleFac=(rhos(krho-1)/rhos(krho))**2
      s=ScaleFac*s
   ENDIF

!  --------------------------------------------------------------
!  asymmetric top matrix elements. type=4
!  --------------------------------------------------------------
ELSEIF(type==4)THEN
   IF(krho==1)THEN
   DO iang=1, nang
      v(iang)=tmursi*(two/(one+sthaph(iang))-one/sthaph(iang)**2)
!  at the above line v=a-b
      v(iang)=fac*v(iang)
   ENDDO
!  form matrix elements over primitive basis.
   CALL matelmb(nang, nbasis(mega), nbasis(megap), nangch, phi1, v, phi2, s, rho, intwt,                         &
                chanl(1,mega), chanl(1,megap), nbasiss(mega), nbasiss(megap), jrot(1,mega), jrot(1,megap))
!  similarity transform.  snew=vect1(trans)*s*vect2 is for surface fcns.
      CALL dgemm('N','N',nbasiss(mega),nbasiss(megap),neigmin(megap),1.d0,s, mxbasiss,vect2,mxbasiss, 0.d0, amatrx,mxbasiss)
      CALL dgemm('T','N',neigmin(mega),nbasiss(mega),neigmin(megap), 1.d0, vect1,mxbasiss,amatrx, mxbasiss, 0.d0, s,mxbasiss)
   ELSE
!  after the first rho on the sector, get others by scaling as these
!  matrix elements just ScaleFac as rho**-2
      ScaleFac=(rhos(krho-1)/rhos(krho))**2
      s=ScaleFac*s
   ENDIF
ELSE
   WRITE(Out_Unit,*)'type=', type
   STOP 'matbasty'
ENDIF
!-----------------------------------------------------------------------
! WRITE out the coupling matrix.
!  tape pmat_unit is 'PotMatrx_ABM'              type=2   case=1
!  tape jtot_unit is 'diag_jtot_ABM'           type=2   case=2
!  tape lam_unit is 'diag_lambda_ABM'          type=2   case=3
!  tape cor_unit is 'coriolis_ABM'             type=3
!  tape asym_unit is 'asym_top_ABM'            type=4
!-----------------------------------------------------------------------

   IF(type==2)THEN
      IF(case==1)THEN
         unit=pmat_unit
         ident = 'Pmatrx     '
      ELSEIF(case==2)THEN
         unit=jtot_unit
         ident = 'diag_jtot  '
      ELSEIF(case==3)THEN
         unit=lam_unit
         ident = 'diag_lambda'
      ELSE
         WRITE(Out_Unit,*)'Error in case',case
         STOP 'matbas'
      ENDIF
   ELSEIF(type==3)THEN
      unit=cor_unit
      ident = 'coriolis   '
   ELSEIF(type==4)THEN
      unit=asym_unit
      ident = 'asym_top   '
   ENDIF

   IF(type/=2)THEN
      WRITE(unit)rhocent,rhos(krho),mega,megap, neigmin(mega),neigmin(megap),ident
      DO jbasis=1,neigmin(megap)
         WRITE(unit)(s(ibasis,jbasis),ibasis=1,neigmin(mega))
      ENDDO
   ELSE
      IF(jtot/=0.or.case==1)THEN
         IF(ABS(rho-rhocent)>3.0*epslon(rho).or.(mega==1 .AND..NOT.cstest).or.case/=1)THEN
            IF(case/=1)THEN
               WRITE(unit)rhocent,rhos(krho),mega,megap, neigmin(mega),neigmin(megap),ident
            ENDIF
            DO jbasis=1,neigmin(megap)
               WRITE(unit)(s(ibasis,jbasis),ibasis=jbasis,neigmin(mega))
            ENDDO
         ENDIF
      ENDIF
   ENDIF
   IF(little)THEN
      WRITE(Out_Unit,*)'routine matbas'
      WRITE(Out_Unit,*)'nbasiss = ', nbasiss(mega), nbasiss(megap)
   ENDIF
   IF(medium.AND..NOT.full)THEN
      WRITE(Out_Unit,*)'Type of coupling matrix = ',ident
      WRITE(Out_Unit,*)'rhocent,rhos(krho),mega,','megap,neigmin(mega),','neigmin(megap),ident'
      WRITE(Out_Unit,*)rhocent,rhos(krho),mega, megap,neigmin(mega),neigmin(megap),ident
      WRITE(Out_Unit,*) 'Part of the coupling matrix at rho = ', rho
      DO ibasis = 1, min(nbasiss(mega),10)
         WRITE(Out_Unit,*)(s(ibasis,jbasis),jbasis = 1, min(nbasiss(megap),10))
      ENDDO
   ENDIF
   IF(full)THEN
      WRITE(Out_Unit,*)'Type of coupling matrix = ',ident
      WRITE(Out_Unit,*)'rhocent,rhos(krho),mega,', 'megap,neigmin(mega),','neigmin(megap),ident'
      WRITE(Out_Unit,*)rhocent,rhos(krho),mega, megap,neigmin(mega), neigmin(megap),ident
      WRITE(Out_Unit,*)' coupling matrix at rho = ', rho
      CALL MxOut(s, mxbasiss, mxbasiss)
!  The last rows and columns of s as just written may contain junk.
ENDIF

ENDDO
ENDDO
ENDDO
  8   CONTINUE
ENDDO

RETURN
ENDSUBROUTINE MatBas
