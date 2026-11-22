SUBROUTINE mat (id, th, ch, f, ff, fbar, dfbar, lambda,  energy, hamil, thetaval, chivals, naph, rholast,          &
                      rhocent, rhonext, numnp, nel, ntheta, nchi,  iunit, xpt, wht, value, nquad, h, p, hcoef)
USE Numeric_Kinds_Module
USE totj_Module
USE csbasis_Module
USE Parms_Module
USE FileUnits_Module
USE Numbers_Module
!
!   $RCSfile: mat.f,v $ $Revision: 1.15 $
!   $Date: 89/11/16 11:03:05 $
!   $State: Stable $
!

!                  P U R P O S E    O F    S U B R O U T I N E

!   This routine calculates potential, Coriolis and asymetric top
!   coupling matrix elements.

!                  I N P U T    A R G U M E N T S

!    naph         Number of surface functions at each mega.
!    rhocent      Distance  at which surface functions were
!                 calculated.
!    rholast      Distance at which the  previous functions were
!                 calculated.
!    rhonext      Distance at which the next surface functions will
!                 be calculated.
!    f            Surface functions calculated at rho.
!    id           List of node numbers at rhocent.
!    th           List of theta indices at rhocent.
!    ch           List of chi indices at rhocent.
!    thetaval     List of theta values.
!    chivals      List of chi values.
!    ntheta       Number of theta values.
!    nchi         Number of chi values.
!    fbar         Temporary storage array.
!    dfbar        Temporary storage array.
!    lambda       Temporary storage array.
!    energy       Temporary storage array.
!    xpt          Abscissas for the Gauss_Legendre quadrature.
!    wht          Weights for the Gauss_Legendre quadrature.
!    nquad        Number of Gauss_Legendre quadrature points.
!    ff
!    nunmp
!    nel
!    iunit
!    value
!    h
!    p
!    hcoef

!                 O U T P U T    A R G U M E N T S

!    hamil        Overlap matrix.

IMPLICIT NONE
CHARACTER(LEN=11) ident
LOGICAL :: little, medium, full, mid_zero
INTEGER :: th(9, nel), ch(9, nel), naph, ithcall, ithsub, megap, i, j, mxelem, ntheta, nchi, imeg, imegp
INTEGER :: nel, numnp, krho, type, iunit, ibasis, ipfac, nquad, unit, case, case_min, case_max
INTEGER :: id(9,nel), lambda(naph)
INTEGER, PARAMETER :: nrho=3
REAL(Kind=WP_Kind) :: rho, thetaval(ntheta), chivals(nchi), fac
REAL(KIND=WP_Kind) :: rholast, rhocent, rhonext
REAL(Kind=WP_Kind) :: xlamp, xlampp, sqti, xpt(*), wht(*), epslon
REAL(Kind=WP_Kind) :: fbar(naph, nquad, nquad), dfbar(naph, nquad, nquad), value(nquad,nquad)
REAL(Kind=WP_Kind) :: h(9, nquad, nquad), p(2, 9, nquad, nquad),  hcoef(9, nquad, nquad)
REAL(Kind=WP_Kind) :: energy(naph), rhos(nrho), hamil(naph, naph), f(numnp, naph), ff(*)
!-----------------------------------------------------------------------
! Determine printing options.
!-----------------------------------------------------------------------
DATA ithcall /0/, ithsub /0/
DATA little /.false./, medium /.false./, full /.false./
DATA mid_zero/.true./
CALL popt ('mat     ', little, medium, full, ithcall, ithsub)

IF(jtot/=0)mid_zero=.false.
!  factors needed
ipfac=(-1)**(jtot+parity)
sqti=1.0d0/sqrt(2.0d0)
type=5
!-----------------------------------------------------------------------
! Determine minimum lambda value.
!-----------------------------------------------------------------------
megamin = 0
IF(ipfac<0) megamin = 1

!-----------------------------------------------------------------------
! Loop over all possible lambda values.
!-----------------------------------------------------------------------
DO 50 imeg = megamin, megamax

mega=imeg
!  angular momentum raising factors.
xlamp=sqrt(jtot*(jtot+1.)-imeg*(imeg+1.))
xlampp=0.
IF(imeg<megamax)THEN
   xlampp=sqrt(jtot*(jtot+1.)-(imeg+1.)*(imeg+2.))
ENDIF
!-----------------------------------------------------------------------
! READ in surface function for imeg.
!-----------------------------------------------------------------------
   IF(imeg/=megamin)THEN
      CALL rdwrsurf (id, th, thetaval, ch, chivals, f, energy, nel, numnp, &
      naph, ntheta, nchi, iunit, 1, mxnode, mxelem, imeg, megamin)
   ENDIF


   DO imegp = imeg, megamax
   megap = imegp

!-----------------------------------------------------------------------
!  Note.  The following is all for imegp>=imeg only.
!  For normalization, all these integrals have an extra factor of 2.
!  --------------------------------------------------------------------
! IF imegp is equal to imeg calculate potential matrix elements.
!-----------------------------------------------------------------------
      IF(imegp==imeg)THEN
         type = 2
         fac = 2.0d0
         ident = 'Pmatrx     '
!-----------------------------------------------------------------------
!  IF imegp=imeg+1, calculate Coriolis coupling matrix elements.
!-----------------------------------------------------------------------
      ELSEIF(imegp==imeg+1)THEN
         type = 3
         fac = 2.0d0*xlamp
         IF(imeg==0)fac=fac*(1+ipfac)*sqti
         ident = 'coriolis   '
!-----------------------------------------------------------------------
!  IF imegp=imeg+2, calculate asymmetric top coupling matrix elements.
!-----------------------------------------------------------------------
      ELSEIF(imegp==imeg+2)THEN
         type = 4
         fac = 2.0d0*0.25d0*xlamp*xlampp
         IF(imeg==0)fac=fac*(1+ipfac)*sqti
         ident = 'asym_top   '
!-----------------------------------------------------------------------
! IF NONE of the above conditions are satisfied there is no coupling
! hence set the coupling matrix elements to zero.
!-----------------------------------------------------------------------
      ELSE
         type=5
         fac=zero
      ENDIF

!-----------------------------------------------------------------------
! READ in surface function for imegp.
!-----------------------------------------------------------------------
IF(imeg/=imegp)THEN
      CALL rdwrsurf (id, th, thetaval, ch, chivals, ff, energy, nel, numnp, &
      naph, ntheta, nchi, iunit, 1, mxnode, mxelem, imegp, megamin)
ENDIF

!-----------------------------------------------------------------------
! ff now contains the eigenfunctions which are orthonormal.
!-----------------------------------------------------------------------
   rhos(1) = (rholast+rhocent)/2
   rhos(2) = rhocent
   rhos(3) = (rhocent+rhonext)/2
WRITE(Pmatrx_FEM_Unit) rhocent,naph,nrho,(rhos(krho),krho=1,nrho),mid_zero
WRITE(Pmatrx_FEM_Unit) (energy(ibasis),ibasis=1,naph)

DO 44 krho = 1,nrho
   rho = rhos(krho)
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

IF(type/=2.AND.krho/=1)GOTO 44
DO 88 case = case_min, case_max

!-----------------------------------------------------------------------
! Calculate matrix elements for three rho values within the sector.
!-----------------------------------------------------------------------

IF(rho/=rhocent.or.jtot/=0.AND.type<5)THEN
   IF(little) WRITE(Out_unit,*)'imeg,imegp,nel,naph =',imeg,imegp, nel, naph
   IF(imeg/=imegp)THEN
      CALL matxelem (nel, naph, f, th, ch, naph, ff, numnp, id, fbar, dfbar, thetaval, chivals, type,  &
                           ntheta, nchi, hamil, rho, rhocent,xpt, wht, value, nquad, h, p, hcoef)
   ELSE
      CALL matxelem (nel, naph, f, th, ch, naph, f,  numnp, id, fbar, dfbar, thetaval, chivals, type,  &
                           ntheta, nchi, hamil, rho, rhocent,  xpt, wht, value, nquad, h, p, hcoef)
      CALL scaleit(naph*naph, hamil, 1, two)
   ENDIF
ELSE
   CALL vsets (naph*naph, hamil, 1, zero)
   mid_zero = .true.
ENDIF

!------------------------------------------------------------------------
! WRITE out coupling matrix.
! tape Pmatrx_FEM_Unit is 'Pmatrx_FEM
! tape Diag_Jtot_FEM_Unit is 'diag_jtot_FEM'
! tape Diag_Lambda_FEM_Unit is 'diag_lambda_FEM'
! tape Asym_Top_FEM is 'asym_top_FEM'
! tape Coriolis_FEM_Unit is 'coriolis_FEM'
!------------------------------------------------------------------------
IF(type==2)THEN
   IF(case==1)THEN
      unit=Pmatrx_FEM_Unit
      ident = 'Pmatrx     '
   ELSEIF(case==2)THEN
      unit=Diag_Jtot_FEM_Unit
      ident = 'diag_jtot  '
   ELSEIF(case==3)THEN
      unit=Diag_Lambda_FEM_Unit
      ident = 'diag_lambda'
   ELSE
      WRITE(Out_unit,*)'Error in case',case
      STOP 'matbas'
   ENDIF
ELSEIF(type==3)THEN
   unit=Asym_Top_FEM_Unit
   ident = 'asym_top   '
ELSEIF(type==4)THEN
   unit=Coriolis_FEM_Unit
   ident = 'coriolis   '
ENDIF

IF(type/=2)THEN
   WRITE(unit) rhocent, rhos(krho), mega, megap,  naph, naph, ident
   DO 100 j  = 1, naph
      WRITE(unit) (hamil(i, j), i = 1, naph)
  100    CONTINUE
ELSE
   IF(jtot/=0.or.case==1)THEN
      IF(ABS(rho-rhocent)>3.0*epslon(rho).or.(mega==1 .AND. .NOT.cstest).or.case/=1)THEN
         IF(case/=1)THEN
            WRITE(unit)rhocent,rhos(krho),mega,megap, naph,naph,ident
         ENDIF
         DO 101 j =1,naph
            WRITE(unit) (hamil(i,j), i=j,naph)
  101          CONTINUE
      ENDIF
   ENDIF
ENDIF

IF(full)THEN
   WRITE(Out_unit, 160) rho, rhocent
   CALL MxOutL(hamil, naph, naph, 0, 'aph  ', 'aph  ')
ENDIF

   88          CONTINUE
   44       CONTINUE
         ENDDO
   50 CONTINUE
RETURN
!
  160 FORMAT(/, 'potential matrix elements: rho, rhocent',2es14.7)
END
