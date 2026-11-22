SUBROUTINE intbas(nxdim,nydim,ngleg,nherm,xptg,thetad,phi,tee,lam,rho)

!  returns the surface function phi(tee,lam) at the given set of Delves
!  angles.
!  Uses the analytic basis function method (ABM).
!  Actually calculates functions.  Does not interpolate.
!     rho       current hyperradius.
!     jtot      total angular momentum.
!     mega      Projection of the total angular momentum on
!               the body-fixed z-axis.
!     parity    parity of the surface functions.
!     vect      contains the coeffs of the surf fcns in the prim. basis.
!     phi       surface function phi(tee,lam) at the quadrature points.
!     nx        no. Gauss_Legendre quad pts.  now same in all channels
!     nxdim     dimension of arrays using nx
!     nherm(3)  no. Gauss_Hermite quad pts in each channel.
!     nydim     dimension ge largest of nherm
!     xptg         Gauss_Legendre quad pts. xptg=cos(delves fat theta)
!     thetad    array of delves little theta in each arrangement.
!
!this routine is called by: AphDel_Old or AphDel_New
!this routine calls
!         popt,readinbb
!-----------------------------------------------------------------------
USE FileUnits_Module
USE Narran_Module
USE Numeric_Kinds_Module
!USE gaussb_Module
USE totj_Module
USE quantb_Module
USE Parms_Module
IMPLICIT NONE
SAVE
INTEGER NXDim, NYDim
REAL(Kind=WP_Kind) rho, one, xptg(mxglegn, narran), thetad(maxhermt, narran), phi(nxdim, nydim, narran)
INTEGER  nherm(narran), ngleg(narran),nthcall,  tee, lam
INTEGER, ALLOCATABLE:: nbasis(:),nbasiss(:)  !nbasis(0:mxmega), nbasiss(0:mxmega)
REAL(Kind=WP_Kind), ALLOCATABLE:: vect(:,:)       ! vect(mxbasis, mxbasis)

!            E X T E R N A L S
EXTERNAL popt, sfunbasb, readinbb, setbassb

!          P A R A M E T E R S
PARAMETER (one=1.0d0)

DATA nthcall /0/
!-----------------------------------------------------------------------
! Determine printing options.
!-----------------------------------------------------------------------
LOGICAL little, medium, full
INTEGER ithsub, ithcall
DATA ithcall /0/, ithsub /0/
DATA little /.false./, medium /.false./, full /.false./
CALL popt ('intbas  ', little, medium, full, ithcall, ithsub)
IF(.NOT.ALLOCATED(nbasis))THEN
ALLOCATE(nbasis(0:mxmega),nbasiss(0:mxmega))
ENDIF
IF(.NOT.ALLOCATED(Vect))THEN
ALLOCATE(Vect(MxBasis,MxBasis))
ENDIF
!  --------------------------------------------------------------
!  do setups and reads only on first CALL
!  --------------------------------------------------------------
IF(nthcall/=0) GOTO 3
nthcall=1
!-----------------------------------------------------------------------
! This routine reads in all the input data
!-----------------------------------------------------------------------
CALL readinbb

!-----------------------------------------------------------------------
! Set up the basis
!-----------------------------------------------------------------------
CALL setbassb (nbasis, ngleg, nherm, nbasiss)
!-----------------------------------------------------------------------
! Calculate the Surface functions at rho.
!-----------------------------------------------------------------------
3 CONTINUE
  
CALL sfunbasb(rho, nbasis(lam), jtot, parity, tee, lam, narran, phi, ngleg, nherm,&
              xptg,thetad,vect,nxdim,nydim,nbasiss(lam))

!  --------------------------------------------------------------
RETURN

ENDSUBROUTINE intbas
