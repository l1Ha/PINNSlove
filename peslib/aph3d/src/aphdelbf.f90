SUBROUTINE aphdelbf(rho, naph, ndel, tstore, umat, temp, nidm)

!
!-----------------------------------------------------------------------
!     This routine is a buffer to AphDel_Old or AphDel_New.
!     It's purpose is to ALLOCATE and DEALLOCATE the arrays needed by 
!     AphDel_Old or AphDel_New on the first energy.
!-----------------------------------------------------------------------
!     Written by B.J. Archer
!
!this routine is called by: aph3d
!this routine calls: popt, AphDel_Old or AphDel_New
!-----------------------------------------------------------------------
USE fileunits_Module
USE Numeric_Kinds_Module
USE Narran_Module
USE memadd_Module
USE Storage_Module
USE rotmax_Module
USE totj_Module
USE Arrch_Module
!USE Gaussb_Module
USE GaussQuady_Module
USE thetas_Module
USE Integrat_Module
USE quantb_Module
USE Parms_Module
IMPLICIT NONE
REAL(Kind=WP_Kind) a, rho, temp, tstore, umat
REAL(Kind=WP_Kind) ,ALLOCATABLE :: beta(:,:), idi(:,:), nid(:)
REAL(Kind=WP_Kind) ,ALLOCATABLE ::  plmst(:,:,:)
REAL(Kind=WP_Kind) ,ALLOCATABLE ::  djmkp(:,:,:,:)
REAL(Kind=WP_Kind) ,ALLOCATABLE ::  phi(:,:,:)
REAL(Kind=WP_Kind) ,ALLOCATABLE ::  betqf(:,:,:)
REAL(Kind=WP_Kind) ,ALLOCATABLE ::  chi1(:,:,:)
REAL(Kind=WP_Kind) ,ALLOCATABLE ::  thaph(:,:,:)
REAL(Kind=WP_Kind) ,ALLOCATABLE ::  snthd(:,:)
INTEGER         ,ALLOCATABLE ::  id(:,:)
REAL(Kind=WP_Kind) ,ALLOCATABLE ::  th(:,:), ch(:,:), f(:)
REAL(Kind=WP_Kind) ,ALLOCATABLE ::  h1(:,:)
INTEGER         ,ALLOCATABLE :: jx(:), jy(:), jc1(:), xnode(:), idni(:)
INTEGER naph, ndel, nxdim, nydim, jrjtmax, nsize, ierr, i, ithcll, ithsub, nidm, vib
LOGICAL little, medium, full
EXTERNAL AphDel_Old, AphDel_New, popt
INTRINSIC min0
!-----------------------------------------------------------------------
DIMENSION temp(*), tstore(*), umat(*)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Set the debug flags.
!-----------------------------------------------------------------------
DATA ithcll/0/, ithsub/0/
DATA little/.false./, medium/.false./, full/.false./
CALL popt('aphdelbf', little, medium, full, ithcll, ithsub)
WRITE(Out_Unit,*)'aph-to-delves projection'
  IF( .NOT. irdindep)THEN

!-----------------------------------------------------------------------
!     Find the dimensions for the arrays.
!-----------------------------------------------------------------------
   jrmax = 0
   DO karran = 1, narran
      DO vib=minvib(karran),maxvib(karran)
         IF(jrmax<jmax(vib,karran)) jrmax = jmax(vib,karran)
      ENDDO
   ENDDO

   WRITE(Out_Unit,*)'jrmax=',jrmax
   WRITE(Msg_Unit,*)'jrmax=',jrmax
   IF(mxl<jrmax+1)THEN
      WRITE(Out_Unit,*)'mxl too small: mxl,jrmax+1=',mxl,jrmax+1
      STOP 'aphdelbf'
   ENDIF

   nxdim=jrmax+1
   IF(mxglegn < jrmax+1)THEN
      WRITE(Out_Unit,*)' mxglegn too small. mxglegn,jrmax+1=', mxglegn, jrmax+1
      STOP 'aphdelbf'
   ENDIF

   IF(nxdim < 4)nxdim=4
   nydim=0
   DO i=1,narran
      IF(nhermt(i) > nydim)nydim=nhermt(i)
   ENDDO
   !revised by X.L.
   nydim=max(nydim,21)
   nsize=3*nxdim*nydim
   WRITE(Msg_Unit,*)'nxdim,nydim,nsize=',nxdim, nydim, nsize
   jrjtmax=min0(jtot, jrmax)


!-----------------------------------------------------------------------
!     Allocate memory for AphDel_Old or AphDel_New if this is the first energy.


WRITE(Out_Unit,*)'aph-to-delves projection'
ALLOCATE(plmst(0:jrjtmax,nxdim,0:jrmax))
ALLOCATE(djmkp(0:jrjtmax,nxdim,nydim,3))
ALLOCATE(phi(nxdim,nydim,narran))
ALLOCATE(betqf(nxdim,nydim,3))
ALLOCATE(chi1(nxdim,nydim,3))
ALLOCATE(thaph(nxdim,nydim,3))
ALLOCATE(snthd(nydim,3))
ALLOCATE(id(9,mxelmnt))
ALLOCATE(th(9,mxelmnt), ch(9,mxelmnt), f(mxnode))
ALLOCATE(h1(9,npoint))
ALLOCATE(jx(npoint), jy(npoint), jc1(npoint), xnode(mxnode), idni(mxnode))
ALLOCATE(beta(nidm,3*mxnode/4), idi(nidm,3*mxnode/4),nid(3*mxnode/4))


ENDIF
!-----------------------------------------------------------------------
!     Project the aph surface functions onto the delves vibrational
!     functions.
!-----------------------------------------------------------------------
DO i=1, narran
   IF(integrat(i))THEN
      c(i)=c(i)/rho
      b(i)=b(i)/rho
   ENDIF
ENDDO


CALL AphDel_New(rho, naph, ndel, tstore, umat, temp, plmst, djmkp, &
   phi, betqf, chi1, thaph, snthd, nxdim, nydim, jrmax, jrjtmax, id, &
   th, ch, f, h1, jx, jy, jc1, xnode, beta, idi, nid, idni, nidm)

DO i=1, narran
   IF(integrat(i))THEN
      c(i)=c(i)*rho
      b(i)=b(i)*rho
   ENDIF
ENDDO
!-----------------------------------------------------------------------
!     Deallocate memory for AphDel_Old or AphDel_New iff this is the first energy.
IF( .NOT. irdindep)THEN


DEALLOCATE(plmst)
DEALLOCATE(djmkp)
DEALLOCATE(phi)
DEALLOCATE(betqf)
DEALLOCATE(chi1)
DEALLOCATE(thaph)
DEALLOCATE(snthd)
DEALLOCATE(id)
DEALLOCATE(th, ch, f)
DEALLOCATE(h1)
DEALLOCATE(jx, jy, jc1, xnode, idni)
ENDIF


!-----------------------------------------------------------------------
RETURN
!------------------***END-aphdelbf***-----------------------------------
ENDSUBROUTINE aphdelbf
