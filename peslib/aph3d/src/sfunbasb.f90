SUBROUTINE sfunbasb(rho, nbasis, jtot, parity, tee, mega, marran, phi, ngleg, nherm, xptg, thetad, vect, nxdim, nydim, nbasiss)

!            P U R P O S E   O F    S U B R O U T I N E
!   This routine calculates the surface functions using an analytic
!   basis set.  This version just generates them for the aph to delves
!   transformation using coefficients determined by the abm version.


!            I N P U T    A R G U M E N T S
!  rho      current hyperradius.
!  jtot     total angular momentum.
!  parity   parity.
!  mega     projection of the total angular momentum on the body-fixed
!           z-axis.

!            O U T P U T    A R G U M E N T S
!  phi    surface functions at the quadrature points.
!  vect     matrix of coeffs of surface fcns expanded in prim. basis.
!
!this routine is called by:
!         intbas
!this routine calls
!         popt,basisb,diagbasb,mxoutd
!-----------------------------------------------------------------------
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
USE fileunits_Module
USE Numeric_Kinds_Module
USE Narran_Module
USE gaussb_Module
USE Numbers_Module
IMPLICIT NONE
SAVE
INTEGER jtot, parity, mega, nbasis, megas, marran, ibasis, jbasis, nherm, ngleg, nthcall, tee, nxdim, nydim, nbasiss
REAL(Kind=WP_Kind)rho,  xptg, thetad, pn, oscil, parityfn, vect, phi, bfaci

DIMENSION xptg(mxglegn, narran), nherm(narran), ngleg(narran), thetad(maxhermt, narran), pn(0:mxl, mxglegn, narran),&
         oscil(0:mxn, maxhermt, narran), vect(nbasiss, nbasiss), phi(nxdim, nydim, narran), parityfn(mxglegn, maxhermt, narran),&
         bfaci(maxhermt,narran)
DATA nthcall /0/, megas /-1/

!-----------------------------------------------------------------------
! Determine printing options.
!-----------------------------------------------------------------------
LOGICAL little, medium, full
INTEGER ithcall, ithsub
DATA ithcall /0/, ithsub /0/
DATA little /.true./, medium /.false./, full /.false./
CALL popt ('sfunbasb', little, medium, full, ithcall, ithsub)
!-----------------------------------------------------------------------
! Calculate the Legendre and Oscillator basis.
!-----------------------------------------------------------------------
IF(Medium)WRITE(Out_Unit,*)'begining of routine sfunbasb'

IF(nthcall==0.or.megas/=mega)THEN
   megas=mega
   nthcall=1
   CALL basisb (xptg, thetad, pn, ngleg, nherm, oscil, parity, parityfn, &
   bfaci, rho, mega)
ENDIF
!-----------------------------------------------------------------------
! Calculate the primitive basis and the surface functions at quad pts.
!-----------------------------------------------------------------------
CALL diagbasb (rho, nbasis, nherm, noscil, oscil, ngleg, nlegndre, pn, &
 parityfn, phi, tee, mega, vect, bfaci, nxdim, nydim,  nbasiss)
!-----------------------------------------------------------------------
!  WRITE out additional information if desired.
!-----------------------------------------------------------
IF(little)THEN
   WRITE(Out_Unit,*)'routine sfunbasb'
   WRITE(Out_Unit,*)'rho, jtot, parity, mega = ',rho,jtot,  parity,mega
ENDIF
IF(medium.AND..NOT.full)THEN
   WRITE(Out_Unit,*)'Part of vect'
   DO  ibasis = 1 , min(10,nbasiss)
      WRITE(Out_Unit,*)(vect(ibasis,jbasis),&
                    jbasis=1,min(10,nbasiss))
   ENDDO
ENDIF
IF(full)THEN
   WRITE(Out_Unit,*)'eigenvectors'
   CALL MxOutD(vect, nbasiss, nbasiss, 0)
ENDIF
IF(Medium)WRITE(Out_Unit,*)'END of routine sfunbasb'
RETURN
END
