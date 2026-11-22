SUBROUTINE hwave(h,z,l,noscil)
!

! makes noscil Hylleras wavefn at the vibrational distance z and
! rotational angular mometnum l
! the particular Hylleras fns used are the 2l+2 fns of Shull & Lowdin,
! J.Chem.Phys. 23, 1362 (1955)
!-----------------------------------------------------------------------
! written by b.j. archer
!-----------------------------------------------------------------------
USE Numeric_Kinds_Module
USE fileunits_Module
USE Arrch_Module
USE Masses_Module
USE Narran_Module
USE etachn_Module
USE das_Module
USE opts_Module
IMPLICIT NONE
REAL(Kind=WP_Kind) alag, ffac, h, gamf, rfac, rhof, z
INTEGER k, l, lk, n, nn, noscil

EXTERNAL alagur, mfac
INTRINSIC exp, sqrt

DIMENSION alag(0:100), h(noscil)

! clear the wavefunctions

DO nn=1,noscil
   h(nn) = 0.0
ENDDO

lk=2*l+2
k=noscil-1
gamf = 2.0*eta(karran)*usys/dscale(karran)
rhof = gamf*z
! make the associated Laguerre polynomials for a fixed upper index lk
! and for a lower index k=0 to noscil-1
CALL alagur(lk,k,rhof,alag)
rfac=rhof**(lk/2)*exp(-0.5*rhof)
k=0

! make the Hylleras wavefunctions
DO 30 nn=1,noscil
n=k+l+1
! make the factorial(k) / fractorial(n+l+1)
CALL mfac(k,n+l+1,ffac)
h(nn) = alag(k)*sqrt(gamf*ffac)*rfac
k=k+1
30    CONTINUE

RETURN
!---------------***END-hwave***----------------------------------------
END
