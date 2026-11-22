SUBROUTINE WavDer_Old(h,hd,hdn,z,l,noscil,scheme,alphas)
!
!---------------------------------------------------------------------
! routine to compute the noscil vibrational wavefunctions and
! their second derivatives at z bohr and for a fixed l
! written by b.j. archer
!---------------------------------------------------------------------
USE Numeric_Kinds_Module
USE fileunits_Module
USE Narran_Module
USE Arrch_Module
USE Masses_Module
USE etachn_Module
USE das_Module
IMPLICIT NONE
REAL(KIND=WP_Kind) alphas, hd, gamf, pref, rhof, z
INTEGER ithcll, ithsub, k, l, nn, noscil, scheme
LOGICAL little, medium, full
EXTERNAL hwave, dhep, popt
REAL(KIND=WP_Kind) h(noscil), hdn(noscil)

DATA ithcll/0/,ithsub/0/
DATA little/.false./,medium/.false./,full/.false./
CALL popt('WavDer_Old  ',little,medium,full,ithcll,ithsub)

IF(scheme==2)THEN

!---------------------------------------------------------------------
! make the Hylleras wavefunctions for k=0,noscil-1, returned in h
!---------------------------------------------------------------------
   CALL hwave(h,z,l,noscil)

!---------------------------------------------------------------------
! make the second derivative of the wavefunction
! hd is the part that does not depend on n, hdn is the part that does
! depend on n
!---------------------------------------------------------------------

   hd = -alphas*REAL(l*(l+1),WP_Kind)/(z*z)
   k=0
   gamf = 2.0d0*eta(karran)*usys/dscale(karran)
   rhof = gamf*z
   pref = alphas*gamf*gamf

   IF(little)THEN
      WRITE(Out_Unit,*)'z=',z,' karran=',karran,'eta=', eta(karran),'gamf=',gamf,'alphas=',alphas
   ENDIF

   DO nn=1,noscil
      hdn(nn) = pref*(REAL(k,WP_Kind)/(rhof*rhof) + REAL(k+l+1,WP_Kind)/rhof  - 0.25d0)
      IF(k>0)hdn(nn) = hdn(nn) - pref*sqrt(REAL(k*(k+2*l+2),WP_Kind))*h(nn-1)/(rhof*rhof*h(nn))
      k=k+1
   ENDDO

ELSE
!---------------------------------------------------------------------
! make the harmonic oscillator wavefunctions and second derivatives
!---------------------------------------------------------------------
   CALL dhep(h,z,noscil)
   hd = -alphas*z**2
   DO nn=1,noscil
      hdn(nn) = alphas*REAL(2*nn-1,WP_Kind)
   ENDDO

ENDIF

RETURN
ENDSUBROUTINE WavDer_Old
