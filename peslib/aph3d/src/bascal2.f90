SUBROUTINE bascal2(nang, thaph, chiaph, nbasis, phi, pn, dpndx, oscil, doscdz, parityfn, bfaci)
USE FileUnits_Module
USE Narran_Module
USE Parms_Module
USE Gaussb_Module
USE Quantb_Module
USE MassFactor2_Module
USE TotJ_Module
USE AzBzCz_Module
USE LBasis_Module
IMPLICIT NONE
INTEGER Nang, Nbasis, karran, iang, n, l, ibasis, marran, isp
REAL(KIND=WP_Kind) chi, schi, cchi, cthaph, sthaph, top, denom, xbig, product
REAL(KIND=WP_Kind) tanth, thlitdel, sinthd, costhd, sinths, cosths, AZSBZC, bfac, Zlit
REAL(KIND=WP_Kind), PARAMETER :: zero = 0.0D0, one = 1.0D0
REAL(KIND=WP_Kind), PARAMETER :: two = 2.0D0, four = 4.0D0
REAL(KIND=WP_Kind) thaph(nang), chiaph(nang), pn(0:mxl,nang,narran)
REAL(KIND=WP_Kind) dpndx(0:mxl,nang,narran)
REAL(KIND=WP_Kind) oscil(0:mxn,nang,narran), doscdz(0:mxn,nang,narran)
REAL(KIND=WP_Kind) parityfn(nang,narran)
REAL(KIND=WP_Kind) bfaci(nang,narran), phi(nang,nbasis), nfac(0:1000)
REAL(KIND=WP_Kind) fnorm(0:1000)
WRITE(Out_Unit,*)'nang=', nang
WRITE(Out_Unit,*)'nbasis=', nbasis
WRITE(Out_Unit,*)'mega=', mega
WRITE(Out_Unit,*)'Parity=', parity
DO karran = 1, narran
   WRITE(Out_Unit,*)
   WRITE(Out_Unit,*)'Arrangement Channel=', karran
   WRITE(Out_Unit,*)'chij(karran,1)=', chij(karran,1)
   WRITE(Out_Unit,*)'nlegndre(karran)', nlegndre(karran)
   WRITE(Out_Unit,*)'noscil(karran)', noscil(karran)
   WRITE(Out_Unit,*)'az(karran),bz(karran),cz(karran)', az(karran), bz(karran), &
                  & cz(karran)
 
   DO iang = 1, nang
       chi = chiaph(iang) - chij(karran,1)
       schi = sin(two*chi)
       cchi = cos(two*chi)
       cthaph = cos(thaph(iang))
       sthaph = sin(thaph(iang))
       top = sqrt(cthaph**2*cchi**2+schi**2)
!----------------------------------------------------------------------
! Calculate Legendre Basis
!----------------------------------------------------------------------
       IF(top/=0D0)THEN
           denom = one/top
           xbig = sthaph*schi*denom
           CALL legendrd(nlegndre(karran),mega,pn(0,iang,karran),               &
                       & dpndx(0,iang,karran),xbig,nfac)
       ELSE
           DO isp = 0, mxl
               pn(isp,iang,karran) = 0D0
           ENDDO
       ENDIF
       product = sthaph*cchi
       IF(product<zero)THEN
           tanth = (one-product)*denom
       ELSE
           tanth = top/(one+product)
       ENDIF
!-----------------------------------------------------------------------
! Calculate Harmonic Oscillator Basis
!-----------------------------------------------------------------------
       thlitdel = atan(tanth)
       sinthd = sin(thlitdel)
       costhd = cos(thlitdel)
       sinths = sinthd**2
       cosths = costhd**2
       azsbzc = one/(az(karran)*sinths+bz(karran)*cosths)
       bfac = four*sinths*cosths*sqrt(azsbzc)
       zlit = az(karran)*tanth - bz(karran)/tanth + cz(karran)
       CALL harmonid(noscil(karran),oscil(0,iang,karran),doscdz(0,iang,karran), &
                   & zlit,fnorm)
       IF(parity==1)THEN
           parityfn(iang,karran) = cos(chi)
       ELSE
           parityfn(iang,karran) = 1.D0
       ENDIF
       IF(bfac/=0D0)THEN
           bfaci(iang,karran) = one/bfac
       ELSE
           bfaci(iang,karran) = 0D0
       ENDIF
   ENDDO
ENDDO
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Basis functions being calculated'
WRITE(Out_Unit,*)'Number of Angles=', nang
WRITE(Out_Unit,*)'ibasis   chanl   nvib   jrot'
WRITE(Out_Unit,*)'jskip=', jskip
DO iang = 1, nang
   ibasis = 0
   DO marran = 1, narran
       DO n = minvib(marran), maxvib(marran)
           !jskip(marran) = 1
           DO l = jmin(n,marran,mega), jmax(n,marran), jskip(marran)
               ibasis = ibasis + 1
               IF(iang==1)WRITE(Out_Unit,*)ibasis, marran, n, l
               IF(ibasis>nbasis)THEN
                   WRITE(Msg_Unit,*)'Error ibasis>nbasis', ibasis, nbasis
                   STOP 'bascal'
               ENDIF
               phi(iang,ibasis) = pn(l,iang,marran)*oscil(n,iang,marran)        &
                                & *parityfn(iang,marran)*bfaci(iang,marran)
           ENDDO
       ENDDO
   ENDDO
ENDDO
WRITE(Out_Unit,*)'Number of Basis Functions = ', ibasis
ENDSUBROUTINE BASCAL2
