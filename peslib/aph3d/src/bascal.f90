SUBROUTINE bascal (thaph, chiaph, nbasis, neigmin, phi, nbasiss, nang)
USE FileUnits_Module
USE Parms_Module
USE Narran_Module
USE Gaussb_Module
USE Quantb_Module
USE MassFactor2_Module
USE Lbasis_Module
USE TotJ_Module
USE AzBzCz_Module
USE Numbers_Module
IMPLICIT NONE
INTEGER nbasis, neigmin, nbasiss, jpar, kbasis, ibasis, karran, nang
INTEGER Marran, iang, isp, n, l
REAL(Kind=WP_Kind) thaph(nang),chiaph(nang),pn(0:mxl,nang,narran), dpndx(0:mxl,nang,narran)
REAL(Kind=WP_Kind) oscil(0:mxn,nang,narran), doscdz(0:mxn,nang,narran)
REAL(Kind=WP_Kind) parityfn(nang,narran), bfaci(nang,narran), phi(nang,mxbasis)
REAL(Kind=WP_Kind) NFAC(1000), fnorm(1000)
REAL(Kind=WP_Kind) chi, schi, cchi, cthaph, sthaph, top, denom, xbig, product, tanth, thlitdel
REAL(Kind=WP_Kind) sinthd, costhd, sinths, cosths, azsbzc, bfac, zlit, sfac
LOGICAL lfirst
DATA lfirst/.true./

!nang=mxang
WRITE(Out_Unit,*)'bascal called'
WRITE(Out_Unit,*)'nbasis,neigmin,nbasiss,nang=',nbasis,neigmin,nbasiss,nang
IF(nang>mxang)THEN
   WRITE(Out_Unit,*)'Error: nang>mxang',nang,mxang
   STOP 'bascal'
ENDIF

IF(lfirst)THEN
   WRITE(Out_Unit,*)'nang=',nang
   WRITE(Out_Unit,*)'nbasis=',nbasis, '   nbasiss=',nbasiss
   WRITE(Out_Unit,*)'mega=',mega
   WRITE(Out_Unit,*)'Parity=',parity
   WRITE(Out_Unit,*)'symmetry=',symmetry
   WRITE(Out_Unit,*)'jeven=',jeven
ENDIF
DO karran = 1, narran
   IF(lfirst)THEN
      WRITE(Out_Unit,*)
      WRITE(Out_Unit,*)'Arrangement Channel=',karran
      WRITE(Out_Unit,*)'chij(karran,1)=',chij(karran,1)
      WRITE(Out_Unit,*)'nlegndre(karran)',nlegndre(karran)
      WRITE(Out_Unit,*)'noscil(karran)',noscil(karran)
      WRITE(Out_Unit,*)'az(karran),bz(karran),cz(karran)', az(karran),bz(karran),cz(karran)
   ENDIF
   DO iang = 1, nang
      chi=chiaph(iang)-chij(karran,1)
      schi   = sin(two*chi)
      cchi   = cos(two*chi)
      cthaph=cos(thaph(iang))
      sthaph=sin(thaph(iang))
      top=sqrt(cthaph**2*cchi**2 + schi**2)
!----------------------------------------------------------------------
! Calculate Legendre Basis
!----------------------------------------------------------------------
      IF(top/=0d0)THEN
         denom=one/top
         xbig =  sthaph*schi*denom
         CALL legendrd (nlegndre(karran), mega, pn(0, iang, karran), dpndx(0,iang,karran), xbig, nfac)
      ELSE
         DO isp=0,mxl
            pn(isp,iang,karran)=0d0
         ENDDO
      ENDIF

      product=sthaph*cchi
      IF(product<zero)THEN
         tanth = (one-product)*denom
      ELSE
         tanth = top/(one+product)
      ENDIF

!-----------------------------------------------------------------------
! Calculate Harmonic Oscillator Basis
!-----------------------------------------------------------------------
      thlitdel=atan(tanth)
      sinthd=sin(thlitdel)
      costhd=cos(thlitdel)
      sinths=sinthd**2
      cosths=costhd**2
      azsbzc=one/(az(karran)*sinths+bz(karran)*cosths)
      bfac=four*sinths*cosths*sqrt(azsbzc)
      zlit=az(karran)*tanth-bz(karran)/tanth+cz(karran)
      CALL harmonid (noscil(karran), oscil(0, iang, karran), doscdz(0,iang,karran), zlit, fnorm)
      IF(parity==1)THEN
         parityfn(iang,karran)=cos(chi)
      ELSE
         parityfn(iang,karran) = 1.d0
      ENDIF

      IF(bfac/=0d0)THEN
         bfaci(iang,karran)=one/bfac
      ELSE
         bfaci(iang,karran)=0d0
      ENDIF
   ENDDO
ENDDO
IF(lfirst)THEN
   WRITE(Out_Unit,*)
   WRITE(Out_Unit,*)'Basis functions being calculated'
   WRITE(Out_Unit,*)'Number of Angles=',nang
   WRITE(Out_Unit,*)'ibasis   chanl   nvib   jrot'
   WRITE(Out_Unit,*)'jskip=',jskip
ENDIF
DO iang = 1, nang
   ibasis = 0
   DO marran = 1, narran
      DO n = minvib(marran), maxvib(marran)
         DO l = jmin(n,marran,mega), jmax(n,marran), jskip(marran)
            ibasis = ibasis + 1
            IF(iang==1.AND.lfirst)THEN
               WRITE(Out_Unit,*)ibasis,marran, n, l
            ENDIF

            IF(ibasis>nbasis)THEN
              WRITE(Msg_Unit,*)'Error ibasis>nbasis',ibasis,nbasis
              STOP 'bascal'
            ENDIF

            phi(iang,ibasis) = pn(l, iang, marran)*oscil(n, iang, marran)*parityfn(iang,marran)*bfaci(iang,marran)
         ENDDO
      ENDDO
   ENDDO
ENDDO
IF(symmetry)THEN
   IF(jeven)THEN
      jpar = 0
   ELSE
      jpar = 1
   ENDIF
   DO ibasis = 1, nbasiss
      IF(chanl(ibasis, mega)==2)THEN
         kbasis = ibasis + nbasis - nbasiss
         IF(mod(jrot(ibasis, mega),2)==jpar)THEN
            sfac = 1.d0
         ELSE
            sfac = -1.d0
         ENDIF
         DO iang = 1, nang
            phi(iang, ibasis) = (phi(iang, ibasis) + sfac*phi(iang, kbasis)) / sqrt(two)
         ENDDO
      ENDIF
      DO iang = 1, nang
         phi(iang, ibasis) = phi(iang, ibasis)*sqrt(two)
      ENDDO
   ENDDO
   ibasis = nbasiss
ENDIF

IF(lfirst)THEN
   WRITE(Out_Unit,*)'Number of Basis Functions = ',ibasis
ENDIF

lfirst = .false.

RETURN
ENDSUBROUTINE bascal
