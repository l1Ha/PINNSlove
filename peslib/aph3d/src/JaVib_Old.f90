SUBROUTINE JaVib_Old(yg,noscil,h,jl,expfac,kin,chi,dchi,ichan,cg)
USE FileUnits_Module
USE nstate_Module
USE Narran_Module
USE parms_Module
USE cherm2_Module
USE etachn_Module
USE das_Module
USE Masses_Module
USE opts_Module
!

! make the jacobi vibrational wavefuncition, chi, and its first
! derivative, dchi, at the delves point yg
!-----------------------------------------------------------------------
! written by b.j. archer
!-----------------------------------------------------------------------

IMPLICIT NONE
INTEGER ichan, jl, k, kin, kv, lk, n, noscil
REAL(KIND=WP_Kind) cg, chi, dchi, dhx, gamf, sf, yg
EXTERNAL alagur, dhep, hwave
intrinsic sqrt
REAL(Kind=WP_Kind) h(100,3), expfac(3), alag(0:100)

IF(scheme==2)THEN

! check cg to make sure the functions are computed at the correct vibrational distance
IF(ABS(1.0d0-cg)>1.5d-8)THEN
   WRITE(Out_Unit,*)'cg/=1.0 => yg/=sf',cg
   STOP 'JaVib_Old'
ENDIF
gamf=2.0d0*eta(ichan)*usys/dscale(ichan)

   call hwave(h(1,ichan),yg,jl,noscil)         ! get the noscil Hylleras functions

lk = 2*jl + 2
sf = yg*gamf
k=noscil-1

! get the k associated Laguerre polynomials
CALL alagur(lk,k,sf,alag)

k=0
chi=0.0d0
dchi=0.0d0
DO kv=1,noscil
   n=k+jl+1
! make the first derivative of the nth Hylleras function
   dhx = n/yg - gamf*0.5d0
   IF(k>0)dhx = dhx - (n+jl+1)*alag(k-1)/(yg*alag(k))
   dhx = dhx*h(kv,ichan)
! make the derivative of the kv wavefunction
   dchi = dchi+chnow(kin+kv,jl,ichan)*dhx
! make the kv wavefunction
   chi = chi+chnow(kin+kv,jl,ichan)*h(kv,ichan)
   k=k+1
ENDDO
chi=chi*expfac(ichan)
dchi=dchi*expfac(ichan)

ELSEIF(scheme==1)THEN
! make the harmonic oscillator basis for chi and dchi
CALL dhep (h(1,ichan),yg,noscil)
chi=0.0d0
! make the kv wavefunction
DO kv=1,noscil
   chi=chi+chnow(kin+kv,jl,ichan)*h(kv,ichan)
ENDDO
! make the first derivative of the kv wavefunction
dchi=0.0d0
DO kv=2,noscil
   dchi=dchi+chnow(kin+kv,jl,ichan)*h(kv-1,ichan)*sqrt(2.d0*(kv-1))
ENDDO
dchi=(dchi-chi*yg)/cg
chi=chi*expfac(ichan)
dchi=dchi*expfac(ichan)
ENDIF

RETURN
ENDSUBROUTINE JaVib_Old
