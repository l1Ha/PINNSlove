SUBROUTINE intdvr (nang, chi1, thaph, phi, naph, rho)
USE Numeric_Kinds_Module
USE FileUnits_Module
USE DVR_Module
USE Parms_Module
USE Numbers_Module
USE FBR_Module
!-----------------------------------------------------------------------
!  constructs the fbr(dvr) aph surface fcns (phi) at given set of
!  angles thaph, chi1.
!  the angles may be those of either the dvr or the abm quadrature.
!-----------------------------------------------------------------------
INTEGER nev(nthetmax), i, ialfa, iang,iaph, ii, inevrd, isfun, jj, ksf
INTEGER l, n, nang, naph, nchdvr, nchi, nchith, nrhord, nthdvr, ntheta, nthnch
INTEGER ialf, j
REAL(Kind=WP_Kind), INTENT(IN) :: chi1(nang), thaph(nang)
REAL(Kind=WP_Kind), INTENT(OUT):: phi(nang,naph)
REAL(Kind=WP_Kind),ALLOCATABLE :: plegn(:,:)       !plegn(nang,ntheta)
REAL(Kind=WP_Kind),ALLOCATABLE :: cos2n(:,:)       !cos2n(nang,nchi)
REAL(Kind=WP_Kind),ALLOCATABLE :: csurfbr(:)       !csurfbr(nang)
REAL(Kind=WP_Kind),ALLOCATABLE :: tthetat(:,:)     !tthetat(ntheta,ntheta)
REAL(Kind=WP_Kind),ALLOCATABLE :: pleg(:)          !pleg(ntheta)
REAL(Kind=WP_Kind),ALLOCATABLE :: falfa(:), fj(:)  !falfa(ntheta), fj(ntheta)
REAL(Kind=WP_Kind) chi, costh, eps, omeps, opeps, rho, rhord, sq2pii, sqpii, theta
REAL(Kind=WP_Kind) sumfalf, sumfj, sumpsi, psifbr
DATA inevrd/0/,isfun/0/,eps/1.0e-3/
!
!-----------------------------------------------------------------------
! Determine printing options.
!-----------------------------------------------------------------------
LOGICAL little, medium, full, there
INTEGER ithsub, ithcall
DATA ithcall /0/, ithsub /0/
DATA little /.false./, medium /.false./, full /.false./
CALL popt ('intdvr  ', little, medium, full, ithcall, ithsub)
!
WRITE(Out_Unit,*)'intdvr called:'
WRITE(Out_Unit,*)'nang,naph,rho=',nang,naph,rho

IF(little)THEN
  WRITE(Out_Unit,*) 'inside intdvr: rho = ',rho
  WRITE(Out_Unit,*) 'nang = ',nang
  WRITE(Out_Unit,*) 'lmax  = ',lmax   ,'mmax = ',mmax
  WRITE(Out_Unit,*) 'nsfunc= ',nsfunc ,'nsym = ',nsym
ENDIF
!
ntheta = lmax + 1
nchi   = mmax
nthnch = ntheta*nchi


ALLOCATE(plegn(nang,ntheta))
ALLOCATE(cos2n(nang,nchi))
ALLOCATE(csurfbr(MAX(nang,nthnch)))
ALLOCATE(tthetat(ntheta,ntheta))
ALLOCATE(pleg(0:ntheta))
ALLOCATE(falfa(ntheta))
ALLOCATE(fj(ntheta))

DO iaph=1,naph
   DO iang=1,nang
      IF(inevrd == 0)THEN
         INQUIRE (file=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Sfunfbr_DVR.bin', exist = there)
         IF(there)THEN
            OPEN(Unit=SFunFBR_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Sfunfbr_DVR.bin', form='unformatted',status='old')
         ELSE
            WRITE(Msg_Unit,*)'Error Sfunfbr_DVR.bin does not exist.'
            WRITE(Out_Unit,*)'Error Sfunfbr_DVR.bin does not exist.'
            STOP 'Stopping in intdvr 1'
         ENDIF
         sqpii = 1.0/sqrt(pi)
         sq2pii = sqrt(2.0/pi)
      ENDIF
      IF(inevrd == 0)THEN
         inevrd = 1
         REWIND SFunFBR_Unit
       4 READ(SFunFBR_Unit,END=991)nrhord,rhord,nthdvr,nchdvr
         WRITE(Out_Unit,*)'reading surfaces functions in sfunfbr at rhord = ',rhord
         WRITE(Out_Unit,*)"NRhord=",nrhord," NthDVR=",nthdvr," NchDVR=",nchdvr
         IF(ntheta/=nthdvr)THEN   ! check to make sure the data is consistent
            WRITE(Out_Unit,*)' lmax(.fbrdata) should match nthesec(.dvrdat) '
            WRITE(Out_Unit,*)'ntheta,nthdvr=',ntheta,nthdvr
            WRITE(*,*)' lmax(.fbrdata) should match nthesec(.dvrdat) '
            WRITE(*,*)'ntheta,nthdvr=',ntheta,nthdvr
            STOP 'Stopping in intdvr 2'
         ENDIF
         IF(nchi/=nchdvr)THEN     ! check to make sure the data is consistent
            WRITE(Out_Unit,*) 'mmax(.fbrdata) should match nchisec(.dvrdat) '
            WRITE(Out_Unit,*)'nchi,nchdvr=',nchi,nchdvr
            WRITE(*,*) 'mmax(.fbrdata) should match nchisec(.dvrdat) '
            WRITE(*,*)'nchi,nchdvr=',nchi,nchdvr
            STOP 'Stopping in intdvr 3'
         ENDIF
         READ(SFunFBR_Unit) (nev(ii),ii=1,nthdvr)
         WRITE(Out_Unit,*)'In IntDVR nev'
         WRITE(Out_Unit,'(10I5)')(nev(ii),ii=1,nthdvr)
         DO ii=1,nthdvr
            READ(SFunFBR_Unit) (tthetat(ii,jj),jj=1,nthdvr)
         ENDDO
         omeps=1.d0-eps
         opeps=1.d0+eps
         IF(omeps*rho > rhord .or. opeps*rho < rhord )THEN
            DO ksf=1,nsfunc
               READ(SFunFBR_Unit)   ! READ the unwanted fbr coefficients
            ENDDO
            GOTO 4
         ENDIF
      ENDIF
      IF(iaph == 1)THEN
         theta = 2.0*thaph(iang)
         costh = cos(theta)
         l = ntheta - 1
         CALL legenfbr (l,pleg,costh)
         DO n=1,ntheta
            plegn(iang,n) = pleg(n-1)   ! legendre polynomials of 2*ThetAPH
         ENDDO
         chi = chi1(iang)
         IF(nsym == 1)THEN             ! a1 chi symmetry
            cos2n(iang,1) = sqpii
            DO n=2,nchi
               cos2n(iang,n) = sq2pii*cos((n + n - 2)*chi)
            ENDDO
         ELSEIF (nsym == 2)THEN         ! a2 chi symmetry
            DO n=1,nchi
               cos2n(iang,n) = sq2pii*sin((n + n - 1)*chi)
            ENDDO
         ELSEIF (nsym == 3)THEN         ! b1 chi symmetry
            DO n=1,nchi
               cos2n(iang,n) = sq2pii*sin((n + n)*chi)
            ENDDO
         ELSEIF (nsym == 4)THEN         ! b2 chi symmetry
            DO n=1,nchi
               cos2n(iang,n) = sq2pii*cos((n + n - 1)*chi)
            ENDDO
         ELSE
            WRITE(Out_Unit,*)'Incorrect nsym: nsym=',nsym
            STOP 'Stopping in IntDVR 4'
         ENDIF
      ENDIF
      IF(iaph /= isfun)THEN
         isfun = iaph
         READ(SFunFBR_Unit) (csurfbr(i),i=1,nthnch) !  READ fbr coefficients
      ENDIF
      nchith = 0              ! Staring to construct fbr wavefunctions
      DO ialfa=1,ntheta
         IF(nev(ialfa) == 0)THEN
            falfa(ialfa) = 0.0
         ELSE
            sumfalf = 0.0
            DO n=1,nchi
               sumfalf = sumfalf + csurfbr(nchith + n)*cos2n(iang,n)
            ENDDO
            falfa(ialfa) = sumfalf
         ENDIF
         nchith = nchith + nchi
      ENDDO
      DO j=1,ntheta
         sumfj = 0.0
         DO ialf=1,ntheta
            sumfj = sumfj + tthetat(j,ialf)*falfa(ialf)
         ENDDO
         fj(j) = sumfj
      ENDDO
      sumpsi = 0.0
      DO j=1,ntheta
         sumpsi = sumpsi + fj(j)*plegn(iang,j)
      ENDDO
      psifbr = sumpsi
      phi(iang,iaph)=psifbr   ! Construct fbr wavefunction
   ENDDO
ENDDO

DEALLOCATE(plegn)
DEALLOCATE(cos2n)
DEALLOCATE(csurfbr)
DEALLOCATE(tthetat)
DEALLOCATE(pleg)
DEALLOCATE(falfa, fj)

RETURN
991 WRITE(Out_Unit,*) 'READ off the END of the fbr coefficient file'
STOP 'Stopping in intdvr 5'
END
