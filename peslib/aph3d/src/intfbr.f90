SUBROUTINE intfbr (nxdim,nydim,nx,ny,chi1,thaph,phi,iaph,ichnl,rho)
!c  purpose:  calculates the fbr aph phi fcn at given set of thaph, ch
!c            from the dvr determined wavefunction.
!-----------------------------------------------------------------------
USE Narran_Module
USE Parms_Module
USE fileunits_Module
USE Numeric_Kinds_Module
USE Integrat_Module
USE todim_Module
USE cook_Module
USE var_Module
USE Numbers_Module
USE fbr_Module

IMPLICIT NONE
SAVE
REAL(Kind=WP_Kind)rho
INTEGER nxdim,nydim,nx,ny,iaph,ichnl
REAL(Kind=WP_Kind)chi1(nxdim,nydim,ichnl), thaph(nxdim,nydim,ichnl)
REAL(Kind=WP_Kind)phi(nxdim,nydim,ichnl)

REAL(Kind=WP_Kind), ALLOCATABLE :: plegn(:,:,:,:)
REAL(Kind=WP_Kind), ALLOCATABLE :: cos2n(:,:,:,:)
REAL(Kind=WP_Kind), ALLOCATABLE :: csurfbr(:),tthetat(:,:),&
                               pleg(:),falfa(:),fj(:)
INTEGER, ALLOCATABLE:: nev(:)

!c
INTEGER inevrd,isfun,ntheta,nchi,nthnch,ic,iy,ix,ii,jj,i,n,j,&
        nrhord,l,nthdvr,nchdvr,ksf,ialfa,nchith,ialf
REAL(Kind=WP_Kind)eps,sq2pii,sqpii,omeps,opeps,theta,costh,chi,&
                sumfalf,sumfj,sumpsi,rhord     !,psifbr
DATA inevrd/0/,isfun/0/,eps/1.0e-3/
!c
!     EXTERNAL plm
!c
!-----------------------------------------------------------------------
!c Determine printing options.
!-----------------------------------------------------------------------
LOGICAL little, medium, full
INTEGER ithsub, ithcall
DATA ithcall /0/, ithsub /0/
DATA little /.false./, medium /.false./, full /.false./
SAVE
CALL popt ('intfbr  ', little, medium, full, ithcall, ithsub)
!c

IF(little)THEN
  WRITE(Out_Unit,*) 'inside intfrb: rho = ',rho
  WRITE(Out_Unit,*) 'nxdim = ',nxdim  ,'  nx = ',nx
  WRITE(Out_Unit,*) 'nydim = ',nydim  ,'  ny = ',ny
  WRITE(Out_Unit,*) 'lmax  = ',lmax   ,'mmax = ',mmax
  WRITE(Out_Unit,*) 'nsfunc= ',nsfunc ,'nsym = ',nsym
ENDIF

!----------------------------------------------------------
!ALLOCATE momory
!----------------------------------------------------------
IF(IAPH==1)THEN
ALLOCATE( plegn(mxl,maxosc,nthetmax,narran))
ALLOCATE(cos2n(mxl,maxosc,nchimax,narran))
ALLOCATE(csurfbr(ndvptmax), tthetat(nthetmax,nthetmax),&
         pleg(nthetmax), nev(nthetmax), falfa(nthetmax),&
         fj(nthetmax))
ENDIF
!----------------------------------------------------------
!c
icook=0
iel=9
nnd5=5
kpri=0
!c
IF(nxdim > mxl .or. nydim > maxosc)THEN
  WRITE(Out_Unit,*) 'nxdim =',nxdim,' > ','mxl =',mxl
  WRITE(Out_Unit,*) '                OR'
  WRITE(Out_Unit,*) 'nydim =',nydim,' > ','maxosc =',maxosc
  STOP 'intfbr'
ENDIF
!c
ntheta = lmax + 1
nchi = mmax
nthnch = ntheta*nchi
!c
!c check array sizes for surface functions
!c
IF(ntheta>nthetmax)THEN
  WRITE(Out_Unit,*)' lmax + 1 (.fbrdata) > nthetmax (dvr.parms) '
  WRITE(Out_Unit,*)' make nthetmax bigger '
  STOP 'intfbr'
ENDIF
IF(nchi>nchimax)THEN
  WRITE(Out_Unit,*)' mmax (fbrdata) > nchimax (dvr.parms) '
  WRITE(Out_Unit,*)' make nchimax bigger '
  STOP 'intfbr'
ENDIF
!c
DO ic=1,ichnl
IF(integrat(ic))THEN
DO iy=1,ny
DO ix=1,nx
!c
IF(inevrd == 0)THEN
  OPEN(Unit=sfunfbr_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/SfunFBR.bin',form='unformatted',status='old')
  sqpii = 1.0/sqrt(pi)
  sq2pii = sqrt(2.0/pi)
ENDIF
!c
!c
IF(inevrd == 0)THEN
inevrd = 1
REWIND sfunfbr_unit
!c
    4 READ(sfunfbr_unit,END=991) nrhord,rhord,nthdvr,nchdvr
WRITE(Out_Unit,*)'reading surfaces functions in ',&
                'sfunfbr at rho = ',rhord
!c
!c check to make sure .fbrdata and .dvrdat are consistent
!c
IF(ntheta/=nthdvr)THEN
  WRITE(Out_Unit,*)' lmax(.fbrdata) should match ',&
                'nthesec(.dvrdat) ',ntheta,nthdvr
  STOP 'intfbr'
ENDIF
IF(nchi/=nchdvr)THEN
  WRITE(Out_Unit,*) 'mmax(.fbdata) should match nchisec(.dvrdat) '
  STOP 'intfbr'
ENDIF
!c
READ(sfunfbr_unit) (nev(ii),ii=1,nthdvr)
DO ii=1,nthdvr
  READ(sfunfbr_unit) (tthetat(ii,jj),jj=1,nthdvr)
ENDDO
!c
omeps = 1.0e0 - eps
opeps = 1.0e0 + eps
IF(omeps*rho > rhord .or.&
   opeps*rho < rhord )THEN
!c READ the unwanted fbr coefficients
  DO ksf=1,nsfunc
    READ(sfunfbr_unit) (csurfbr(i),i=1,nthnch)
  ENDDO
GOTO 4
ENDIF
ENDIF
!c
IF(iaph == 1)THEN
  theta = 2.0*thaph(ix,iy,ic)
  costh = cos(theta)
  l = ntheta - 1
  CALL legenfbr (l,pleg,costh)
  DO n=1,ntheta
    plegn(ix,iy,n,ic) = pleg(n)
  ENDDO

!*       DO 10 n=1,ntheta
!*         l = n - 1
!*         pleg(n) = plm_fun(l,0,costh)
!*  10   CONTINUE
!c
  chi = chi1(ix,iy,ic)
!c
  IF(nsym == 1)THEN
!c a1 symmetry
    cos2n(ix,iy,1,ic) = sqpii
    DO n=2,nchi
      cos2n(ix,iy,n,ic) = sq2pii*cos((n + n - 2)*chi)
    ENDDO
  ELSEIF (nsym == 2)THEN
!c a2 symmetry
    DO n=1,nchi
      cos2n(ix,iy,n,ic) = sq2pii*sin((n + n - 1)*chi)
    ENDDO
  ELSEIF (nsym == 3)THEN
!c b1 symmetry
    DO n=1,nchi
      cos2n(ix,iy,n,ic) = sq2pii*sin((n + n)*chi)
    ENDDO
  ELSEIF (nsym == 4)THEN
!c b2 symmetry
    DO n=1,nchi
      cos2n(ix,iy,n,ic) = sq2pii*cos((n + n - 1)*chi)
    ENDDO
  ELSE
     WRITE(Out_Unit,*)'Incorrect nsym: nsym=',nsym
     STOP 'IntDVR 6'
  ENDIF
!c
ENDIF
!------------------------------------------------------------------
IF(iaph /= isfun)THEN
  isfun = iaph
  READ(sfunfbr_unit) (csurfbr(i),i=1,nthnch)
ENDIF
!------------------------------------------------------------------
nchith = 0
DO ialfa=1,ntheta
!c
  IF(nev(ialfa) == 0)THEN
    falfa(ialfa) = 0.0
    GOTO 95
  ENDIF
!c
  sumfalf = 0.0
  DO n=1,nchi
    sumfalf = sumfalf + csurfbr(nchith + n)*cos2n(ix,iy,n,ic)
  ENDDO
!c
  falfa(ialfa) = sumfalf
!c
   95   CONTINUE
  nchith = nchith + nchi
ENDDO

!------------------------------------------------------------------
DO j=1,ntheta
  sumfj = 0.0
  DO ialf=1,ntheta
    sumfj = sumfj + tthetat(j,ialf)*falfa(ialf)
  ENDDO
  fj(j) = sumfj
ENDDO
!------------------------------------------------------------------
sumpsi = 0.0
DO j=1,ntheta
  sumpsi = sumpsi + fj(j)*plegn(ix,iy,j,ic)
ENDDO
 phi(ix,iy,ic)=sumpsi
!      psifbr = sumpsi
!*     WRITE(Out_Unit,160) theta,chi,psifbr
!* 160 FORMAT(2x,e12.5,2x,e12.5,2x,e15.8)
!------------------------------------------------------------------
!      phi(ix,iy,ic)=psifbr
 ENDDO
 ENDDO
 ENDIF
 ENDDO
!----------------------------------------------------------
!DEALLOCATE momory
!----------------------------------------------------------
IF(iaph==80)THEN
WRITE(Msg_Unit,*)'watch out in intfbr naph=80!!!!'
DEALLOCATE( plegn)
DEALLOCATE(cos2n)
DEALLOCATE(csurfbr, tthetat, pleg, nev, falfa, fj)
ENDIF
!----------------------------------------------------------


RETURN
!c
  991 WRITE(Out_Unit,*) 'READ off the END of the fbr coefficient file'
STOP 'intfbr'
ENDSUBROUTINE IntFBR
