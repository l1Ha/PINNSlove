SUBROUTINE upsiln (vecold,tstore,xsq,xk,rhoval,chinuj)
!
USE FileUnits_Module
USE Narran_Module
USE rotmax_Module
USE melm_Module
USE NumNuj_Module
USE new_Module
USE Qall_Module
USE GaussQuady_Module
USE opts_Module
USE pow_Module
USE chltot_Module
USE Spectro_Module
USE Arrch_Module
USE Storage_Module
USE Regins_Module
USE Masses_Module
USE TotalEng_Module
USE Oops_Module
USE region_Module
USE cherm1_Module
USE cherm2_Module
USE qnumsp_Module
USE rhonum_Module
USE NState_Module
USE Parms_module
USE DiatomicPot_MODULE, ONLY: DiatomicPot
USE Convrsns_Module, ONLY: autoev
IMPLICIT NONE
LOGICAL little,medium,full
INTEGER i, ii, ind1, ind2, index, iold, ithcll, ithcnl, ithsub, ip1, ip2, ip3
INTEGER jl, jlmax, jlp1, jold, mgh, nend, nosq, nstrt, nu, nudim, nuhigh
INTEGER nulow, nuold, nvp1, i1, i2, inew, innew, inold, jlmin, ms, nst
REAL(Kind=WP_Kind)al, rhoval, rmax, sqj, v2, xksq, xxk, expen, eint
REAL(Kind=WP_Kind) chinuj(nvbrthrt)
REAL(Kind=WP_Kind) xsq(nvibrot+500), xk(nvibrot+500)
REAL(Kind=WP_Kind) tstore(*), vecold(*)
!REAL(Kind=WP_Kind) tstore(nvibrot**2), vecold(nvibrot**2)
REAL(Kind=WP_Kind) esave(nstate+500)
REAL(Kind=WP_Kind) enrgy(nstate+500)
INTEGER,PARAMETER::NBasisFac=3

DATA ithcll/0/,ithsub/0/
DATA little/.false./,medium/.false./,full/.false./
CALL popt('upsiln  ',little,medium,full,ithcll,ithsub)
!WRITE(Out_Unit,*)'Entering Upsiln'
!-----------------------------------------------------------------------

IF(nrigid)RETURN
IF(.Not.Allocated(Numax))ALLOCATE(numax(nvibrot),numin(nvibrot))
jrmax=0
DO 160 karran=1,narran
IF(nchanl(1)==0.AND.karran==1) GOTO 160
IF(nchanl(2)==0.AND.karran==2) GOTO 160
IF(nchanl(3)==0.AND.karran==3) GOTO 160
IF(little) WRITE(Out_Unit,170) karran,rhoval
rhonum=rhoval
rmax=finish
IF(karran==1)THEN
   Ip1=2
   Ip2=3
   nstrt=1
   nend=nchanl(1)
ELSEIF(karran==2)THEN
   Ip1=3
   Ip2=1
   nstrt=nchanl(1)+1
   nend=nchanl(1)+nchanl(2)
ELSE
   Ip1=1
   Ip2=2
   nstrt=nchanl(1)+nchanl(2)+1
   nend=nchanl(1)+nchanl(2)+nchanl(3)
ENDIF
DiatomicPot=TRIM(AtomicSymbol(ip1))//TRIM(AtomicSymbol(ip2))
!-------------------------------------------------------------------
!  calculate alpha
!-------------------------------------------------------------------
IF(nrigid) alpha(karran)=1.0d+00
IF(.NOT.nrigid) alpha(karran)=ralpha(karran)*sqrt(we(karran)*usys)*re(karran)
IF(medium) WRITE(Out_Unit,240)
IF(medium) WRITE(Out_Unit,250) ralpha(karran),alpha(karran),rx(karran),re(karran),we(karran),wexe(karran),be(karran), de(karran),al
!-----------------------------------------------------------------------
! Determine the minimum/maximum little j used in this channel.
! Also find the maximum rotational quantum number in all channels.
!-----------------------------------------------------------------------
jlmax=0
jlmin=10000
DO ii=nstrt,nend
   IF(jrot3(ii)>jlmax) jlmax=jrot3(ii)
   IF(jrot3(ii)<jlmin) jlmin=jrot3(ii)
   IF(medium) WRITE(Out_Unit,*)'jrot3(', ii,')= ', jrot3(ii)
ENDDO
IF(jlmax>jrmax)jrmax=jlmax
IF(medium)WRITE(Out_Unit,*)'jrmax=',jrmax,' nstrt,nend=',nstrt,nend
!-----------------------------------------------------------------------
! store the maximum vibration for each little j in numax
! store the minimum vibration for each little j in numin
!-----------------------------------------------------------------------
jlmax=jlmax+1
jlmin=jlmin+1
DO jlp1=jlmin,jlmax
   jl=jlp1-1
   numax(jlp1)=0
   numin(jlp1)=10000
   DO i=nstrt,nend
      IF(jl==jrot3(i))THEN
         IF(mvib3(i)>numax(jlp1)) numax(jlp1)=mvib3(i)
         IF(mvib3(i)<numin(jlp1)) numin(jlp1)=mvib3(i)
         IF(medium) WRITE(Out_Unit,*)'mvib3(', i,')= ', mvib3(i)
      ENDIF
   ENDDO
ENDDO

ioops=.true.
IF(iregion=='jacobi ') ioops=.false.
rhoj=rhoval
IF(karran==1) index=1
IF(karran==2) index=1+nnuj(1)
IF(karran==3) index=1+nnuj(1)+nnuj(2)
DO jlp1=jlmin,jlmax
   IF(numin(jlp1)/=10000)THEN
      jl=jlp1-1
      nuhigh=numax(jlp1)
      nulow=numin(jlp1)
      nudim=nuhigh-nulow+1
      mgh=nhermt(karran)*nudim
      jrott=jl
      nosq=(NBasisFac*noscil(karran))*(NBasisFac*noscil(karran))
      DO i=1,nosq
         chlast(i,jl,karran)=chnow(i,jl,karran)
      ENDDO
      IF(.NOT.irdindep)THEN
         CALL VibFun_New(NBasisFac*noscil(karran),jl,nhermt(karran),tstore,chnow(1,jl,karran),chinuj(index),nulow,nuhigh,nudim,enrgy,   &
            alpha(karran),rx(karran),xpth(1,karran),wpth(1,karran), re(karran),xsq(nstrt),xk(nstrt),scheme)
      ENDIF
      IF(jlp1==jlmin)THEN
         WRITE(Out_Unit,*)
         WRITE(Out_Unit,*)"karran=",karran
         WRITE(Out_Unit,*)"DiatomicPot(upsiln)=", DiatomicPot
      ENDIF
      ithcnl=karran
      DO nvp1=1,nudim
         WRITE(Out_unit,'(2(A,i5),2(A,ES23.15))')"nu=",nvp1-1,"  jl=",jl,"  Delves_Energy(Ha)=",enrgy(nvp1),"  Delves_Energy(eV)=",enrgy(nvp1)*autoev
         nu=nvp1+nulow-1
         DO i=nstrt,nend
            IF(jl==jrot3(i).AND.nu==mvib3(i).AND.ithcnl==kchan(i))THEN
               esave(i)=enrgy(nvp1)
            ENDIF
         ENDDO
      ENDDO
!-----------------------------------------------------------------------
! store internal vibrational-rotational energies in energy.
!-----------------------------------------------------------------------
      index=index+mgh
   ENDIF
ENDDO
!-----------------------------------------------------------------------
! reorder the basis functions.
!-----------------------------------------------------------------------
IF(.NOT.irdindep)THEN
DO i=nstrt,nend
   eint=esave(i)
   xsq(i)=usys2*(Etot-eint)
   xksq3(i)=xsq(i)
   xk(i)=sqrt(ABS(xsq(i)))
ENDDO
IF(little)WRITE(Out_Unit,107)iregion,xksq3(1)
  107 FORMAT(1x,'for iregion=  ',a7,'  xksq3(1)=',f10.5)
ENDIF
IF(irdindep.or.iwrindep)THEN
  !Write(Out_unit,*)'Calling IO_Basis', karran, nstrt, nchanl(karran),'upsiln'
  !Write(Msg_unit,*)'Calling IO_Basis', karran, nstrt, nchanl(karran),'upsiln'
  CALL IO_Basis(nchanl(karran), xsq, xk, esave, nstrt-1)
ENDIF
IF(karran==1) ind1=1
IF(karran==2) ind1=1+nnuj(1)
IF(karran==3) ind1=1+nnuj(1)+nnuj(2)
IF(karran==1) ind2=nnuj(1)
IF(karran==2) ind2=nnuj(1)+nnuj(2)
IF(karran==3) ind2=nnuj(1)+nnuj(2)+nnuj(3)
index=0
DO i=ind1,ind2
   index=index+1
   vecold(index)=chinuj(i)
ENDDO
iold=0
DO 131 jlp1=jlmin,jlmax
jold=jlp1-1
nuhigh=numax(jlp1)
nulow=numin(jlp1)
IF(nulow==10001)GOTO 131
DO 130 nvp1=nulow+1,nuhigh+1
nuold=nvp1-1
iold=iold+1
DO 100 i=nstrt,nend
IF(jrot3(i)/=jold) GOTO 100
IF(mvib3(i)/=nuold) GOTO 100
inew=nuj3(i)
GOTO 110
  100 CONTINUE
WRITE(Out_Unit,180) jrot3(i),jold,mvib3(i),nuold
STOP 'upsiln'
  110 inold=nhermt(karran)*(iold-1)
innew=nhermt(karran)*(inew-1)
DO i=1,nhermt(karran)
   i1=innew+i
   i2=inold+i
   chinuj(i1)=vecold(i2)
ENDDO
  130 CONTINUE
  131  CONTINUE
!-----------------------------------------------------------------------
! IF desired print reordered unweighted wavefunctions.
!-----------------------------------------------------------------------
IF(full)THEN
WRITE(Out_Unit,190)
ms=0
IF(karran==1) nst=1
IF(karran==2) nst=1+nnuj(1)
IF(karran==3) nst=1+nnuj(1)+nnuj(2)
CALL MxOutD(chinuj(nst),nhermt(karran),iold,ms)
IF(melmnt=='space') WRITE(Out_Unit,210)
IF(melmnt=='probf') WRITE(Out_Unit,220)
IF(melmnt=='molbf') WRITE(Out_Unit,230)
DO 150 i=nstrt,nend
eint=esave(i)
xksq=usys2*(Etot-eint)
xxk=sqrt(ABS(xksq))
sqj=jrot3(i)*(jrot3(i)+1)
v2=mvib3(i)+0.5d+00
IF(scheme==1)THEN
   expen=we(karran)*v2 -wexe(karran)*v2**2 +weye(karran)*v2**3+be(karran)*sqj -de(karran)*sqj*sqj -alfe(karran)*v2*sqj+wzero(karran)
ELSEIF(scheme==2)THEN
   expen=0.0
ENDIF
IF(melmnt=='space') WRITE(Out_Unit,200) mvib3(i),jrot3(i),lorb3(i),eint,expen,xksq,xxk,nuj3(i)
IF(melmnt=='probf') WRITE(Out_Unit,200) mvib3(i),jrot3(i),mega3(i),eint,expen,xksq,xxk,nuj3(i)
IF(melmnt=='molbf') WRITE(Out_Unit,200) mvib3(i),lorb3(i),mega3(i),eint,expen,xksq,xxk,nuj3(i)
  150 CONTINUE
ENDIF
160 CONTINUE

!WRITE(Out_Unit,*)'Leaving Upsiln'
    !DEALLOCATE(numax,numin)
RETURN
!----------------***END-upsiln***------------------------------------
!
  170 FORMAT(//,'  delves surface function calculation for channel= ',i1,'   rho-j= ',f10.5)
  180 FORMAT(/,'***error***: msrot/=jold or msvib/=nuold',4i5)
  190 FORMAT(/,' unweighted vibrational eigenfunctions in upsiln')
  200 FORMAT(1x,3i5,3x,1pe14.7,3x,3es15.7,i6)
  210 FORMAT(/,'  mvib jrot lorb             eint        experiment  ',4x,'k-squared',14x,'k',1x,' nsuj in upsiln')
  220 FORMAT(/,'  mvib jrot mega             eint        experiment  ',4x,'k-squared',14x,'k',1x,' nsuj in upsiln')
  230 FORMAT(/,'  mvib lorb lamd             eint        experiment  ',4x,'k-squared',14x,'k',1x,' nsuj in upsiln')
  240 FORMAT(/,1x,8x,"ralpha",8x,"alpha",11x,"rx",10x,"re",11x,"we",10x,"wexe",10x,"be",9x,"de",10x,"alfe")
  250 FORMAT(/,10es13.6)
ENDSUBROUTINE Upsiln
