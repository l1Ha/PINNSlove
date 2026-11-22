SUBROUTINE qlevel (negvl, vect, haml, xsq, h, eint, nnuj, chinuj)

!-----------------------------------------------------------------------
! this routine determine the quantum numbers for a particular
! symmetry and the unweighted vibrational wavefunctions
! for that symmetry.
! written by g. a. parker. modified by jim lill.
!this routine is called by: input
!this routine calls: popt, figure, VibFun_Old, VibFun_New, eigfun
!-----------------------------------------------------------------------
USE Numeric_Kinds_Module
USE fileunits_Module
USE potzero_Module
USE das_Module
USE Arrch_Module
USE Masses_Module
USE Qnumber_Module
USE quantb_Module
USE new_Module
!USE Gaussb_Module
USE GaussQuady_Module !, ONLY: wpth, xpth
USE opts_Module
USE pow_Module
USE melm_Module
USE Spectro_Module
USE TotalEng_Module
USE totj_Module
USE approx_Module
USE elemnt_Module
USE qnumsp_Module
USE convrsns_Module
USE rhonum_MODULE
USE Regins_Module
USE Parms_Module
USE DiatomicPot_MODULE, ONLY: DiatomicPot
LOGICAL little, medium, full
INTEGER Negvl, i, ii, index, ithcll, ithsub, jl, jlmax, jlmin, jlp1, kthcall, mgh, mprnt
INTEGER nprnt, nudim, nuhigh, nulow, nvp1
INTEGER msvib(nstate), msrot(nstate), mslrb(nstate), msega(nstate), nsuj(nstate), mstype(nstate)
INTEGER nnuj(narran)
INTEGER ib, ib1, im1, imin, inew, iold, is2, jold, mt, neg1, nu, nuold
INTEGER i1, i2, icount, innew, inold, is, n, nsgvl, nulevl
INTEGER jrots, lorbs, meg, megas, megmax, megs, ms, mvibs
REAL(Kind=WP_Kind) emin, esave, expen, sqj, v2
REAL(Kind=WP_Kind) chinuj(nvbrthrt)
REAL(Kind=WP_Kind) eint(nvbrthrt)   !eint(nvibrot)
REAL(Kind=WP_Kind) energs(nstate+500)
REAL(Kind=WP_Kind) vect(nvbrthrt)   !vect(nvibrot**2)
REAL(Kind=WP_Kind) haml(nvbrthrt)   !haml(nvibrot**2)
REAL(Kind=WP_Kind) xsq(nvibrot+500), h(nvibrot+500)
INTEGER,PARAMETER::NBasisFac=3


DATA kthcall /1/
DATA ithcll/0/,ithsub/0/
DATA little/.false./,medium/.false./,full/.false./

CALL popt('qlevel  ',little, medium, full, ithcll, ithsub)
!little=.true.
!medium=.true.
!full=.true.
!!determine rho
rhonum=endaph
IF(little)WRITE(Out_Unit,380) jtot,karran
mprnt=minvib(karran)
nprnt=maxvib(karran)
IF(little)WRITE(Out_Unit,370) mprnt,nprnt
!-----------------------------------------------------------------------
! determine quantum numbers.
!-----------------------------------------------------------------------
CALL figure (jmin(0,karran,0),jmax(0,karran),mvib,jthrot,lorb,megaz,negvl,minvib(karran),maxvib(karran))
!-----------------------------------------------------------------------
! check to see if negvl is bigger than the dimensions
! statements allow.
!-----------------------------------------------------------------------
IF(negvl>nvibrot)THEN
   WRITE(Out_Unit,390) negvl
   STOP 'qlevel'
ENDIF

DO i=1,negvl
   msrot(i)=jthrot(i)
   mslrb(i)=lorb(i)
   msega(i)=megaz(i)
   msvib(i)=mvib(i)
ENDDO
nsgvl=negvl
!-----------------------------------------------------------------------
! determine the maximum little j used.
! determine the minimum little j used
!-----------------------------------------------------------------------
jlmax=0
jlmin=10000
DO ii=1,negvl
  IF(jthrot(ii)>jlmax) jlmax=jthrot(ii)
  IF(jthrot(ii)<jlmin) jlmin=jthrot(ii)
ENDDO
!-----------------------------------------------------------------------
! store the maximum vibration for each little j in numax
! store the minimum vibration for each little j in numin
! note numax and numin are scratch space, not used outside of qlevel
!-----------------------------------------------------------------------
jlmax=jlmax+1
jlmin=jlmin+1
DO jlp1=jlmin,jlmax
   jl=jlp1-1
   numax(jlp1)=0
   numin(jlp1)=10000
   DO i=1,negvl
      IF(jl==msrot(i))THEN
         IF(msvib(i)>numax(jlp1)) numax(jlp1)=msvib(i)
         IF(msvib(i)<numin(jlp1)) numin(jlp1)=msvib(i)
      ENDIF
   ENDDO
ENDDO
!-----------------------------------------------------------------------
! calculate wavefunctions and store them in chinuj
!-----------------------------------------------------------------------
index=1
DO jlp1=jlmin,jlmax
   IF(numin(jlp1)/=10000)THEN
      jl=jlp1-1
      nuhigh=numax(jlp1)
      nulow=numin(jlp1)
      IF(little)WRITE(Out_Unit,*)' nulow=', nulow
      nudim=nuhigh-nulow+1
      mgh=nhermt(karran)*nudim
      IF(nhermt(karran)<1)THEN
         WRITE(Out_Unit,*)'nhermt(karran)<1:nhermt,karran=',nhermt(karran),karran
         STOP 'QLevel'
      ENDIF

      CALL VibFun_New(NBasisFac*noscil(karran),jl,nhermt(karran),haml,vect, chinuj(index), &
         nulow,nuhigh,nudim,energs,alpha(karran),rx(karran),xpth(1,karran), &
         wpth(1,karran),re(karran),xsq,h,scheme)
      IF(ioops.eq.0)THEN
         WRITE(Out_Unit,*)
         WRITE(Out_Unit,*)"karran=",karran
         WRITE(Out_Unit,*)"DiatomicPot(qlevel)=", DiatomicPot
         DO nvp1=1,noscil(karran)
            WRITE(Out_unit,'(2(A,i5),2(A,ES23.15))')"nu=",nvp1-1,"  jl=",jl,"  Jacobi_Energy(Ha)=",energs(nvp1),"  Jacobi_Energy(eV)=",energs(nvp1)*autoev
         ENDDO
      ENDIF
      CALL EigFun(nudim,nhermt(karran),energs,chinuj(index),xpth(1,karran),  &
           wpth(1,karran),jl,karran,re(karran),rx(karran),alpha(karran),haml,h)
      !-----------------------------------------------------------------------
      ! store internal vibrational-rotational energies in energy.
      !-----------------------------------------------------------------------
      DO nvp1=1,nudim
          nu=nvp1+nulow-1
         DO i=1,negvl
            IF(jl==msrot(i).AND. nu==msvib(i))THEN
               energy(i)=energs(nvp1)
               !WRITE(Msg_Unit,*)i,energy(i)
            ENDIF
         ENDDO
      ENDDO
      index=index+mgh
   ENDIF
ENDDO
!-----------------------------------------------------------------------
! order quantum numbers according to energy.
!-----------------------------------------------------------------------
neg1=negvl-1
DO ib=1,neg1
   emin=energy(ib)
   imin=ib
   ib1=ib+1
   DO i=ib1,negvl
      IF(emin>energy(i))THEN
         emin=energy(i)
         imin=i
      ENDIF
   ENDDO
   esave=energy(ib)
   mvibs=msvib(ib)
   lorbs=mslrb(ib)
   jrots=msrot(ib)
   megas=msega(ib)
   energy(ib)=energy(imin)
   msvib(ib)=msvib(imin)
   mslrb(ib)=mslrb(imin)
   msrot(ib)=msrot(imin)
   msega(ib)=msega(imin)
   energy(imin)=esave
   msvib(imin)=mvibs
   msrot(imin)=jrots
   mslrb(imin)=lorbs
   msega(imin)=megas
ENDDO
!-----------------------------------------------------------------------
! WRITE out the complete basis set used.
!-----------------------------------------------------------------------
DO i=1,negvl
!-----------------------------------------------------------------------
! determine the symmetry type
!    mstype(I)=1 j-little even positive omegas.
!    mstype(I)=2 j-little odd positive omegas.
!    mstype(I)=3 j-little even negative omegas.
!    mstype(I)=4 j-little odd negative omegas.
!-----------------------------------------------------------------------
   IF(2*(msrot(i)/2)==msrot(i))THEN
      IF(msega(i)>=0)THEN
         mstype(i)=1
      ELSE
         mstype(i)=3
      ENDIF
   ELSE
      IF(msega(i)>=0)THEN
         mstype(i)=2
      ELSE
         mstype(i)=4
      ENDIF
   ENDIF
ENDDO
!-----------------------------------------------------------------------
! order quantum numbers according to omegas.
!-----------------------------------------------------------------------
megmax=0
DO i=1,negvl
   IF(iabs(msega(i))>megmax) megmax=iabs(msega(i))
ENDDO
megmax=megmax+1
icount=1
DO meg=1,megmax
   megs=meg-1
   is2=2
   IF(meg==1) is2=1
   DO is=1,is2
      IF(is==2) megs=-megs
      DO i=1,negvl
         IF(megs==msega(i))THEN
            eint(icount)=energy(i)
            mvib(icount)=msvib(i)
            jthrot(icount)=msrot(i)
            lorb(icount)=mslrb(i)
            megaz(icount)=msega(i)
            mtype(icount)=mstype(i)
            icount=icount+1
         ENDIF
      ENDDO
   ENDDO
ENDDO
!-----------------------------------------------------------------------
! order quantum numbers according to symmetry type.
!-----------------------------------------------------------------------
icount=1
DO mt=1,4
   DO i=1,negvl
      IF(mt==mtype(i))THEN
         energy(icount)=eint(i)
         msvib(icount)=mvib(i)
         msrot(icount)=jthrot(i)
         mslrb(icount)=lorb(i)
         msega(icount)=megaz(i)
         mstype(icount)=mtype(i)
         icount=icount+1
      ENDIF
   ENDDO
ENDDO
!-----------------------------------------------------------------------
! create a label for each nondegenerate basis function.
!-----------------------------------------------------------------------
icount=1
DO i=1,negvl
   DO ii=1,i-1
      IF(msvib(ii)==msvib(i).AND.msrot(ii)==msrot(i)) GOTO 190
   ENDDO
   nsuj(i)=icount
   icount=icount+1
   GOTO 210
  190 im1=i-1
   DO ii=1,im1
      IF(msvib(ii)==msvib(i).AND.msrot(ii)==msrot(i))THEN
         nsuj(i)=nsuj(ii)
      ENDIF
   ENDDO
  210 CONTINUE
ENDDO

DO mt=1,4
   icount=0
   DO i=1,negvl
      IF(mstype(i)==mt)icount=icount+1
   ENDDO
ENDDO
!-----------------------------------------------------------------------
! reorder the basis functions.
!-----------------------------------------------------------------------
index=index-1
DO i=1,index
   vect(i)=chinuj(i)
ENDDO
iold=0
DO jlp1=jlmin,jlmax
   jold=jlp1-1
   nuhigh=numax(jlp1)
   nulow=numin(jlp1)
   DO nvp1=nulow+1,nuhigh+1
      nuold=nvp1-1
      iold=iold+1
      DO 240 i=1,negvl
         IF(msrot(i)/=jold) GOTO 240
         IF(msvib(i)/=nuold) GOTO 240
         inew=nsuj(i)
         GOTO 250
  240 CONTINUE
      WRITE(Out_Unit,310) msrot(i),jold,msvib(i),nuold
      STOP 'qlevel'
  250 inold=nhermt(karran)*(iold-1)
      innew=nhermt(karran)*(inew-1)
      IF(inew==iold) GOTO 270
      DO i=1,nhermt(karran)
         i1=innew+i
         i2=inold+i
         chinuj(i1)=vect(i2)
      ENDDO
  270 CONTINUE
   ENDDO
ENDDO
nnuj(karran)=nhermt(karran)*iold
!   nnuj(karran) = iold
!-----------------------------------------------------------------------
! if desired print reordered unweighted wavefunctions.
!-----------------------------------------------------------------------
IF(full)THEN
   WRITE(Out_Unit,320)
   ms=0
   CALL MxOutD(chinuj,nhermt(karran),iold,ms)
ENDIF
!-----------------------------------------------------------------------
! store quantum numbers for a particular symmetry
! in jthrot, mvib and megaz and the associated energies
! in eint.
!-----------------------------------------------------------------------
negvl=0
IF((medium.AND.kthcall<=narran).or.full)THEN
IF(melmnt=='space') WRITE(Out_Unit,340)
IF(melmnt=='probf') WRITE(Out_Unit,350)
IF(melmnt=='molbf') WRITE(Out_Unit,360)
ENDIF
DO i=1,nsgvl
   negvl=negvl+1
   jthrot(negvl)=msrot(i)
   lorb(negvl)=mslrb(i)
   mvib(negvl)=msvib(i)
   megaz(negvl)=msega(i)
   lorb(negvl)=mslrb(i)
   eint(negvl)=energy(i)
   xsq(negvl)=usys2*(Etot-eint(negvl))
   h(negvl)=sqrt(ABS(xsq(negvl)))
   mtype(negvl)=mstype(i)
   nuj(negvl)=nsuj(i)
   sqj=jthrot(negvl)*(jthrot(negvl)+1)
   v2=mvib(negvl)+0.5d+00
   nulevl=nuj(negvl)
   DO ii=1,karran-1
     nulevl=nulevl+nnuj(ii)/nhermt(ii)
   ENDDO
   IF(scheme==1)THEN
      expen= we(karran)*v2 -wexe(karran)*v2**2 +weye(karran)*v2**3  &
            +be(karran)*sqj -de(karran)*sqj*sqj &
            -alfe(karran)*v2*sqj+wzero(karran)
   ELSEIF(scheme==2)THEN
      n= mvib(negvl) + jthrot(negvl) + 1
      expen = -usys/(dscale(karran)**2)/(2.d0*n**2)+vzero
   ENDIF
   IF((medium.AND.kthcall<=narran).or.full)THEN
      IF(melmnt=='space') WRITE(Out_Unit,330) mvib(negvl),jthrot(negvl),&
          lorb(negvl),eint(negvl),expen,xsq(negvl),h(negvl),nulevl,eint(negvl)*autoev
      IF(melmnt=='probf')WRITE(Out_Unit,330) mvib(negvl),jthrot(negvl), &
          megaz(negvl),eint(negvl),expen,xsq(negvl),h(negvl),nulevl,eint(negvl)*autoev
      IF(melmnt=='molbf')WRITE(Out_Unit,330) mvib(negvl),lorb(negvl), &
          megaz(negvl),eint(negvl),expen,xsq(negvl),h(negvl),nulevl,eint(negvl)*autoev
   ENDIF
ENDDO
kthcall = kthcall + 1
IF(negvl==0)THEN
   WRITE(Out_Unit,400)
   STOP 'qlevel'
ENDIF
RETURN

  310 FORMAT(//,'***error***: msrot/=jold or msvib/=nuold',4i5)
  320 FORMAT(//,' unweighted vibrational eigenfunctions in qlevel')
  330 FORMAT(1x,3i5,3x,1pe14.7,3x,3es15.7,i6,4x,es14.7)
  340 FORMAT(//,'  mvib jrot lorb             eint        experiment  ',4x,'k-squared',14x,'k',1x,' nsuj',10x,'eint(eV) in qlevel')
  350 FORMAT(//,'  mvib jrot megaz            eint        experiment  ',4x,'k-squared',14x,'k',1x,' nsuj',10x,'eint(eV) in qlevel')
  360 FORMAT(//,'  mvib lorb lamd             eint        experiment  ',4x,'k-squared',14x,'k',1x,' nsuj',10x,'eint(eV) in qlevel')
  370 FORMAT(//,"vibrational states are nu=",i5,3x,"to nu=",i5)
  380 FORMAT(//,"jtot=",i5,'  channel=  ',i1,'  jacobi region')
  390 FORMAT(//,"***** error nvibrot is too small   negvl=",i10)
  400 FORMAT(//,"***** error no basis possible")
ENDSUBROUTINE QLevel
