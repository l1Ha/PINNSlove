SUBROUTINE AphDel_Old(rho, naph, ndelve, tstore, umat, temp,&
      plmst,djmkp,phi,betqf,chi1,thaph,snthd,&
      nxdim,nydim,jrmax,jrjtmax,id,th,ch,f,&
      h1,jx,jy,jc1,xnode,beta,idi,nid,idni,nidm)
!-----------------------------------------------------------------------
!this routine is called by:
!         aphdelbf
!this routine calls
!         popt,numgen,factor,lgnd11,aphdef,legendre,clebp,intfbr,intbas,wigner,mxoutl,delchk,aphchk
!-----------------------------------------------------------------------
USE fileunits_Module
USE Numeric_Kinds_Module
USE Convrsns_Module
USE Narran_Module
USE Integrat_Module
USE sfntyp_Module
USE totj_Module
USE Numbers_Module
USE thetas_Module
USE qlam_Module
USE Qall_Module
!USE Gaussb_Module
USE GaussQuady_Module
USE opts_Module
USE chiang_Module
USE VFunc_Module
USE Storage_Module
USE rhosur_Module
USE TotalEng_Module
USE Masses_Module
USE quantb_Module
IMPLICIT NONE
SAVE
LOGICAL AphDel_All
LOGICAL little, medium, full
INTEGER ithcll, ithsub
INTEGER tea, lam, karran
INTEGER lamsave, jtsav,ipsav
INTEGER naph, ndelve, jrmax, nxdim, nydim, nidm, jrjtmax, indexdel
INTEGER jx(npoint), jy(npoint), jc1(npoint), xnode(mxnode), idni(mxnode)
INTEGER id(9,mxelmnt),ngor(narran)
INTEGER nx,k,i,ix,in,iy,jr,ic,ja,jb,jc,kp,md,me,kmax,idel,iaph, kpoint,lr,nujm
REAL(Kind=WP_Kind) rho, rhosav
REAL(Kind=WP_Kind) beta(nidm,3*mxnode/4), idi(nidm,3*mxnode/4),nid(3*mxnode/4)
REAL(Kind=WP_Kind) plmtemp(0:mxl)
REAL(Kind=WP_Kind) plmst(0:jrjtmax,nxdim,0:jrmax)
REAL(Kind=WP_Kind) djmkp(0:jrjtmax,nxdim,nydim,3)
REAL(Kind=WP_Kind) phi(nxdim,nydim,narran)
REAL(Kind=WP_Kind) betqf(nxdim,nydim,3)
REAL(Kind=WP_Kind) chi1(nxdim,nydim,3)
REAL(Kind=WP_Kind) thaph(nxdim,nydim,3)
REAL(Kind=WP_Kind) snthd(nydim,3)
REAL(Kind=WP_Kind) th(9,mxelmnt), ch(9,mxelmnt), f(mxnode)
REAL(Kind=WP_Kind) tstore(ndelve,ndelve), temp(naph,naph), bfaci
REAL(Kind=WP_Kind) umat(ndelve,naph), bigtemp,  eaph, edel, ejac
REAL(Kind=WP_Kind) upsil
REAL(Kind=WP_Kind) h1(9,npoint)
REAL(Kind=WP_Kind), ALLOCATABLE :: wx(:),x(:),wy(:)
REAL(Kind=WP_Kind), ALLOCATABLE :: thetd(:,:), xptg(:,:)
REAL(Kind=WP_Kind), ALLOCATABLE :: sqci(:), a(:)
REAL(Kind=WP_Kind) sumx,sumk,sum0,pi2,sq2i,rac,cnorm

DATA ithcll/0/,ithsub/0/
DATA little/.false./,medium/.false./,full/.false./
DATA lamsave/-1000/, rhosav/-1.0d0/,jtsav/-1/,ipsav/-1/
CALL popt('AphDel_Old  ',little,medium,full,ithcll,ithsub)
!----------------------------------------------------------------------
!  this routine transforms the wigner r-matrix in the
!  (aph) basis to the wigner r-matrix in the delves basis.
!-----------------------------------------------------------------------
!  requires subprograms lgnd11,lgndrx,cleb,legendre, and trnsfm
! also requires function djmk and its routines numgen and factor
!  ---------------------------------------------------------------------
!  rho.......hyperradius (bohr)
!  naph......number of coupled channels in the aph frame.
!  lambda....lambda quantum numbers in the aph frame.
!  ndelve......number of coupled channels in the delves frame. hopefully
!       naph is approximately equal to ndelve.
!  jtot......total angular momentum quantum number.
!  kchan.....kchan(idel) is arrangement channel (1, 2 or 3) with which
!                state idel is associated.
!  mvib......vibrational quantum numbers in the delves frame.
!  jrot3......rotational quantum numbers in the delves frame.
!  lorb3......orbital angular momentum quantum numbers in the delves
!      frame.
!  tstore....temporary storage matrix.  used here to store unitarity
!            test and also clebsch-gordan coeffs.
!  temp......temporary storage matrix used in unitarity tests.
!  umat......orthogonal u-matrix that transforms the r-matrix in the
!       aph frame to the r-matrix in the delves frame.
!  nxdim.....number of quadrature points in legendre quad on rotor angle
!  nhermt..number of quadrature points in hermite quad on delves angle.
!  nuj3....(nuj3(idel)-1)*nhermt locates start of upsil in chinuj array.
!  plmst.....array of plm(k to jrmx, nxdim, jrmx)
!  djmkp.....array of definite parity djmk(k to jrmx, nxdim,nhermt)
!  phi.......phi(nxdim,nhermt) are aph surf fcns of fixed jtot and iaph
!            at the delves quadrature points.
!  upsil.....upsil(nhermt,ndelve) are delves vibrational functions.
!----------------------------------------------------------------------
!      plmst            array for storing associated legendre
!                       polynomials.  This array is needed for the
!                       APH to delves transformation only.
!      djmkp            array for storing Wigner rotation functions.
!                       This array is needed for the APH to Delves
!                       transformation only.
!      betqf            array for storing beta angles. This array
!                       is needed for the APH to Delves transformation
!                       only.
!      chi1             array for storing chi angles.  This array is
!                       needed for the APH to Delves transformation
!                       only.
!      thaph            array for storing the APH theta angles.  This
!                       array is needed for the APH to Delves
!                       transformation only.
!      snthd            array for storing the sin of the theta Delves
!                       This array is needed for the APH to Delves
!                       transformation only.
!      nx               number of Gauss_Legendre quadrature points.
!      ny               number of Gauss_Laguerre quadrature points.
!      jrmax            value of the largest rotational quantum
!                       number used in the Delves region.
!
!  ------------------------------------------------------------------
!
!allocata storage for several matrix
!
ALLOCATE(chanl(mxbasis, 0:mxmega), nvib(mxbasis),jrot(mxbasis,0:mxmega))
IF(.Not.ALLOCATED(wx))ALLOCATE(wx(96),x(96),wy(96))
IF(.Not.ALLOCATED(thetd))ALLOCATE(thetd(maxhermt,narran),xptg(mxglegn,narran))
IF(.Not.ALLOCATED(sqci))ALLOCATE(sqci(3),a(200))

WRITE(Out_Unit,*)'oliver size of chanl,nvib,jrot=',size(chanl),size(nvib),size(jrot)
WRITE(Out_Unit,*)'oliver size of thetd,xptg=',size(thetd),size(xptg)
WRITE(Out_Unit,'(//," aph to delves projection at rho=",f10.5)')rho
IF(rho==rhosav) GOTO 2
!  ------------------------------------------------------------------
!  start calculations done on the first pass only
!  ------------------------------------------------------------------
IF(nstate<naph.or.nstate<ndelve)THEN
   WRITE(Out_Unit,*) 'nstate=', nstate, ' is < naph or ndelve=', naph, ndelve
   WRITE(Out_Unit,*) 'AphDel_Old: increase nstate.parms'
   STOP 'aphdel1'
ENDIF

WRITE(Out_Unit,*)'ndelve=',ndelve
nx=nxdim
WRITE(Out_Unit,*)'nxdim,nydim=',nxdim,nydim
rhosurf=rho
sq2i=1.d0/sqrt(2.d0)
WRITE(Out_Unit,'(/," chif1=",3f10.7)')(chif1(ic),ic=1,3)
pi2=HalfPi
!
!  CALL routines to generate numbers needed by djmk.
CALL numgen
CALL factor
!  --------------------------------------------------------------
!  construct points and weights of Gauss_Legendre quadrature.
!  this is now same in all arrangements but needs to be generalized.
!  --------------------------------------------------------------
CALL lgnd11(nx,nx,wx,x)

DO karran=1,narran
   ngor(karran)=nx
   DO ix=1,ngor(karran)
      xptg(ix,karran)=x(ix)
   ENDDO
ENDDO
!  --------------------------------------------------------------
!  get points and weights of hermite quadrature
!  temp this requires all nhermts to be equal
!  --------------------------------------------------------------
IF(nhermt(1)/=nhermt(2).or.nhermt(2)/=nhermt(3))THEN
  WRITE(Out_Unit,*)'AphDel_Old requires nhermt to be ','equal for all channels'
  WRITE(Out_Unit,*)'nhermt=',nhermt
  STOP 'AphDel_Old'
ENDIF

IF(integrat(1))THEN
  in=1
ELSEIF(integrat(2))THEN
  in=2
ELSEIF(integrat(3))THEN
  in=3
ENDIF

IF(little)WRITE(Out_Unit,*)'scheme= ',scheme
DO iy=1,nhermt(in)
   IF(scheme==1)THEN
      !  modify the weights if the exponentials don't match those of the Hermite or Laguerre quadrature exactly.
      wy(iy)=wpth(iy,in)*exp(0.5d0*xpth(iy,in)**2)
      IF(medium)WRITE(Out_Unit,'(A,3ES15.7)')'wpth,xpth,wy= ', wpth(iy,in),xpth(iy,in), wy(iy)
   ELSEIF(scheme==2)THEN
      wy(iy)=wpth(iy,in)
   ENDIF
!  calc delves theta and fcns of it at vib quad pts in each arrangement

   DO ic=1,3
      IF(.NOT.integrat(ic))GOTO 13 !tempmod
      thetd(iy,ic)=c(ic)*xpth(iy,ic)+b(ic) !!!it has been calculated in the new way
      sqci(ic)=1./sqrt(c(ic))
      IF(thetd(iy,ic)<0.0)THEN
         WRITE(Out_Unit,'(1x," for iy,ic=",2i5," thetd,c,b,y=",4e13.6)')iy,ic,thetd(iy,ic),c(ic),b(ic),xpth(iy,ic)
         thetd(iy,ic)=0.   !  for negative thetd set it barely back on the range.
      ENDIF

      IF(thetd(iy,ic)>pi2)THEN
         WRITE(Out_Unit,'(1x," for iy,ic=",2i5," thetd,c,b,y=",4e13.6)')iy,ic,thetd(iy,ic),c(ic),b(ic),xpth(iy,ic)
         STOP 'AphDel_Old'
      ENDIF
      !  --------------------------------------------------------------
      !  calc theta, chi, and betaqf at each of the quadrature points.
      !  --------------------------------------------------------------
      CALL aphdef(nx,iy,ic,nxdim,nydim,x,snthd,thaph,chi1,chif1,betqf,thetd(iy,ic))
13 ENDDO
ENDDO

2     CONTINUE
WRITE(Out_Unit,*)'nx,jrjtmax,jrmax=',nx,jrjtmax,jrmax
!------------------------------------------------------------------
!  start calculation done at each jtot
!  ----------------------------------------------------------------
!  test parity and construct quantities needed repeatedly
!  parity check for jtot=0 case
IF(jtot==0.AND.parity/=0) GOTO 999
!  construct matrix of associated legendre polynomials,
!  normalized to unity on x=(-1,1).
WRITE(*,*) "Calling Legendre from AphDel_Old"
DO  ix=1,nx
   DO  k=0,jrjtmax
      CALL legendre (jrmax,k,plmtemp,x(ix))
      DO  jr=k,jrmax
         plmst(k,ix,jr)=plmtemp(jr)
      ENDDO
   ENDDO
ENDDO
!  construct array of clebsch-gordan coeffs
DO idel=1,ndelve
   jr=jrot3(idel)
   IF(Medium)WRITE(Out_Unit,*)'idel,jrot3(idel),lorb3(idel)', idel,jrot3(idel),lorb3(idel)
   ja=2*jr
   jb=2*jtot
   jc=2*lorb3(idel)
   kmax=min0(jr,jtot)
   DO k=0,kmax
      kp=k+1
      md=2*k
      me=-md
      CALL clebp(ja,jb,jc,md,me,0,1,rac)
      !  use phase to avoid calculating a square root.
      !  store clebsch-gordan coeffs in tstore.
      tstore(kp,idel)=rac*(-1)**(jr-k)
   ENDDO
ENDDO
CALL flusher(Out_Unit)

WRITE(Msg_Unit,*) 'tstore done. start loop to naph=', naph
WRITE(Out_Unit,*) 'tstore done. start loop to naph=', naph
!  -----------------------------------------------------------------
!  begin construction of umat with loop on aph quantum number index
!  -----------------------------------------------------------------
DO iaph=1,naph
   !  -------------------------------------------------------------
   !  get the aph surface functions phi by interpolation.
   !  method used depends on type of surface functions.
   !  -------------------------------------------------------------
   lam=lambda(iaph)
   tea=tee(iaph)
   IF(sfuntype=='FEM')THEN
      !  CALL interp(nxdim,nydim,nx,nhermt(in),chi1,thaph,phi,id,iaph,3, th,ch,f,beta,idi,nid,nidm,rho,h1,jx,jy,jc1,kpoint, xnode,idni
      WRITE(Msg_Unit,*)'FEM for surface functions is not wanted!'
      STOP 'AphDel_Old'
   ELSEIF(sfuntype=='DVR')THEN
      IF(iaph==1)CALL readinfbr  ! Change by X.L. 1-2-2006
      CALL intfbr(nxdim,nydim,nx,nhermt(in),chi1,thaph,phi,iaph,3,rho)
   ELSEIF(sfuntype=='ABM')THEN
      CALL intbas(nxdim,nydim,ngor,nhermt,xptg,thetd,phi,tea,lam, rho)
   ELSE
      WRITE(Out_Unit,*)' sfuntype not specified as one ',' of allowed types'
      WRITE(Out_Unit,'(1x,a3)')sfuntype
      STOP 'AphDel_Old'
   ENDIF

   !  -------------------------------------------------------------
   !  construct wigner little d fcns of good parity, djmkp(k,ix,iy,ic).
   !  k is delves index.  The other indexes specify betqf angle between
   !  frames.  j=J and m=lam (aph index) are assumed fixed.
   !  If memory permits, adding lam index allows moving this outside loop.
   !  reconstruct only when necessary.
   !  -------------------------------------------------------------
   IF(lam==lamsave) GOTO 31
   cnorm=sq2i
   IF(lam==0) cnorm=0.5d0
   DO  ic=1,3
      IF(integrat(ic))THEN
         CALL wigner(nx,nhermt(ic),nxdim,nydim,jrjtmax,betqf,djmkp, ic,jtot,lam,a,cnorm,parity)
      ENDIF
   ENDDO
   lamsave=lam
31 CONTINUE
!  -------------------------------------------------------------
!  start main loop on delves quantum number index
!  -------------------------------------------------------------
   DO idel=1,ndelve
      umat(idel,iaph)=0.d0
      ic=kchan(idel)
      jr=jrot3(idel)
      lr=lorb3(idel)
      kmax=min0(jr,jtot)
      nujm=(nuj3(idel)-1)*nhermt(ic)
      IF(full)THEN
         WRITE(Out_Unit,*)'idel=',idel,'   nuj3=',nuj3(idel), '   nhermt=',nhermt(ic),'   nujm=',nujm
      ENDIF
      IF((-1)**(jr+lr+parity)/=1) GOTO 998
      !  the above is the final parity check
      !      IF(medium)THEN
      !         WRITE(Out_Unit,*)'for ic=',ic,' iaph=',iaph,' idel=',idel
      !         WRITE(Out_Unit,*)'nhermt=',nhermt(ic),' jr=',jr,' lr=',lr,' nujm=',nujm
      !         WRITE(Out_Unit,*)'nxdim=',nxdim,' parity=',parity,' kmax=',kmax
      !      ENDIF
      IF(full)THEN
         WRITE(Out_Unit,*)'iy upsil c wy  snthd sumx sumx*snthd umat(idel,iaph)'
      ENDIF
      !  ---------------------------------------------------------------
      !  start integration loops
      !  ----------------------------------------------------------------
      DO iy=1,nhermt(ic)
         sumx=0.
         DO  ix=1,nx
            sumk=0.
            sum0=tstore(1,idel)*plmst(0,ix,jr)*djmkp(0,ix,iy,ic)
            !------------------------------------------------------------------
            !  note.  this version goes directly from BF APH to SF Delves.
            !  IF indexes are rearranged, it could GOTO BF Delves, omitting this
            !  inner loop on k and vectorizing better.  THEN, the resulting umat
            !  can be converted to the present one by similarity transformation
            !  using a matrix of C-G coeffs similar to tstore.
            !------------------------------------------------------------------
            DO  k=1,kmax
               kp=k+1
               sumk=sumk+tstore(kp,idel)*plmst(k,ix,jr)*djmkp(k,ix,iy,ic)
            ENDDO
            sumk= sumk+sumk+sum0
            sumx=sumx+phi(ix,iy,ic)*wx(ix)*sumk
            !------------------------------------------------------------------
            !     WRITE(Out_Unit,*)'sumk,sum0,phi,wx,sumx='
            !     WRITE(Out_Unit,*)sumk,sum0,phi(ix,iy,ic),wx(ix),sumx
            !------------------------------------------------------------------
         ENDDO
         !  jacobian is presently for integration over delves theta.
         !  in=iy+nujm
         upsil=chinuj(iy+nujm)*sqci(ic)
         umat(idel,iaph)=umat(idel,iaph)+upsil*c(ic)*wy(iy)*snthd(iy,ic)*sumx*2
         IF(full)THEN
            IF(iy==1)WRITE(Out_Unit,'(1x,6(A,I5))')'idel=',idel,'   channel=',ic,'   mvib=',mvib3(idel),'   jrot=',jr,'   lorb=',lr,'   iaph=',iaph
            IF(ABS(umat(idel,iaph))>0.001) WRITE(Out_Unit,951)iy,upsil,c(ic),wy(iy), snthd(iy,ic),sumx,2*sumx*snthd(iy,ic),umat(idel,iaph)
         ENDIF
 951     FORMAT(1x,i4,7ES14.6)
      ENDDO
   ENDDO
ENDDO
!  -------------------------------------------------------------
!  END loops.  completes construction of umat.
!  -------------------------------------------------------------

IF(little)THEN
   WRITE(Out_Unit,'(//," umat")')
   CALL MxOutL(umat,ndelve,naph,0,'space','aph  ')
ENDIF
!----------------------------------------------------------------
!used to compare the ABM code Energy to Delves energy
! by X.L. 12-08-06
!---------------------------------------------------------------
IF(.True.)THEN !TmpModGregParker
   OPEN(Unit=Eaph_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/eaph.txt')
   OPEN(Unit=Compare_APH_Jac_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/Compare_APH_Jac.txt')
   WRITE(Compare_APH_Jac_Unit,*)'rho=',rho
   WRITE(Compare_APH_Jac_Unit,'(3A10,7A5,3A15,A8)')'iaph','idel','MxOvrlp','ch','elec','nvib','jrot','lorb','mega','nuj3','E_Delves(eV)','E_APH(eV)','E_Jacobi(eV)','%error'
      WRITE(Msg_Unit,*)'iaph=',iaph
      READ(Eaph_Unit,*)eaph
      bigtemp=0.0d0
      DO idel=1,ndelve
         IF(ABS(umat(idel,iaph))>bigtemp)THEN
            bigtemp=ABS(umat(idel,iaph))
            indexdel=idel
         ENDIF
      ENDDO
      WRITE(Out_Unit,*)'iaph=',iaph,'idel=',indexdel,bigtemp
      edel=(Etot-xksq3(indexdel)/usys2)*autoev
      ejac=Energy3_Jacobi(indexdel)*autoev
      eaph=Energy3_APH(iaph)*autoev
      WRITE(Compare_APH_Jac_Unit,'(A5,i5,A5,i5,F10.7,6i5,i5,3es15.7,F7.2,A1)')' iaph=',iaph,' idel=',indexdel, bigtemp, &
            kchan(indexdel), elect3(indexdel), mvib3(indexdel), jrot3(indexdel), &
            lorb3(indexdel), mega3(indexdel),  nuj3(indexdel),  edel, eaph, ejac, abs((eaph-edel)/edel)*100d0,'%'
   CLOSE(Eaph_Unit)
   CLOSE(Compare_APH_Jac_Unit)
ENDIF

!  ------------------------------------------------------------------
!  test orthogonality of umat
!  ------------------------------------------------------------------
!  the first direction tests the adequacy of the delves basis
CALL delchk(naph,temp,ndelve,umat)
IF(little)THEN
   WRITE(Out_Unit,'(//,"utranspose*u"/)')
   CALL MxOutL(temp,naph,naph,0,'aph  ','aph  ')
ENDIF

AphDel_All=.False.  ! Controls the amount of output in APH_To_Delves.txt
DO  iaph=1,naph
   IF(ABS(temp(iaph,iaph))>1.1)THEN
      WRITE(Out_Unit,7768)iaph,temp(iaph,iaph)
   ELSEIF(ABS(temp(iaph,iaph))<0.9)THEN
      WRITE(Out_Unit,7768)iaph,temp(iaph,iaph)
   ELSE
      IF(AphDel_All)WRITE(Out_Unit,7765)iaph,temp(iaph,iaph)
   ENDIF
ENDDO

 7768 FORMAT(/,'***warning**** for iaph=',i5,' the delves basis is not adequate',5x,f10.5)
 7765 FORMAT(/,10x,'Fantastic for iaph=',i5, 'the delves basis is adequate',5x,f10.5)
!  the following direction tests the aph basis.
CALL aphchk(ndelve,tstore,naph,umat)
IF(little)THEN
   WRITE(Out_Unit,'(//,"u*utranspose",/)')
   CALL MxOutL(tstore,ndelve,ndelve,0,'space','space')
ENDIF

DO idel=1,ndelve
   IF(ABS(tstore(idel,idel))>1.1)THEN
       WRITE(Out_Unit,8768)idel,kchan(idel), mvib3(idel),jrot3(idel), lorb3(idel),tstore(idel,idel)
   ELSEIF(ABS(tstore(idel,idel))<0.4)THEN
       WRITE(Out_Unit,8768)idel,kchan(idel), mvib3(idel),jrot3(idel), lorb3(idel),tstore(idel,idel)
   ELSE
       IF(AphDel_All)WRITE(Out_Unit,8769)idel,kchan(idel), mvib3(idel),jrot3(idel), lorb3(idel),tstore(idel,idel)
   ENDIF
ENDDO

 8768 FORMAT(/,'warning APH basis incomplete: idel=',i4,' channel=',i4,' mvib=',i4,' jrot=',i4,' lorb=',i4,f8.5)
 8769 FORMAT(/,'Fantastic APH basis complete: idel=',i4,' channel=',i4,' mvib=',i4,' jrot=',i4,' lorb=',i4,f8.5)

DO  i = 1, ndelve
   elect3(i) = 0
ENDDO

WRITE(Out_Unit,*)'naph,ndelve=', naph, ndelve

IF(little)THEN
   WRITE(Out_Unit,*)'Etot=', Etot
ENDIF

!  save umat on disc umat_unit
WRITE(Out_Unit,*)'Quantum Numbers used in APH to Delves Transformation'
DO I=1,ndelve
   energy3_delves(i)=Etot-xksq3(i)/usys2
ENDDO
CALL Quant_Out (ndelve, kchan, elect3, mvib3, jrot3, lorb3, energy3_delves, xksq3,'AphDel_Old')
OPEN(Unit=umat_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Umat.bin',form='unformatted', status='unknown')
WRITE(umat_unit)naph, ndelve
WRITE(umat_unit)(Etot-xksq3(i)/usys2,i=1,ndelve)
WRITE(umat_unit)kchan
WRITE(umat_unit)elect3
WRITE(umat_unit)mvib3
WRITE(umat_unit)jrot3
WRITE(umat_unit)lorb3
WRITE(umat_unit)umat
CLOSE(Unit=umat_unit,status='keep')

jtsav=jtot
ipsav=parity
rhosav=rho


!
!deallocata storage for several matrix
!
DEALLOCATE(chanl, nvib,jrot)
!DEALLOCATE(wx,x,wy)
!DEALLOCATE(thetd,xptg)
!DEALLOCATE(sqci,a)
!

RETURN
!  ------------------------------------------------------------------
!  access the following only for parity error
!  ------------------------------------------------------------------
998   WRITE(Out_Unit,62)jtot,jr,lr,parity
62    FORMAT(/,'parity error.  for jtot=',i3,' it tried jrot3=',i3, ' and lorb3=',i3,' with parity=',i3)
STOP 'AphDel_Old'
999   WRITE(Out_Unit,61)jtot,parity
61    FORMAT(/,'parity error.  for jtot=',i3,' it tried parity=',i3)
STOP 'AphDel_Old'
ENDSUBROUTINE AphDel_Old
