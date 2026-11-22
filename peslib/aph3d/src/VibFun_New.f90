SUBROUTINE VibFun_New (noscil,jang,nhermt,haml,vect,chinuj,nulow,nuhigh,nudim,eneg,alpha,rx,z,w,re,Evib,h,scheme)
!-----------------------------------------------------------------------
! this routine determines the vibrational eigenfunctions and
! eigenvalues of a diatomic molecule in a fixed rotational
! state using a harmonic oscillator basis set to diagonalize the hamiltonian.
! on entering
! noscil is the number of harmonic oscilators to be used
! jang is the orbital angular momentum quantum of the diatomic molecule.
! nhermt is the number of Gauss_Hermite integration points.
! nulow and nuhigh label the minimun and maximun vibrational levels
! rx is the ratio of the harmonic oscillator position
! to the equilibrium position of the diatomic molecule.
! alpha is
!-----------------------------------------------------------------------
! on RETURN
! vect contains the eigenvectors.
! chinuj contains the unweighted eigenfunctions.
! z Gauss_Hermite integration points.
! w Gauss_Hermite wieghts.
!-----------------------------------------------------------------------
! this program was written by g. a. parker.
! this routine is called by: qlevel
! this routine calls: popt
!-----------------------------------------------------------------------
USE Numbers_Module
USE Narran_Module
USE FileUnits_Module 
USE Arrch_Module
USE Spectro_Module
USE Oops_Module
USE Das_Module
USE Convrsns_Module
USE PotZero_Module
USE Masses_Module
USE RhoNum_Module
USE DiatomicPot_Module, ONLY: DiatomicPot
IMPLICIT NONE
INTEGER scheme, noscil, jang, nulow, nuhigh, nhermt, nudim
INTEGER IthCll, IthSub, n, np1, i, mprnt, nprnt, np, j, jnew, k
INTEGER im1, inew, ms, Info, lwork
LOGICAL little,medium,full
!-----------------------------------------------------------------------
REAL(Kind=WP_Kind), ALLOCATABLE :: Work(:)
REAL(Kind=WP_Kind) vect(noscil,noscil), chinuj(nhermt,nudim),vect2(noscil,noscil)
REAL(Kind=WP_Kind) eneg(noscil), Evib(noscil), h(noscil),eneg2(noscil),Atheta
REAL(Kind=WP_Kind) haml(noscil,noscil), w(nhermt), z(nhermt), hdn(maxhermt), zlit(nhermt), thetd(nhermt,Narran), hd2, hd2n(maxhermt)
REAL(Kind=WP_Kind) Alpha, Re,  rx 
REAL(Kind=WP_Kind) Sqrj, hs, vph, alphas, alphai, hsinv, r, x, pot, rre
REAL(Kind=WP_Kind) cterm, wi, effpot, term, perr, enegwn, trans, enegev, err, enegi
REAL(Kind=WP_Kind) sumtemp, rotm, potbc, hd,az(3),bz(3),cz(3),factor1,factor2,rhonow,slit(nhermt),a(3),zmc,ax(3),bfaci(nhermt,Narran)
INTEGER Ichange
DATA ithcll/0/,ithsub/0/,Ichange/2/
DATA little/.false./,medium/.true./,full/.false./
CALL popt('vibfun_new  ', little, medium, full, ithcll, ithsub)
!little = .true.
!medium = .true.
!full = .true.
lwork=3*Noscil-1
ALLOCATE(work(lwork))
IF(little)WRITE(Out_Unit,*)"VibFun_New Called"


!---------------------------------------------------------------------
!Ichange=0 is for the old way using gaussx
!       =1 is the new way for Jacobi region
!       =2 is the new way for Delves part
!---------------------------------------------------------------------
IF(.not.ioops)THEN
   Ichange=1
ELSE
   Ichange=2
ENDIF

IF(nulow/=0)THEN
   WRITE(Out_Unit,*)'nulow=',nulow
   WRITE(Msg_Unit,*)'VibFun_New: nulow/=0  nulow=',nulow,' karran=',karran
   STOP 'VibFun_New nulow'
ENDIF
IF (nrigid)THEN
   IF((nhermt/=1).or.(nulow/=0).or.(nuhigh/=0).or.(nudim/=1))THEN
      WRITE(Out_Unit,310) nhermt,nulow,nuhigh,nudim
      WRITE(Msg_Unit,*)'VibFun_New: nrigid constraint ',nhermt,nulow,nuhigh,nudim
      STOP 'VibFun_New nrigid'
   ENDIF
   alpha=1.d0
   re=1.d0/sqrt(usys2*be(karran) )
   z(1)=0.d0
   CALL dhep (h,z(1),0)
   w(1)=1.d0/(h(1)*h(1))
ENDIF
!print*,'oliver vib_2'
!-----------------------------------------------------------------------
! determine zeros and weights
!-----------------------------------------------------------------------
IF(.NOT.nrigid)THEN
   IF(nhermt<=noscil)THEN
   WRITE(Out_Unit,*)'nhermt must be larger than noscil',nhermt,noscil
   WRITE(Msg_Unit,*)'VibFun_New: nhermt<=noscil ',nhermt,noscil,' karran=',karran
   STOP 'VibFun_New nhermt'
ENDIF

IF(scheme==1)THEN
   CALL hermit (nhermt,z,w,0)
   IF(Ichange==1)THEN   !!calculate az,bz,cz in Jacobi jacobi region
      call basis_jacobi(z(1), az(karran), bz(karran), cz(karran))   
   ELSEIF(Ichange==2)THEN   !!calculate theta_del and az,bz,cz in delves region
      call basis_delves(z,thetd,az,bz,cz,karran)
   ELSE
      Write(Msg_Unit,*)'Ichange/= 1 or 2 : ichange=',ichange
      STOP 'VibFun_New2'
   ENDIF
      
ELSEIF(scheme==2)THEN
   CALL laguer(z,w,nhermt,1)
ENDIF
ENDIF
!-----------------------------------------------------------------------
! if ioops=.false. it is calculated in Jacobi basis
! we have to change az,bz,cz
!-----------------------------------------------------------------------
IF(Ichange==2)rhonow=rhonum

IF((.not.ioops).AND.Ichange==1)THEN
   a(1)=az(karran)
   a(3)=-bz(karran)
   DO i= 1, nhermt
      zmc = z(i) - cz(karran)
      a(2)=-zmc
      zmc=-(a(2)+sign(1.0d0,a(2))*sqrt(a(2)**2-4.0d0*a(1)*a(3)))
      zmc=zmc/2.0d0
      ax(1)=zmc/a(1)
      ax(2)=a(3)/zmc
      slit(i)=max(ax(1),ax(2))
   ENDDO
ENDIF

sqrj=jang*(jang+1)
IF(Ichange==0)hs=1.0d0/(usys2*re*re)
IF(Ichange==1.or.Ichange==2)hs=1.0d0/usys2
!-----------------------------------------------------------------------
! calculate experimental energies for a comparison.
!-----------------------------------------------------------------------
IF(full)WRITE(Out_Unit,*)'DO 30 , nulow, nuhigh',nulow,nuhigh
DO np1=nulow+1,nuhigh+1
   IF(scheme==1)THEN         ! non-coulombic potential
      vph=REAL(np1,WP_Kind)-0.5d0
      Evib(np1)=we(karran)*vph-wexe(karran)*vph**2+weye(karran)*vph**3+be(karran)*sqrj-de(karran)*sqrj*sqrj-alfe(karran)*vph*sqrj+wzero(karran)
   ELSEIF(scheme==2)THEN    ! compute the hydrogenic energies
      n=np1+jang
      Evib(np1) = -(usys/dscale(karran)**2)/(2.d0*n**2)+vzero
   ENDIF
ENDDO

hd2n=0.0d0
hdn=0.0d0
IF(alpha==0.d0)THEN
   write(*,*)"alpha=",alpha
   stop "alpha==0 in VibFun_New"
endif
alphas=alpha*alpha
alphai=1.0d0/alpha
hsinv=1.0d0/hs
IF (medium)THEN       ! IF desired print the input parameters.
   WRITE(Out_Unit,250) alpha,rx
   IF(.NOT.nrigid)THEN
      WRITE(Out_Unit,260) noscil
      WRITE(Out_Unit,270) nhermt
      mprnt=nulow
      nprnt=nuhigh
      WRITE(Out_Unit,280) mprnt,nprnt
   ENDIF
   WRITE(Out_Unit,170) hs,hsinv
   WRITE(Out_Unit,240) jang
ENDIF

IF(.NOT.nrigid)THEN
   haml=0.d0
!-----------------------------------------------------------------------
! preform the integrations
!-----------------------------------------------------------------------
DO i=1,nhermt
IF(Ichange==0)THEN
   CALL WavDer_New(h,hd,hdn,hd2,hd2n,z(i),jang,noscil,scheme,alphas)
   hd2=0.0d0
   hd2n=0.0d0
   r=alphai*z(i)+rx
   x=r-1.0d0
   rre=r*re
ELSEIF(Ichange==1)THEN
    CALL WavDer_New(h,hd,hdn,hd2,hd2n,z(i),jang,noscil,scheme,One)
    rre=slit(i)
    r=slit(i)
ELSEIF(Ichange==2)THEN
   CALL WavDer_New(h,hd,hdn,hd2,hd2n,z(i),jang,noscil,scheme,One)
   rre=rhonow*sin(thetd(i,karran))
   r=rre!*cos(thetd(i,karran))
ELSE
   STOP 'wrong in vibfun_new'
ENDIF

   pot=potbc(rre)
IF(Ichange==1)THEN
   hd2=0.0d0
   hdn=0.0d0
   Atheta=-bz(karran)/(az(karran)*slit(i)**3+bz(karran)*slit(i))
   factor1=az(karran)+bz(karran)/slit(i)**2
   hd=Atheta-factor1*z(i)
   hd2n=-hd2n*factor1
ELSEIF(Ichange==3)THEN
   hd2=0.0d0
   hdn=0.0d0
   Atheta=0.0d0
   factor1=az(karran)+bz(karran)/slit(i)**2
   hd=Atheta-factor1*z(i)
   hd2n=-hd2n*factor1
ENDIF

IF(ICHANGE==2)THEN
   Atheta=-2.0d0/tan(2.0d0*thetd(i,karran))+sin(2.0d0*thetd(i,karran))*(az(karran)-bz(karran))/2.0d0/(az(karran)*sin(thetd(i,karran))**2+bz(karran)*cos(thetd(i,karran))**2)
   factor1=az(karran)/(cos(thetd(i,karran)))**2+bz(karran)/(sin(thetd(i,karran)))**2
   factor2=0.0d0
   hd=Atheta-factor1*z(i)
   hd2n=(-1.0d0)*factor1*hd2n
   hd2=0.0d0
   bfaci(i,karran)=sin(2.0d0*thetd(i,karran))/sqrt(az(karran)*(sin(thetd(i,karran)))**2+bz(karran)*(cos(thetd(i,karran)))**2)/2.0d0
ENDIF

      IF (ioops)THEN
IF(Ichange==0)THEN
   cterm=hd+hd2+sqrj*re**2/(rhoj*sin(rre/rhoj))**2+pot*hsinv
ELSE IF(ICHANGE==1)THEN
   cterm=sqrj/(rhoj*sin(rre/rhoj))**2+pot*hsinv
ELSE IF(ICHANGE==2)THEN
   cterm=sqrj/(rhoj*sin(thetd(i,karran))*cos(thetd(i,karran)))**2+pot*hsinv-2.0d0/4.0d0/(rhonow)**2
ELSE
   STOP 'wrong in vibfun_new'
ENDIF
      ELSE
IF(Ichange==0)THEN
         cterm=hd+hd2+sqrj/r**2+pot*hsinv
ElSEIF(Ichange==1)THEN
         cterm=sqrj/r**2+pot*hsinv
ELSEIF(Ichange==2)THEN
         cterm=sqrj/r**2+pot*hsinv-1.0d0/4.0d0/(rhonow)**2
ENDIF
      ENDIF
      wi=w(i)
      effpot=sqrj/r**2/hsinv+pot
      IF(medium) WRITE(Out_Unit,160) i,rre,pot,effpot
      DO n=1,noscil
IF(Ichange==0)THEN
        term=wi*(hdn(n)+cterm)*h(n)
ELSEIF(Ichange==1)THEN
         term=wi*cterm*h(n)
ELSEIF(IChange==2)THEN
        term=wi*cterm*h(n)
ENDIF
         DO np=1,noscil
            haml(n,np)=haml(n,np)+h(np)*term
	    IF(Ichange==2)THEN
	       IF(n>1.AND.np>1)THEN
	        haml(n,np)=haml(n,np)+(wi*(hd*h(n)+hd2n(n)*h(n-1))*(hd*h(np)+hd2n(np)*h(np-1)))/(rhonow)**2
	       ELSE IF(n>1.AND.np==1)THEN
                haml(n,np)=haml(n,np)+(wi*(hd*h(n)+hd2n(n)*h(n-1))*(hd*h(np)))/(rhonow)**2
	       ELSE IF(n==1.AND.np>1)THEN
                haml(n,np)=haml(n,np)+(wi*(hd*h(n))*(hd*h(np)+hd2n(np)*h(np-1)))/(rhonow)**2
		ELSE IF(n==1.AND.np==1)THEN
                haml(n,np)=haml(n,np)+(wi*(hd*h(n))*(hd*h(np)))/(rhonow)**2
	       ENDIF
	    ELSEIF(Ichange==1)THEN
               IF(n>1.AND.np>1)THEN
	        haml(n,np)=haml(n,np)+(wi*(hd*h(n)+hd2n(n)*h(n-1))*(hd*h(np)+hd2n(np)*h(np-1)))
	       ELSE IF(n>1.AND.np==1)THEN
                haml(n,np)=haml(n,np)+(wi*(hd*h(n)+hd2n(n)*h(n-1))*(hd*h(np)))
	       ELSE IF(n==1.AND.np>1)THEN
                haml(n,np)=haml(n,np)+(wi*(hd*h(n))*(hd*h(np)+hd2n(np)*h(np-1)))
		ELSE IF(n==1.AND.np==1)THEN
                haml(n,np)=haml(n,np)+(wi*(hd*h(n))*(hd*h(np)))
	    ENDIF
	   ENDIF
      ENDDO
   ENDDO
ENDDO

 IF(medium)THEN
    WRITE(Out_Unit,*)'haml(ev)=',noscil
    DO n=1,Min(noscil,5)
       WRITE(Out_Unit,'(1x,5es15.7)')(haml(n,np)*hs*autoev, np=1,Min(noscil,5))
    ENDDO
 ENDIF

!-----------------------------------------------------------------------
! diagonalize the hamiltonian to determine the eigenvectors.
!-----------------------------------------------------------------------
CALL DSYEV('V','U',noscil,haml,noscil,eneg,work,lwork,info)
vect=haml
IF(Info/=0)THEN
   WRITE(Out_Unit,*)'Error in diagonalization routine: Info=',Info
   WRITE(Msg_Unit,*)'VibFun_New: diag error Info=',Info,' noscil=',noscil,' karran=',karran, &
     ' wzero=',wzero(karran),' nulow=',nulow,' nuhigh=',nuhigh
   STOP 'VibFun_New'
ENDIF
!-----------------------------------------------------------------------
! determine the unweighted eigenfunctions at the gaussian integration points.
!-----------------------------------------------------------------------
DO i=1,nhermt
   IF(Ichange==1.or.Ichange==2)THEN
      CALL WavDer_New(h,hd,hdn,hd2,hd2n,z(i),jang,noscil,scheme,One)
      IF(Ichange==2)h=h/bfaci(i,karran)
   ELSEIF(Ichange==0)THEN
      CALL WavDer_New(h,hd,hdn,hd2,hd2n,z(i),jang,noscil,scheme,alphas)
   ENDIF

   DO j=nulow+1,nuhigh+1
      jnew=j-nulow
      sumtemp=0.0d0
      DO k=1,noscil
         sumtemp=sumtemp+vect(k,j)*h(k)
      ENDDO
      chinuj(i,jnew)=sumtemp
   ENDDO
ENDDO

ELSE
   STOP 'VibFun_New wrong'
   rotm=usys*re**2
   eneg(1)=(.5d0*we(karran)+jang*(jang+1)/(2.d0*rotm))*hsinv
   chinuj(1,1)=h(1)
ENDIF
!-----------------------------------------------------------------------
! IF desired compare with the experimental energies.
!-----------------------------------------------------------------------
IF(little) WRITE(Out_Unit,300)
DO i=nulow+1,nuhigh+1
   im1=i-1
   enegi=eneg(i)*hs
   err=enegi-Evib(i)
   perr=err
   IF(Evib(i)/=0.d0) perr=100.0d0*err/MAX(ABS(Evib(i)),ABS(enegi))
   enegwn=enegi/cmm1toau
   inew=i-nulow
   eneg(inew)=enegi
   trans=0.0
   enegev=enegwn*cmm1toau*autoev
   IF(i>nulow+1)THEN
      trans=(eneg(inew)-eneg(inew-1))/cmm1toau
   ENDIF
   IF(ABS(perr)>10.0d0.and.i<nuhigh-1)THEN
     WRITE(Out_Unit,*)
     WRITE(Out_Unit,*)"karran=",karran
     WRITE(Out_Unit,*)"DiatomicPot(VibFun_New)= ",TRIM(DiatomicPot)
     WRITE(Out_Unit,*)'Error: VibFun_New has a problem, perr=',perr
     WRITE(Out_Unit,300)
     WRITE(Out_Unit,180) im1, jang, enegwn, enegi, Evib(i), err, perr, trans, enegev
        
     WRITE(Msg_Unit,*)
     WRITE(Msg_Unit,*)"karran=",karran
     WRITE(Msg_Unit,*)"DiatomicPot(VibFun_New)= ",TRIM(DiatomicPot)
     WRITE(Msg_Unit,*)'Error: VibFun_New has a problem, perr=',perr
     WRITE(Msg_Unit,300)
     WRITE(Msg_Unit,180) im1, jang, enegwn, enegi, Evib(i), err, perr, trans, enegev
     
     WRITE(Msg_Unit,*)"Stopping in VibFun_New"
     !STOP 'Error in VibFun_New'
   ELSE
      IF(little)WRITE(Out_Unit,180)im1, jang, enegwn, enegi, Evib(i), err, perr, trans, enegev
   ENDIF
ENDDO

IF (.NOT.nrigid)THEN
   IF(medium)THEN   ! IF desired print eigenvectors and eigenfunctions.
      ms=0
      WRITE(Out_Unit,230)
      CALL MxOutD(vect,noscil,noscil,ms)
      WRITE(Out_Unit,190)
      DO k=1,nhermt
         x=alphai*z(k)+rx-1.0d0
         r=(1.d0+x)*re
         WRITE(Out_Unit,200) z(k),x,r,w(k)
      ENDDO
      WRITE(Out_Unit,220)
      CALL MxOutD(chinuj,nhermt,nudim,ms)
   ENDIF
ENDIF
DEALLOCATE(work)
RETURN
  160 FORMAT(1x,'i=',i5,'  r=',1pe14.7,'  pot=',1pe14.7,'  effpot=',1pe14.7)
  170 FORMAT(/,"hs=",e14.7,5x,"hsinv=",e14.7)
  180 FORMAT(1x,2i3,f15.5,4f15.10,f13.5,f15.10)
  190 FORMAT(/,6x,"z",18x,"x",18x,"r",16x,"weight")
  200 FORMAT(1x,e14.7,5x,e14.7,5x,e14.7,5x,e14.7)
  210 FORMAT(/,"weighted wavefunctions")
  220 FORMAT(/,"non-weighted wavefunctions")
  230 FORMAT(/,"vibrational part of vib-rot eigenvectors")
  240 FORMAT(//,"j=",i5)
  250 FORMAT(/,"alpha=",f10.5,5x,"rx=",f10.5)
  260 FORMAT(/,"number of oscillator basis functions =",i5)
  270 FORMAT(1x,"number of integration points =",i5)
  280 FORMAT(1x,"vibrational states are nu=",i5,3x,"to nu=",i5)
  290 FORMAT(///,"the vibrational energies")
  300 FORMAT(1x,"  nu  j  e-calc(w.n.)",3x,"e-calc(a.u.)",4x,"e-exp(a.u.)",  &
               10x,"error",2x,"percent error",8x,"trans",5x,"e-calc(ev)")
  310 FORMAT(1x,"parameters not valid for rigid rotor option,"//,'nhermt',i5,'  nulow',i5,' nuhigh',i5,'  nudim',i5)
ENDSUBROUTINE VibFun_New
