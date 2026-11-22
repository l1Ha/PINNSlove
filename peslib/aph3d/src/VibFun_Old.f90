SUBROUTINE VibFun_Old (noscil,jang,nhermt,haml,vect,chinuj,nulow,nuhigh,nudim,eneg,alpha,rx,z,w,re,Evib,h,scheme)
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
USE Narran_Module
USE FileUnits_Module 
USE Arrch_Module
USE Spectro_Module
USE Oops_Module
USE Das_Module
USE Convrsns_Module
USE PotZero_Module
USE Masses_Module
IMPLICIT NONE
LOGICAL little,medium,full
INTEGER scheme, noscil, nhermt, nudim, ithcll, ithsub, nulow, nuhigh, inew
INTEGER jang, np1, n, mprnt, nprnt, npos, i, nnp, np, j, jnew, k, im1, index, ms
REAL(Kind=WP_Kind) vect(noscil,noscil), chinuj(nhermt,nudim)
REAL(Kind=WP_Kind) eneg(noscil), Evib(noscil), h(noscil+1)
REAL(Kind=WP_Kind) haml(nhermt*nudim), w(nhermt), z(nhermt), hdn(maxhermt)
REAL(Kind=WP_Kind) alpha, re, sqrj, hs, vph, alphas, alphai, hsinv, rx
REAL(Kind=WP_Kind) rre, pot, potbc, cterm, wi, effpot, term, r, x, sum, rotm
REAL(Kind=WP_Kind) err, perr, enegwn, trans, enegev, hd, enegi
DATA ithcll/0/,ithsub/0/
DATA little/.false./,medium/.false./,full/.false./
CALL popt('vibfun_old  ',little,medium,full,ithcll,ithsub)
!little = .true.
!medium = .true.
!full = .true.
IF(little)WRITE(Out_Unit,*)"VibFun_Old Called"

IF(nulow/=0)THEN
   WRITE(Out_Unit,*)'nulow=',nulow
   STOP 'vibfun_Old'
ENDIF

IF(nrigid)THEN
   IF(nhermt/=1.or.nulow/=0.or.nuhigh/=0.or.nudim/=1)THEN
      WRITE(Out_Unit,310) nhermt,nulow,nuhigh,nudim
      STOP 'vibfun_Old'
   ENDIF
   alpha=1.d0
   re=1.d0/sqrt(usys2*be(karran) )
   z(1)=0.d0
   CALL dhep (h,z(1),0)
   w(1)=1.d0/(h(1)*h(1))
ELSE
   IF(nhermt<=noscil)THEN
      WRITE(Out_Unit,*)'nhermt must be larger than noscil',nhermt,noscil
      STOP 'vibfun_Old'
   ENDIF
   IF(scheme==1)THEN
      CALL hermit (nhermt,z,w,0) ! determine Hermite zeros and weights
   ELSEIF(scheme==2)THEN
      CALL laguer(z,w,nhermt,1)  ! determine Laguerre zeros and weights
   ENDIF
ENDIF


sqrj=jang*(jang+1)
hs=1.0d0/(usys2*re*re)
IF(full)THEN
   WRITE(Out_Unit,*)'make sqrj, jang',jang
   WRITE(Out_Unit,*)'make hs, usys2, re',usys2,re
   WRITE(Out_Unit,*)'DO np1=nulow, nuhigh',nulow,nuhigh
ENDIF
!-----------------------------------------------------------------------
! calculate experimental energies for a comparison.
!-----------------------------------------------------------------------
DO np1=nulow+1,nuhigh+1
   IF(scheme==1)THEN         ! non-coulombic potential
      vph=REAL(np1,WP_Kind)-0.5d0
      Evib(np1)=we(karran)*vph -wexe(karran)*vph**2 +weye(karran)*vph**3+be(karran)*sqrj &
            -de(karran)*sqrj*sqrj-alfe(karran)*vph*sqrj + wzero(karran)
   ELSEIF(scheme==2)THEN     ! compute the hydrogenic energies
      n=np1+jang
      Evib(np1) = -(usys/dscale(karran)**2)/(2.d0*n**2)+vzero
   ENDIF
ENDDO

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

IF(nrigid)THEN
   rotm=usys*re**2
   eneg(1)=(.5d0*we(karran)+jang*(jang+1)/(2.d0*rotm))*hsinv
   chinuj(1,1)=h(1)
ELSE
   npos=noscil*(noscil+1)/2
   DO i=1,npos
      haml(i)=0.0d0     ! initialize the hamiltonian matrix to zero
   ENDDO
   DO i=1,nhermt        ! preform the integrations
      CALL WavDer_Old(h,hd,hdn,z(i),jang,noscil,scheme,alphas)
      r=alphai*z(i)+rx
      x=r-1.0d0
      rre=r*re
      pot=potbc(rre)
      IF(ioops)THEN
         cterm=hd+sqrj*re**2/(rhoj*sin(rre/rhoj))**2+pot*hsinv
      ELSE
         cterm=hd+sqrj/r**2+pot*hsinv
      ENDIF
      wi=w(i)
      IF(medium) WRITE(Out_Unit,160) i,rre,pot,effpot
      nnp=0
      DO n=1,noscil
         term=wi*(hdn(n)+cterm)*h(n)
         DO np=1,n
            nnp=nnp+1
            haml(nnp)=haml(nnp)+h(np)*term
         ENDDO
      ENDDO
   ENDDO
   IF(medium)THEN
      WRITE(Out_Unit,*)'haml(ev)=',noscil
      DO n=1,Min(noscil,5)
         WRITE(Out_Unit,'(1x,5es15.7)')(haml(n+(np-1)*Noscil)*hs*autoev, np=1,Min(noscil,5))
      ENDDO
   ENDIF

   !-----------------------------------------------------------------------
   ! diagonalize the hamiltonian to determine the eigenvectors.
   !-----------------------------------------------------------------------
   !CALL tdiagrw(haml, eneg, vect, h, noscil, noscil, noscil*noscil)
   CALL tdiag(haml,eneg,vect,h,noscil)
   !-----------------------------------------------------------------------
   ! determine the unweighted eigenfunctions at the gaussian
   ! integration points.
   !-----------------------------------------------------------------------
   DO i=1,nhermt
      CALL WavDer_Old(h,hd,hdn,z(i),jang,noscil,scheme,alphas)
      DO j=nulow+1,nuhigh+1
         jnew=j-nulow
         sum=0.0d0
         DO k=1,noscil
            sum=sum+vect(k,j)*h(k)
         ENDDO
        chinuj(i,jnew)=sum
      ENDDO
   ENDDO
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
   IF(Evib(i)/=0.d0) perr=100.0d0*err/Evib(i)
   enegwn=enegi/cmm1toau
   inew=i-nulow
   eneg(inew)=enegi
   trans=0.0
   enegev=enegwn*cmm1toau*autoev
   IF(i>nulow+1)THEN
      trans=(eneg(inew)-eneg(inew-1))/cmm1toau
   ENDIF
   IF(little)WRITE(Out_Unit,180)im1, jang, enegwn, enegi, Evib(i), err, perr,  trans, enegev
   IF(ABS(perr)>30.0d0)THEN
     WRITE(Out_Unit,'(A,F10.3)')'Error: vibfun_Old has a problem, perr=',perr
     WRITE(Out_Unit,300)
     WRITE(Out_Unit,180) im1,jang,enegwn,enegi,Evib(i),err,perr,trans, enegev
     !STOP 'vibfun_Old'
   ENDIF
ENDDO

!
! IF desired print eigenvectors and eigenfunctions.
!-----------------------------------------------------------------------
IF(.NOT.nrigid)THEN
   IF(medium)THEN
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
   index=0
   DO i=1,nudim
      DO j=1,nhermt
         index=index+1
         haml(index)=chinuj(j,i)*sqrt(w(j))
      ENDDO
   ENDDO
   IF(medium)THEN
      WRITE(Out_Unit,210)
      CALL MxOutD(haml,nhermt,nudim,ms)
   ENDIF
ENDIF

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
ENDSUBROUTINE VibFun_Old
