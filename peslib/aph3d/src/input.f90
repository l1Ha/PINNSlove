SUBROUTINE input (nt, vecold, tstore, xsq, xk, rmax, chinuj, nnuj)

!-------------------------------------------------------------------
! written by g. a. parker
!this routine is called by:
!         JacBasis
!this routine calls
!         popt,qlevel,IO_Basis
!-----------------------------------------------------------------------
USE Numeric_Kinds_Module
USE fileunits_Module
USE Narran_Module
USE quantb_Module
!USE Gaussb_Module
USE GaussQuady_Module
USE pow_Module
USE Spectro_Module
USE chltot_Module
USE Storage_Module
USE Qall_Module, ONLY: xksq3, xksq3_Delves, xksq3_Jacobi, xksq3_APH, mvib3, jrot3, lorb3, mega3, nuj3, kchan, Energy3_Jacobi
USE engpro_Module
USE Arrch_Module
USE thetas_Module
USE das_Module
USE Oops_Module
USE qnumber_Module
USE melm_Module
USE tlogd_Module
USE Masses_Module
USE elemnt_Module
USE approx_Module
USE iopt_Module
USE etachn_Module
USE convrsns_Module
USE opts_Module
USE totj_Module
USE TotalEng_Module
USE Parms_Module
IMPLICIT NONE
!------------------------------------------------------------------
SAVE
LOGICAL little, medium, full, ifirst
INTEGER I, Ithcll, Ithsub, j, jmx, mmlmnt, mprnt, nchanacc, nprnt, mycall
INTEGER nt, ntt, nnujdh, nnuj(narran)
REAL(Kind=WP_Kind) chinuj(nvbrthrt)
REAL(Kind=WP_Kind) xsq(nvibrot+500), xk(nvibrot+500)
REAL(Kind=WP_Kind) tstore(nvibrot**2+500), vecold(nvibrot**2+500)
REAL(Kind=WP_Kind), ALLOCATABLE :: eint(:)
REAL(Kind=WP_Kind) potbc, rmax
DATA ifirst/.true./
DATA ithcll/0/,ithsub/0/
DATA little/.true./,medium/.false./,full/.false./
DATA mmlmnt /1/    !, mycall/0/
censud=.false.
engsud=.false.
iossud=.false.
pseudo=.false.
! mmlmnt is not currently read in, but it should be
WRITE(Out_Unit,*)'Arrangement Channel = ',karran
CALL popt('input   ',little,medium,full,ithcll,ithsub)
!-------------------------------------------------------------------
!    Etot       the total energy in hartree atomic units.
!    jtot       the total angular momentum.
!    minvib     the lowest vibrational level plus one.
!    maxvib     the highest vibrational level plus one.
! if nsymc is .true. make sure the masses are have the correct symmetry.
!-----------------------------------------------------------------------
!little=.true.
!medium=.true.
!full=.true.
IF(.NOT.ALLOCATED(Eint))ALLOCATE(Eint(nvibrot+500))
IF(nsymc(karran))THEN
IF(karran==1.AND.mass(2)/=mass(3)) GOTO 10
IF(karran==2.AND.mass(3)/=mass(1)) GOTO 10
IF(karran==3.AND.mass(1)/=mass(2)) GOTO 10
GOTO 20
   10 WRITE(Out_Unit,120) mass(1),mass(2),mass(3)
STOP 'input'
ENDIF
   20 CONTINUE
!-----------------------------------------------------------------------
! set parity=0 for even parity and set parity=1 for odd parity.
!-----------------------------------------------------------------------
IF(mmlmnt==1) melmnt='space'
IF(mmlmnt==2) melmnt='probf'
IF(mmlmnt==3) melmnt='molbf'
IF(mmlmnt<1.or.mmlmnt>3) melmnt='space'
IF(medium)THEN
WRITE(Out_Unit,180) Etot,jtot,minvib(karran),maxvib(karran)
!-------------------------------------------------------
!    jmin       an array containing the lowest rotational state
!               for each vibrational level.
!    jmax       an array containing the highest rotational state
!               for each vibrational level.
!-------------------------------------------------------
WRITE(Out_Unit,140)
WRITE(Out_Unit,200) (jmin(j,karran,0),j=minvib(karran), maxvib(karran))
WRITE(Out_Unit,190) (jmax(j,karran),j=minvib(karran), maxvib(karran))
!-------------------------------------------------------
!    mass(1)      the mass of the atom in carbon 12 mass units.
!    mass(2)      the mass of one of the atoms in the diatomic
!               molecule.
!    mass(3)      the mass of the other atom in the diatomic
!               molecule.
!-------------------------------------------------------
WRITE(Out_Unit,210) mass(1),mass(2),mass(3)
!-------------------------------------------------------
!    npar       this is true if parity decoupling is assumed.
!    nsymc       this is true if the molecule is homonuclear.
!    jeven       true for even rotational states. false for
!               odd rotation states.
!-------------------------------------------------------
WRITE(Out_Unit,220) npar,nsymc(karran),jeven
!-------------------------------------------------------
!    noscil      number of harmonic oscillators to be used in the
!                expansion of the vibrational wave functions.
!    nglegn      number of Gauss_Legendre quadrature points.
!    nhermt      number of Gauss_Hermite quadrature points.
!    ralpha      controls the curvature of the parabola.
!    rx          controls the position of the parabola.
!   npow        number of term to include in the
!                vibration expansion of the potential.
!-------------------------------------------------------
WRITE(Out_Unit,230) noscil(karran),nglegn(karran), nhermt(karran), ralpha(karran), rx(karran), npow(karran)
!-------------------------------------------------------
! scheme=1  use the harmonic oscillator basis set for the
!        vibrational asymptotic wavefunctions
! scheme=2  use the Hylleras function basis set
!-------------------------------------------------------
WRITE(Out_Unit,*)'scheme=',scheme
!-------------------------------------------------------
! the experimental coefficients for the energy. the parameters
! we,wexe,weye,be,de,alfe are in inverse centimeters and re is in bohr.
! the above parameters are defined in herzberg's book 'spectra of
! diatomic molecules.
!-------------------------------------------------------
WRITE(Out_Unit,240) re(karran),wecm(karran),wexecm(karran), weyecm(karran),becm(karran), decm(karran),alfecm(karran)
ENDIF
nrigid=.false.
IF(nhermt(karran)==1) nrigid=.true.
IF(scheme==1)THEN
  we=wecm*cmm1toau
  wexe=wexecm*cmm1toau
  weye=weyecm*cmm1toau
  be=becm*cmm1toau
  de=decm*cmm1toau
  alfe=alfecm*cmm1toau
ELSEIF(scheme==2)THEN
  we=1.0d0
ENDIF
!IF(scheme==1)THEN
!  we(karran)=wecm(karran)*cmm1toau
!  wexe(karran)=wexecm(karran)*cmm1toau
!  weye(karran)=weyecm(karran)*cmm1toau
!  be(karran)=becm(karran)*cmm1toau
!  de(karran)=decm(karran)*cmm1toau
!  alfe(karran)=alfecm(karran)*cmm1toau
!ELSEIF(scheme==2)THEN
!  we(karran)=1.0d0
!ENDIF
!-------------------------------------------------------
! compute the mass scaling factors and the system reduced mass
!-------------------------------------------------------

xktot=usys2*Etot
ioops=.false.
IF(ifirst)THEN
   ifirst=.false.
ENDIF
!-------------------------------------------------------
IF(kenergy==1)THEN
   re(karran)=re(karran)/dscale(karran)
   WRITE(*,*)'DBG input.f90 ch=',karran,' re_scaled=',re(karran),' dscale=',dscale(karran), &
             ' rhoj=',rhoj,' ioops=',ioops,' usys2=',usys2,' karran=',karran
   wzero(karran) = potbc(re(karran))
   WRITE(*,*)'DBG input.f90 ch=',karran,' potbc=',wzero(karran)
   WRITE(Out_Unit,*)'Channel=',karran,'   wzero=',wzero(karran)
ENDIF
!-------------------------------------------------------
IF(full)THEN
WRITE(Out_Unit,280) mass(1),mass(2),mass(3),usys
WRITE(Out_Unit,130) dscale
WRITE(Out_Unit,240) re(karran),we(karran), wexe(karran),weye(karran), be(karran),de(karran),alfe(karran)
!-------------------------------------------------------
!    melmnt=space  potential matrix elements are calculated
!                  in space-fixed coordinates.
!    melmnt=probf  potential matrix elements are calculated
!                  in projectile frame body-fixed coordinates.
!    melmnt=molbf  potential matrix elements are calculated
!                  in molecular frame body-fixed coordinates.
!-------------------------------------------------------
WRITE(Out_Unit,250) melmnt, iwrindep, irdindep
!-------------------------------------------------------
!    censud      true for centrifugal sudden approximations.
!    engsud      true for energy sudden approximations.
!    iossud      true for infinite order sudden approximations.
!    pseudo      pseudo sudden approximations.
!-------------------------------------------------------
WRITE(Out_Unit,260) censud,engsud,iossud,pseudo
IF(censud.AND.engsud) engsud=.false.
IF(censud.AND.iossud) censud=.false.
IF(engsud.AND.iossud) engsud=.false.
IF(pseudo) WRITE(Out_Unit,350)
IF(censud) WRITE(Out_Unit,360) lbar,mbar
IF(engsud) WRITE(Out_Unit,370) jbar,mbar
IF(iossud) WRITE(Out_Unit,380) lbar,jbar
IF(iossud.or.censud.or.engsud)THEN
IF(icase==1) WRITE(Out_Unit,150)
IF(icase==2) WRITE(Out_Unit,160)
IF(icase==3) WRITE(Out_Unit,170)
ENDIF
IF(melmnt=='space') ibody=.false.
IF(maxvib(karran)-minvib(karran)+1>noscil(karran) ) noscil(karran)=maxvib(karran)-minvib(karran)+1
WRITE(Out_Unit,310)
WRITE(Out_Unit,340) melmnt,ibody,ispace,numder
WRITE(Out_Unit,290)
WRITE(Out_Unit,320) npar,nsymc(karran),jeven, nrigid,noscil(karran),nhermt(karran),nglegn(karran),npow(karran)
ENDIF
!-------------------------------------------------------------------
!  calculate alpha
!-------------------------------------------------------------------
IF(nrigid) alpha(karran)=1.0d+00
IF(.NOT.nrigid) alpha(karran)=ralpha(karran)* sqrt(weau(karran)*usys)*re(karran)
IF(full)THEN
WRITE(Out_Unit,300)
WRITE(Out_Unit,330)ralpha(karran), alpha(karran),rx(karran),re(karran),we(karran), wexe(karran),be(karran),de(karran),alfe(karran)
ENDIF
jmx=0
DO  i=minvib(karran),maxvib(karran)
    IF(jmx<jmax(i,karran)) jmx=jmax(i,karran)
ENDDO
!-------------------------------------------------------------------
! calculate vibrational wavefunction and quantum numbers.
!-------------------------------------------------------------------
symmetry=nsymc(karran)
c(karran)=re(karran)/alpha(karran)
b(karran)=re(karran)*rx(karran)
!IF(.NOT.irdindep)THEN
nchanacc=0
nnujdh=0
DO  i=1,karran-1
  nnujdh=nnujdh+nnuj(i)/nhermt(i)
  nchanacc=nchanacc+nchanl(i)
ENDDO
CALL qlevel (nt, vecold, tstore, xsq(nchanacc+1), xk(nchanacc+1), eint(nchanacc+1), nnuj, chinuj)
!ENDIF
!write(*,*)karran  !tmpmodgregparker
!write(*,*)eint(1) !tmpmodgregparker
!write(*,*)xsq(1)  !tmpmodgregparker
!write(*,*)xk(1)   !tmpmodgregparker
!write(*,*)nnuj(1) !tmpmodgregparker
!write(*,*)nt      !tmpmodgregparker

!-------------------------------------------------------------------
! The following variables are indexed from 1 to the sumrof nt for all
! the channels, where nt is the number of
! vibrational-rotational states in
! channel karran; note that nt and nchanl(karran) are identical

! KCHAN is the label for the arrangement channel
! MVIB3 is the vibrational state
! JROT3 is the rotational state
! LORB3 is the orbital state
! MEGA3 is the z-axis projection of JTOT
! NUJ3 is an index into the chinuj array for each rot-vib wavefunction
! XKSQ3 is the wavenumber
!-------------------------------------------------------------------

!IF(iwrindep)THEN
   DO  i=1,nt
      kchan(i+nchanacc)=karran
      nuj3(i+nchanacc)=nuj(i)+nnujdh
      mvib3(i+nchanacc)=mvib(i)
      jrot3(i+nchanacc)=jthrot(i)
      lorb3(i+nchanacc)=lorb(i)
      mega3(i+nchanacc)=megaz(i)
      xksq3(i+nchanacc)=usys2*(Etot-eint(i+nchanacc))
      xsq(i+nchanacc)=xksq3(i+nchanacc)  !tmpmodGregParker
      xk(i+nchanacc)=SQRT(ABS(xsq(i+nchanacc)))    !tmpmodGregParker
   ENDDO
   !mycall=mycall+1
   !IF(mycall>=19)THEN
   !   WRITE(*,*)"MyCall=",mycall
   !ENDIF
   !IF(karran==3)THEN
   !IF(xksq3(1)/=xksq3(36))THEN
   !   write(*,*)"ts2",xksq3(1),xksq3(36)
   !   !STOP "ts2"
   !ENDIF
   !IF(xsq(1)/=xsq(36))THEN
   !   write(*,*)"ts3",xsq(1),xsq(36)
   !   !STOP "ts3"
   !ENDIF
   !IF(xk(1)/=xk(36))THEN
   !   write(*,*)"ts4",xk(1),xk(36)
   !   !STOP "ts4"
   !ENDIF
   !ENDIF
   DO i=1,nt+nchanacc
       Energy3_Jacobi(i)=eint(i)
       xksq3(i)=usys2*(Etot-eint(i))
   ENDDO
   !IF(karran==3)THEN
   !WRITE(*,'(a,e23.15,a,i5)')"etot=",etot, "  mycall=", mycall
   !WRITE(*,'(A,6e23.15)')"eint=",eint(1),eint(36),eint(71),Energy3_Jacobi(1),Energy3_Jacobi(36),Energy3_Jacobi(71)
   !WRITE(*,'(A,6e23.15)')"xsq1 =",xsq(1) ,xsq(36) ,xsq(71) ,xksq3(1)         ,xksq3(36)         ,xksq3(71)
   !   if(eint(1)/=eint(36).or.xksq3(1)/=xksq3(36))THEN
   !      stop "error in routine input1"
   !   endif
   !   if(Energy3_Jacobi(1)/=Energy3_Jacobi(36).or.xsq(1)/=xsq(36))THEN
   !      stop "error in routine input2"
   !   endif
   !ENDIF
!ENDIF
ntt=nt+nchanacc
IF(ntt>nstate)THEN
   WRITE(Out_Unit,*) ' ntt=', ntt, ' is greater than nstate=', nstate
   WRITE(Out_Unit,*) ' input: increase nstate.parms'
   WRITE(Msg_Unit,*) ' ntt=', ntt, ' is greater than nstate=', nstate
   WRITE(Msg_Unit,*) ' input: increase nstate.parms'
   STOP 'input2'
ENDIF
!Write(Out_unit,*)'Calling IO_Basis', karran, nchanl(karran),'input'
!Write(Msg_unit,*)'Calling IO_Basis', karran, nchanl(karran),'input'
CALL IO_Basis(nt,xsq,xk,eint,nchanacc) !newtmpmodgregparker
DO i=1,nt+nchanacc
       xsq(i)=xksq3(i)
ENDDO
!xsq=xksq3 !tmpmodgregparker

   !IF(karran==3)THEN
      !WRITE(*,'(A,6e23.15)')"eint=",eint(1),eint(36),eint(71),Energy3_Jacobi(1),Energy3_Jacobi(36),Energy3_Jacobi(71)
      !WRITE(*,'(A,6e23.15)')"xsq3 =",xsq(1) ,xsq(36) ,xsq(71) ,xksq3(1)         ,xksq3(36)         ,xksq3(71)
      !if(eint(1)/=eint(36).or.xksq3(1)/=xksq3(36))THEN
      !   stop "error in routine input3"
      !endif
      !if(Energy3_Jacobi(1)/=Energy3_Jacobi(36).or.xsq(1)/=xsq(36))THEN
      !   stop "error in routine input4"
      !endif
   !ENDIF
IF(full)THEN
WRITE(Out_Unit,280) mass(1),mass(2),mass(3),usys
WRITE(Out_Unit,300)
WRITE(Out_Unit,330)ralpha(karran), alpha(karran),rx(karran),re(karran),we(karran), wexe(karran),be(karran),de(karran),alfe(karran)
WRITE(Out_Unit,400) jtot
mprnt=minvib(karran)
nprnt=maxvib(karran)
WRITE(Out_Unit,390) mprnt,nprnt
WRITE(Out_Unit,420) (jmin(i,karran,0),i=minvib(karran), maxvib(karran))
WRITE(Out_Unit,410) (jmax(i,karran),i=minvib(karran), maxvib(karran))
WRITE(Out_Unit,440) nt
WRITE(Out_Unit,270) Etot
WRITE(Out_Unit,*) ' nchanacc=', nchanacc, '  jrot3='
WRITE(Out_Unit,*) (jrot3(i+nchanacc), i=1,nt)
ENDIF
RETURN
!----------------***END-input***------------------------------------

  100 FORMAT(1x,30a4)
  110 FORMAT(//,10x,30a4)
  120 FORMAT(/,'***error*** nsymc is true but the masses DO not',' have the correct symmetry',3f10.5)
  130 FORMAT(/,'dscale=',3es14.7)
  140 FORMAT(/,' jmin and jmax values for printing ')
  150 FORMAT('  exact diagonal jtot(jtot+1)+j(j+1)-2omega*omega')
  160 FORMAT('  diagonal is lorb*(lorb+1)')
  170 FORMAT('  diagonal is lbar*(lbar+1)')
  180 FORMAT(/,'Etot=',e14.7,5x,'jtot=',i5,5x,'minvib=',i5,5x,'maxvib=',i5)
  190 FORMAT(/,'jmax=',20i5)
  200 FORMAT(/,'jmin=',20i5)
  210 FORMAT(/,'mass(1)=',f20.5,5x,'mass(2)=',f20.5,5x,'mass(3)=',f20.5)
  220 FORMAT(/,'npar=',l2,5x,'nsymc=',l2,5x,'jeven=',l2)
  230 FORMAT(/,'noscil=',i5,5x,'nglegn=',i5,5x,'nhermt=',i5,5x,'ralpha=',f10.5,5x,'rx=',f10.5,'     npow=',i5)
  240 FORMAT(/,'re=',es10.3,5x,'we=',es10.3,5x,'wexe=',es10.3,5x,  &
         'weye=',es10.3,5x,'be=',es10.3,5x,'de=',es10.3,5x,'alfe=',es10.3)
  250 FORMAT(/,'melmnt=',a5,5x,'iwrindep=',l2,5x,'irdindep=',l2)
  260 FORMAT(/,'censud=',l2,5x,'engsud=',l2,5x,'iossud=',l2,5x,'pseudo=',l2)
  270 FORMAT(1x,"Etot=",e20.10)
  280 FORMAT(/,"mass(1)=",e14.7,5x,"mass(2)=",e14.7,5x,"mass(3)=",e14.7,5x, "usys=",e14.7,5x)
  290 FORMAT(/,'npar nsymc jeven nrigid noscil nhermt nglegn npow')
  300 FORMAT(/,6x,"ralpha",8x,"alpha",11x,"rx",10x,"re",11 x,"we",10x,"wexe",10x,"be",9x,"de",10x,"alfe")
  310 FORMAT(/,"melmnt ibody ispace numder")
  320 FORMAT(/,l4,l5,l5,l5,l7,i7,i7,i7,i5)
  330 FORMAT(/,10es13.6)
  340 FORMAT(/,a5,l6,l7,l7,l7,l7)
  350 FORMAT(/,'pseudo approximation for tesing only')
  360 FORMAT(/,'centrifugal sudden approximation.  lbar=',i6,'  mbar=> ',i6)
  370 FORMAT(/,'energy sudden approximation.  jbar=',i6,'  mbar=',i6)
  380 FORMAT(/,'infinite order sudden approximation.  lbar=',i6,'jbar=',i6)
  390 FORMAT(/,"vibrational states are nu=",i5,3x,"to nu=",i5)
  400 FORMAT(/,"jtot=",i5)
  410 FORMAT(/,"jmax(i)=",20i5)
  420 FORMAT(/,"jmin(i)=",20i5)
  440 FORMAT(/,"the total number of states is:",i5)
ENDSUBROUTINE input
