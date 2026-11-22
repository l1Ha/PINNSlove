SUBROUTINE FEM
USE Numeric_Kinds_Module
Use FileUnits_Module
USE Parms_Module   !USE Narran_Module
USE point_mat_Module
USE point_msher_Module
USE point_ovr_Module
USE point_picture_Module
USE point_quad_Module
USE point_rmsphi_Module
USE point_RSF_Module
USE point_sfunc_Module
USE elemnt_Module
USE Storage_Module
USE opts_Module
USE SfnDis_Module
USE Boundary_Module
USE restart_Module
USE mdata_Module
USE main_Module
USE ngrid_Module
USE ovrextra_Module
USE rhosur_Module
USE fuzzy_Module
USE totj_Module
USE philis_Module
USE savnam_Module
USE quad_Module
USE covlp_Module
USE crays_Module
USE potzero_Module
USE convrsns_Module
USE Numbers_Module
!
!   $RCSfile: main_FEM.f,v $ $Revision: 1.15 $
!   $Date: 89/11/16 11:03:02 $
!   $State: Stable $
!
!-----------------------------------------------------------------------
!     This is the main program for the APH 3-Dimensional Reactive
!     scattering code.
!-----------------------------------------------------------------------
CHARACTER(LEN=10) today, hour, curzone
INTEGER i, imega, iunit, mxelem, ndisce, nel, neq, mxtheta, mxchi
INTEGER nfreqsav, nmega, nnc, nodemax, numnp, nwk, nwork
INTEGER Iunit1, megmin, nelem, nelemx, idt(4000)
INTEGER dtvalues(8)
REAL(Kind=WP_Kind) cpu1, cpu2, deltarho, overhead, rholast, rhonext
REAL(Kind=WP_Kind) TSfunc, rho, tmat, tmesher, tovr, trmsphi, TPlot

!      I N T E G E R S
INTEGER jtsym, irho, pmax
!DATA nidm/5/
!-----------------------------------------------------------------------
!     READ in DATA:
!     rhomin      Value of the first rho for calculating
!                 surface functions. (units of Bohr)
!     rhomax      Value of the last rho for calculating
!                 surface functions. (units of Bohr)
!     deltarho    Length of the first sector. (unit of Bohr)
!     ithrho      IF ithrho is less than or equal to 1 a new
!                 calculation will be started.  IF ithrho is
!                 greater than 1 it is assumed that ithrho
!                 steps have already been calculated
!                 and the calculations will start at the ithrho+1
!                 step.
!     lmesher     LOGICAL variable set equal to .true. IF a mesh
!                 is to be prescribed.
!     lsfunc      LOGICAL variable set equal to .true. IF surface
!                 functions are to be calculated.
!     lmatelm     LOGICAL variable set equal to .true. IF potential
!                 energy matrix elements are to be calculated.
!     loverlap    LOGICAL variable set equal to .true. IF overlap
!                 integrals are to be calculated.
!     lplot       LOGICAL variable set equal to .true. IF plots
!                 are to be made
!-------------------------------------------------------------
!     Determine the time and todays date.
!-------------------------------------------------------------
WRITE(Out_Unit,*)'Called    FEM'
CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/FEM.txt',form='formatted', status='unknown')
WRITE(Out_unit,'(" trackertool: ithrho rho tmesher tsfunc trmsphi tplot  tovr tmat numnp nel neq ENG1 ENG2 ENG3 ENG4 ENG5 ENG6")')
WRITE(Out_unit,*)'Finite Element Program'
WRITE(Out_unit,*)'$Revision: 1.15 $'
WRITE(Out_unit,*)'$Date: 89/11/16 11:03:02 $'
WRITE(Out_unit,*)'$State: Stable $'
WRITE(Out_Unit,*)'Entering FEM'
IthRho=1
!-------------------------------------------------------------
! Determine the time and todays date.
!-------------------------------------------------------------
CALL Date_And_Time(today, hour, curzone, dtvalues)
CALL Print_Date_And_Time(Today, Hour, CurZone, DtValues, Out_Unit)
CALL title
!     -------------------------
!     Get overhead for CPUTIME.
!     -------------------------

CALL cputime (cpu1)
CALL cputime (cpu2)
overhead = cpu2 - cpu1


!-----------------------------------------------------------------------
!     Establish all conversion factors and physical constants
!-----------------------------------------------------------------------
CALL constant

!-----------------------------------------------------------------------
!     READ in all the input DATA.
!-----------------------------------------------------------------------
megamax=0
CALL readfem(nrho)

!-----------------------------------------------------------------------
!     OPEN the input and output files.
!-----------------------------------------------------------------------
CALL openfem(jtot)

!-----------------------------------------------------------------------
!     determine the smallest mega and the number of megas
!-----------------------------------------------------------------------
IF(megamax>jtot)THEN
   WRITE(Out_unit,*) '***error***: megamax bigger than jtot'
   WRITE(Out_unit,*) 'megamax=',megamax,' jtot=',jtot
   STOP 'main_FEM'
ENDIF

IF(mod(jtot,2)/=0)THEN
   jtsym=1
ELSE
   jtsym=0
ENDIF

IF(jtsym+parity==1)THEN
   megamin=1
   nmega=megamax
ELSE
   megamin=0
   nmega=megamax+1
ENDIF

!-----------------------------------------------------------------------
! Determine number of theta and chi points.
!-----------------------------------------------------------------------
ntheta=nlt*(2**ndiv)*2+1
symmetry=nsymc(1)
IF(.NOT.symmetry)THEN
   nchi=(nlc(1)+nlc(2)+nlc(3))*2**(ndiv+2)+1
ELSE
   nchi=(nlc(1)+nlc(2)+nlc(3))*2**(ndiv+1)+1
ENDIF
mxelmnt = (ntheta-1)*(nchi-1)/4
mxnode = min(ntheta*nchi,nodmax)
WRITE(Out_unit,*)'mxnode=',mxnode,'  mxelmnt=',mxelmnt
ALLOCATE(IntTheta(9*mxelmnt), IntChi(9*mxelmnt))
!-----------------------------------------------------------------------
! Allocate memory for quadrature arrays.
!-----------------------------------------------------------------------
WRITE(Out_unit,*)'nq=',nq
WRITE(Out_unit,*)'nquad=',nquad
CALL mem_quad (nquad, ntheta, nchi, nq)

!-----------------------------------------------------------------------
!     This part of the code allows a restart assuming a normal ending
!     of the previous run.
!-----------------------------------------------------------------------
IF(ithrho/=1)THEN
!-----------------------------------------------------------------------
!     READ in the nel and numnp from the last rho, used for checks.
!     Warning:  not set up yet for mega>megamin!!!!
!-----------------------------------------------------------------------
   IF(megamax/=megamin)THEN
      WRITE(Out_unit,*)'not set up to restart IF megamax/=megamin', megamin, megamax
      STOP 'aamain'
   ENDIF

   READ(OldMsh_Unit)nel,nel1,numnp
   CALL mem_rmsphi (nel, numnp, nq, ntheta, nchi, "ALLOCATE")
   CALL readfrst(ithrho, nq, rhoend, id, IntTheta, thetaval, IntChi,  chivals, &
        phi, energy, nel, numnp, ntheta, nchi, nodemax,  mxelem,     &
        megamin, megamin)
   rhosurf=rhoend
   nfreqsav=nfreq
ENDIF

!     *******************
!     *       APH       *
!     * BASIS FUNCTIONS *
!     *******************


DO 21 ithrho = 1, nrho
irho=ithrho
!WRITE(Out_Unit,*)'rhomin...=',rhomin,deltarho1,deltrho2,ithrho,rhosurf,rholast,rhonext
!CALL rhocalc(rhomin, deltarho1, deltarho2, ithrho, rhosurf, rholast, rhonext)
!WRITE(Out_Unit,*)'rhomin...=',rhomin,deltarho1,deltrho2,ithrho,rhosurf,rholast,rhonext
CurrentSector=IthRho+IthStart_FEM-1
RhoSurf=Basis_Dist(CurrentSector)
RhoNext=Basis_Dist(CurrentSector+1)
WRITE(Msg_Unit,*)
WRITE(Msg_Unit,'("CurrentSector=",I5," at Rho=",ES22.15,A25)')CurrentSector,Basis_Dist(CurrentSector),Basis_Type(CurrentSector)
WRITE(Out_Unit,'("CurrentSector=",I5," at Rho=",ES22.15,A25)')CurrentSector,Basis_Dist(CurrentSector),Basis_Type(CurrentSector)
WRITE(Out_Unit,*)'rhonext,rhosurf=',rhonext,rhosurf
deltarho = rhonext - rhosurf
tolfn=tolfn1

IF(ithmax /= 0 .AND. ithrho > ithmax) GOTO 21
IF(deltarho<=zero)THEN
   WRITE(Out_unit,*)'deltarho<=0, deltarho=',deltarho
   WRITE(Out_unit,*)'ithrho=',ithrho,' irho=',irho
   WRITE(Out_unit,*)'rhosurf,rholast,rhonext=',rhosurf,rholast,rhonext
!   STOP 'aamain'
ENDIF

irho=ithrho
nel1=nel
numnp1=numnp
IF(ithrho>1)THEN
   tolfn=tolfn2
ENDIF

!----------------------------------------------------------------------
!     copy old phirms to allow it to be used while new phirms is being
!     written.
!----------------------------------------------------------------------
IF(ithrho==1)THEN
   CALL mem_rmsphi (mxelmnt, mxnode, nq, ntheta, nchi, "ALLOCATE")
ENDIF
IF(lmesher.AND.ithrho>1)THEN
   CALL savrms(phirms, ntheta, nchi, megamin, megamax)
ENDIF

DO imega=megamin,megamax
   WRITE(Msg_Unit,'("starting rho = ",f30.15," mega = ",i3)')rhosurf, mega
   nfreq = nfreqsav
   mega=imega
!-----------------------------------------------------------------------
!     Determine the Finite element mesh at rhosurf.
!-----------------------------------------------------------------------
   IF(lmesher)   THEN
      REWIND (Msher_Bin_Unit)
!----------------------------------------------------------------------
!     IF first rho initialize phirms to one.
!     phirms is the probability density and is used in generating the
!     mesh.
!----------------------------------------------------------------------
      CALL cputime(cpu1)
      IF(ithrho==1)THEN
         WRITE(Out_unit,*)'ntheta=',ntheta,'  nchi=',nchi
         CALL vsets (ntheta*nchi, phirms, 1, One)
      ENDIF
      CALL mem_msher(ntheta, nchi, mxelmnt, mxnode, 'ALLOCATE')
      mxtheta=ntheta
      mxchi=nchi
      CALL msher(rhosurf, mxelmnt, mxnode, mxtheta, mxchi, lnod, idt, scr1, scr2, iscr1, iscr2, Virtualx, nid,    &
                 beta, nodj, id, th, ch, idteq, idiv, lcor, mapr, nodj, idivo, lcoro, mapnod, sinlist, philist,   &
                 iwork33, ithrho, numnp, ndisce, nidm, neq, ndiv, nel, ntheta, nchi)
      CALL mem_msher(ntheta, nchi, mxelmnt, mxnode, 'release')
      WRITE(Msg_Unit,'( " msher done ")')
      CALL cputime (cpu2)
      tmesher = cpu2 - cpu1 - overhead
      WRITE(Out_unit,166) tmesher
 166        FORMAT(/,'tmesher=',f10.3)

   ENDIF
!-----------------------------------------------------------------------
!     Calculate the Surface functions at rhosurf.
!-----------------------------------------------------------------------
   IF(lsfunc)    THEN
!-----------------------------------------------------------------------
!     During the first CALL to sfunc nwk is zero which signals sfunc to
!     READ in the input DATA and determine the storage parameter nwk.
!-----------------------------------------------------------------------
      nwk=0
      pmax=100
      IF(nq>neq) nq=neq
      nnc=nq*(nq+1)/2

!      CALL memory('phirms',ntheta*nchi,m_phirms,'release')
      DEALLOCATE(Phirms)
      CALL cputime (cpu1)
      DO i=1,2
         REWIND (Msher_Bin_Unit)
         CALL mem_sfunc (nel, numnp, 2*ndisce, nq, nidm, neq, nwk, nnc, nwork, nlancz, 'ALLOCATE', pmax)
         WRITE(Out_unit,*)'nwk,nidm=',nwk,nidm
         IF(i==1)THEN
            nc=nq !TMPModGregParker
            ALLOCATE(ar(nnc), br(nnc), vec(nc,nc), eigv(nc), d(nc), tt(neq), w(neq))
            ALLOCATE(rtolv(nc), evc1(nc), evc2(mxnode), ww(neq))
         ENDIF
         CALL sfunc(id, th, ch, phi, yz, lm, nid, idi, beta, nel, nodj, maxa, mht, numnp, ndisce, nidm,            &
                    neq, nfreq, nwk, nnc, nq, bghaml, bghaml(1+nwk), ar, br, vec, eigv, d, tt, w, nloc, rtolv,    &
                    evc1, evc2, energy, ww, nsit, nlancz,  nblock, maxj, nwork, mxnode, pmax)
         DEALLOCATE(ar)
         !DEALLOCATE(ar, br, vec, eigv, d, tt, w)
         !DEALLOCATE(rtolv, evc1, evc2, ww)
      ENDDO
      CALL mem_msher(ntheta, nchi, mxelmnt, mxnode, 'rel2')
      CALL mem_sfunc (nel, numnp, ndisce, nq, nidm, neq, nwk, nnc, nwork, nlancz, 'release', pmax)
      WRITE(Msg_Unit,'(" sfunc:")')

      CALL cputime (cpu2)
      tsfunc = cpu2 - cpu1 - overhead
      WRITE(Out_unit,FMT=266) tsfunc
 266  FORMAT(/,'tsfunc=',f10.3)
!-----------------------------------------------------------------------
!     Calculate the probability density.
!-----------------------------------------------------------------------
      CALL cputime (cpu1)
!      CALL memory('phirms',ntheta*nchi,m_phirms,'ALLOCATE')
      ALLOCATE(PhiRMS(NTheta,NChi))
      CALL angles (thetaval, chivals, ntheta, nchi)
      CALL rmsphi (ntheta, nchi, phirms, nel, id, IntTheta, IntChi, nq, numnp, phi, thetaval, chivals)
      WRITE(Msg_Unit,'(" rmsphi:")')
      CALL cputime (cpu2)
      trmsphi = cpu2 - cpu1 - overhead
      WRITE(Out_unit,FMT=366) trmsphi
 366  FORMAT(/,'trmsphi=',f10.3)


      IF((ithrho/2)*2==ithrho)THEN
         iunit=33
      ELSE
         iunit=34
      ENDIF
      mxelem=1715               !Tempfix
      CALL rdwrsurf (id, IntTheta, thetaval, IntChi, chivals, phi, energy, nel, &
           numnp, nq, ntheta, nchi, iunit, 2,  mxnode, mxelem, imega, megmin)
!----------------------------------------------------------------------
!     Plot the surface functions IF desired.
!----------------------------------------------------------------------
      IF(lplot)THEN
         IF(nsymc(1))THEN
            nelem = nel
            nelemx = 2*nel
         ELSE
            nelem = nel
            nelemx = nel
         ENDIF

         CALL cputime(cpu1)
         CALL mem_picture (nelemx, 'ALLOCATE')
         CALL picture (id, IntTheta, IntChi, phi, phirms, igloss, kolr, center, &
              rad, nelem, numnp, nelemx, ntheta, nchi,  nidm, chivals, &
              thetaval, nq)
         CALL mem_picture (nelemx, 'release')
         WRITE(Msg_Unit,'(" picture:")')
         CALL cputime (cpu2)
         tplot = cpu2 - cpu1 - overhead
         WRITE(Out_unit,FMT=466) tplot
 466     FORMAT(/,'tplot=',f10.3)
      ENDIF
   ENDIF
!-----------------------------------------------------------------------
!     Calculate overlap integrals of the surface functions at rhosurf
!     with the surface functions at the previous rho (rholast).
!-----------------------------------------------------------------------
   IF(loverlap.AND.ithrho/=1)  THEN
      CALL waph(phi, nq, numnp)

      IF(iunit==34)THEN
         iunit1=33
      ELSE
         iunit1=34
      ENDIF

      CALL cputime (cpu1)
      lovlp = .true.
          mxelem=1715              !Tempfix
      CALL mem_ovr (numnp, nq, nquad, nel1, numnp1, ntheta, nchi, 'ALLOCATE')
      CALL ovr (id,  IntTheta,  IntChi, phi, fb,  fbp, idlast, thlast, chlast,  philast, imega, energy1, sovrlap, thetaval,         &
                 chivals, nq, rholast, rhosurf, irho, ntheta, nchi,  numnp1, numnp, nel1, nel, iunit, iunit1, megamin,    &
                 ovrlapmx, f3, IntChi, IntTheta, igrid, xpt, wht, value,  nquad, h, p, hcoef)
      CALL mem_ovr (numnp, nq, nquad, nel1, numnp1, ntheta, nchi, 'release')
      lovlp = .false.
      WRITE(Msg_Unit,'(" ovr:")')
      CALL cputime (cpu2)
      tovr = cpu2 - cpu1 - overhead
      WRITE(Out_unit,FMT=566) tovr
 566  FORMAT(/,'tovr=',f10.3)

   ENDIF
ENDDO
!-----------------------------------------------------------------------
!     Calculate matrix elements at several rho values in this sector.
!-----------------------------------------------------------------------
IF(lmatelm)   THEN
   CALL cputime (cpu1)
   CALL mem_mat (numnp, nq, nquad, 'ALLOCATE')
   CALL mat(id,  IntTheta,  IntChi, phi,  ff,  fbar, dfbar,  lambda,  energy, hamil, thetaval, chivals, nq, rholast, rhosurf, &
            rhonext, numnp, nel, ntheta, nchi, iunit, xpt, wht, value, nquad, h, p, hcoef)
   CALL mem_mat(numnp,nq,nquad,'release') !tempmod
   WRITE(Msg_Unit,'(" mat:")')
   CALL cputime (cpu2)
   tmat = cpu2 - cpu1 - overhead
   WRITE(Out_unit,666) tmat
 666     FORMAT(/,'tmat=',f10.3)

ENDIF
!-----------------------------------------------------------------------
!     loop back to the beginning unless this was already last rho.
!-----------------------------------------------------------------------
WRITE(Msg_Unit,'(" END")')
WRITE(Out_unit,823)ithrho,rhosurf,tmesher,tsfunc,trmsphi,tplot,tovr,tmat,numnp,nel,neq,(autoev*(energy(i)-vzero),i=1,6)
 823  FORMAT(1x,'trackertool:',i5,f10.5,6f10.3,3i5,6f10.3)
CALL flusher(6)
CALL flusher(Wrmod_Unit)
 21   CONTINUE

IF(loverlap)THEN
   WRITE(FEM_Ovrlp_Unit) ithrho, rhosurf ,rhonext
   WRITE(Out_unit,*)'ithrho,rho,rhonext',ithrho, rho, rhonext
ENDIF

CALL CloseFEM
CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/BasisForScatt.txt',Form='FORMATTED',ACCESS='APPEND')
WRITE(Out_Unit,*)'Completed FEM'
WRITE(Out_Unit,*)

RETURN
END
