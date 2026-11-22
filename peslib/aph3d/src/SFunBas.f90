SUBROUTINE SFunBas(rho, nbasis, jtot, parity, mega, nang, marran, chanl, nangch, philast, phipr, s, v,      &
                   ThetAPH, chiaph, weight, vect, hamil, Ovrlp,  hnew, xvect, loverlap, ithrho, neigmin,    &
                   eigen, megamin, iunit, sthaph, cthaph,  nbasiss, jrot, megamax, chitau)
!

!            P U R P O S E   O F    S U B R O U T I N E
!   This routine calculates the surface functions using an analytic
!   basis set.


!            I N P U T    A R G U M E N T S
!  rho      current hyperradius.
!  jtot     total angular momentum.
!  parity   parity.
!  mega     projection of the total angular momentum on the body-fixed
!           z-axis.
!  nbasis   the total number of raw primitive basis functions.
!  nbasiss  the number of symmetry-adapted primitive basis fcns.
!           IF symmetry=.false., nbasiss=nbasis.

!            O U T P U T    A R G U M E N T S
!  philast  used to store derivatives of primitive basis.  On output
!          it contains basis of previous rho at present quad pts.
!  phipr    primitive basis at the quadrature points.
!  vect     matrix of coeffs of surface fcns expanded in prim. basis.

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
USE Numeric_Kinds_Module
USE Parms_Module
USE FileUnits_Module
USE Masses_Module
USE Converge_Module
USE Gaussb_Module
USE Convrsns_Module
USE AzBzCz_Module
USE Boundary_Module
USE Numbers_Module
IMPLICIT NONE
INTEGER jtot, parity, mega, nbasis, neigen, nang, Info, k
INTEGER marran, ibasis, iang, ithrho, mbad, jeigen, jbasis
INTEGER neigmin, megamin, iunit, nbasiss, megamax
INTEGER nangch(narran+1),  jrot(nbasis), chanl(nbasis), Potntl_Min_loc(1)
LOGICAL loverlap
REAL(KIND=WP_Kind) cpu1, cpu2, ctime, overhead
SAVE xptg, wptg, xpth, wpth, temp1, temp2, oscil
SAVE dpndx, doscdz, pn, nfac, fnorm, ScaleFac, iarran
SAVE eigenav
REAL(Kind=WP_Kind) rho, Potntl_Min
REAL(Kind=WP_Kind) ScaleFac, eigenav
REAL(Kind=WP_Kind) weight(mxang), eigenev(mxbasiss), eigen(nbasiss)
REAL(Kind=WP_Kind) Ovrlp(nbasiss,nbasiss), hamil(nbasiss,nbasiss)
REAL(Kind=WP_Kind) philast(nang, nbasis),s(nbasiss,nbasiss), v(nang)
REAL(Kind=WP_Kind) vect(nbasiss, nbasiss), phipr(nang, nbasis)
REAL(Kind=WP_Kind) ThetAPH(nang), chiaph(nang)
REAL(Kind=WP_Kind) sthaph(nang), cthaph(nang), chitau(mxang)
REAL(Kind=WP_Kind) hnew(nbasiss*nbasiss), xvect(nbasiss*nbasiss)
INTEGER, ALLOCATABLE:: iarran(:)
REAL(Kind=WP_Kind), ALLOCATABLE:: xptg(:,:), wptg(:,:)
REAL(Kind=WP_Kind), ALLOCATABLE:: xpth(:,:), wpth(:,:)
REAL(Kind=WP_Kind), ALLOCATABLE:: pn(:,:,:), dpndx(:,:,:), nfac(:)
REAL(Kind=WP_Kind), ALLOCATABLE:: oscil(:,:,:),doscdz(:,:,:), fnorm(:)
REAL(Kind=WP_Kind), ALLOCATABLE:: temp2(:), temp1(:)
REAL(Kind=WP_Kind) parityf2(mxang,narran), parityfn(mxang, narran)
REAL(Kind=WP_Kind) bfaci(mxang,narran), dzdthd(mxang,narran)
REAL(Kind=WP_Kind) dlnb(mxang,narran), dxdchi(mxang,narran)
REAL(Kind=WP_Kind) dxdtha(mxang,narran), dthddtha(mxang,narran)
REAL(Kind=WP_Kind) dthddchi(mxang,narran)
EXTERNAL popt, basis, primitiv, ovlaps, potelmnt, kintheta
EXTERNAL kinchi, calcpot
EXTERNAL rdwrsurb


!-----------------------------------------------------------------------
! Determine printing options.
!-----------------------------------------------------------------------
LOGICAL little, medium, full
INTEGER ithcall, ithsub
DATA ithcall /0/, ithsub /0/
DATA little /.false./, medium /.false./, full /.false./
CALL popt ('sfunbas ', little, medium, full, ithcall, ithsub)

IF(.NOT.ALLOCATED(xptg))THEN
   ALLOCATE(iarran(mxang))
   ALLOCATE(xptg(mxglegn, narran), wptg(mxglegn, narran))
   ALLOCATE(xpth(mxhermt, narran), wpth(mxhermt, narran))
   ALLOCATE(pn(0:mxl, mxang, narran), dpndx(0:mxl, mxang, narran), nfac(0:2*mxl))
   ALLOCATE(oscil(0:mxn, mxang, narran),doscdz(0:mxn, mxang, narran), fnorm(0:mxn))
   ALLOCATE(temp2(mxbasiss),temp1(mxbasiss))
ENDIF
!-----------------------------------------------------------------------
! Calculate the Legendre and Oscillator basis.
!-----------------------------------------------------------------------
CALL cpu_time(cpu1)
CALL cpu_time(cpu2)
overhead = cpu2 - cpu1
CALL cpu_time(cpu1)
CALL basis (xptg, wptg, xpth, wpth, dpndx, doscdz, ThetAPH, sthaph, cthaph, chiaph, pn, oscil,   &
            nang, nangch, parity, parityfn, parityf2, bfaci, dzdthd, dlnb, dxdchi, dxdtha,       &
            dthddtha, dthddchi, rho, mega, megamin, 1, weight, nfac, fnorm, chitau, iarran)
IF(rho<2.2d0)THEN
   OPEN(Unit=88, File='sfunbas_dbg.txt', Status='UNKNOWN', Position='APPEND')
   WRITE(88,'("DBG after basis: theta min/max=",2ES13.5," |sin| min=",ES13.5," nzeros=",I5," nang=",I6)') &
        MINVAL(ThetAPH), MAXVAL(ThetAPH), MINVAL(ABS(sthaph)), COUNT(ABS(sthaph)<1.d-10), nang
   CLOSE(88)
ENDIF
CALL cpu_time(cpu2)
ctime = cpu2 - cpu1 - overhead
IF(Medium)WRITE(Out_Unit,*)'   Time in sfunbas(basis) = ', ctime

!-----------------------------------------------------------------------
! Calculate the primitive basis.
!-----------------------------------------------------------------------
CALL cpu_time(cpu1)
CALL primitiv (nang, nbasis, mega, oscil, pn, parityfn, weight, phipr, wptg, wpth, bfaci)
CALL cpu_time(cpu2)
ctime = cpu2 - cpu1 - overhead
IF(Medium)WRITE(Out_Unit,*)'   Time in sfunbas(primitiv) = ', ctime

!-----------------------------------------------------------------------
! Determine overlap matrix elements since the basis is not orthonormal.
!-----------------------------------------------------------------------
CALL cpu_time(cpu1)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)"Overlap Matrix Elements in MatBas"
CALL ovlaps (nang, nbasis, chanl, nangch, phipr, Ovrlp, rho, intwt, nbasiss, jrot)
IF(rho<2.2d0)THEN
   OPEN(Unit=88, File='sfunbas_dbg.txt', Status='UNKNOWN', Position='APPEND')
   WRITE(88,'("DBG ovrlp max|.|=",ES12.4," diag min/max=",2ES12.4)') &
        MAXVAL(ABS(Ovrlp(1:nbasiss,1:nbasiss))), &
        MINVAL([(Ovrlp(ibasis,ibasis),ibasis=1,nbasiss)]), MAXVAL([(Ovrlp(ibasis,ibasis),ibasis=1,nbasiss)])
   WRITE(88,'("DBG phipr max|.|=",ES12.4)') MAXVAL(ABS(phipr(1:nang,1:nbasis)))
   CLOSE(88)
ENDIF
CALL cpu_time(cpu2)
ctime = cpu2 - cpu1 - overhead
IF(Medium)WRITE(Out_Unit,*)'   Time in sfunbas(ovlaps) = ', ctime
!-----------------------------------------------------------------------
!  Calculate the interaction potential at all of the quadrature points.
!-----------------------------------------------------------------------
CALL cpu_time(cpu1)
CALL calcpot(nang, rho, ThetAPH, chiaph, v)
CALL cpu_time(cpu2)
ctime = cpu2 - cpu1 - overhead
IF(Medium)WRITE(Out_Unit,*)'   Time in sfunbas(calcpot) = ', ctime

!-----------------------------------------------------------------------
! Calculate potential matrix elements.
!-----------------------------------------------------------------------
CALL cpu_time(cpu1)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)"Potential Matrix Elements in MatBas"
CALL potelmnt (nang, nbasis, chanl, nangch, phipr, v, hamil, rho, intwt, nbasis, nbasiss,  nbasiss, jrot)
 Potntl_Min=MinVal(V)
Potntl_Min_Loc=MinLoc(V)
CALL cpu_time(cpu2)
ctime = cpu2 - cpu1 - overhead
IF(Medium)WRITE(Out_Unit,*)'   Time in sfunbas(potelmnt) = ', ctime

!-----------------------------------------------------------------------
! Calculate the derivitive of the primitiv basis with respect to theta.
!-----------------------------------------------------------------------
CALL cpu_time(cpu1)
CALL kintheta (nang, nbasis, mega,  oscil, pn, parityfn,  weight, philast, wptg, wpth,     &
                       bfaci, dpndx, doscdz, dzdthd, dlnb, dxdtha, dthddtha)
CALL cpu_time(cpu2)
ctime = cpu2 - cpu1 - overhead
IF(Medium)WRITE(Out_Unit,*)'   Time in sfunbas(kintheta) = ', ctime

!-----------------------------------------------------------------------
! set v = 2/usys*rho**2
!-----------------------------------------------------------------------
ScaleFac=two/(usys*rho**2)
v=ScaleFac

!-----------------------------------------------------------------------
! Calculate matrix elements which have the partial with respect to
! theta.
!-----------------------------------------------------------------------
CALL cpu_time(cpu1)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)"Potential Matrix Elements (partial with respect to theta) in MatBas"
CALL potelmnt (nang, nbasis, chanl, nangch, philast, v, s, rho, intwt, nbasis, nbasiss,  nbasiss, jrot)
CALL cpu_time(cpu2)
ctime = cpu2 - cpu1 - overhead
IF(Medium)WRITE(Out_Unit,*)'   Time in sfunbas(potelmnt) = ', ctime

!-----------------------------------------------------------------------
! Accumulate Hamiltonian matrix elements in hamil
!-----------------------------------------------------------------------
CALL cpu_time(cpu1)
hamil=hamil+s
IF(rho<2.2d0)THEN
   OPEN(Unit=88, File='sfunbas_dbg.txt', Status='UNKNOWN', Position='APPEND')
   WRITE(88,'("DBG hamil max|.|=",ES12.4," (after kin theta)")') MAXVAL(ABS(hamil(1:nbasiss,1:nbasiss)))
   WRITE(88,'("DBG s     max|.|=",ES12.4," v(1)=",ES12.4)') MAXVAL(ABS(s(1:nbasiss,1:nbasiss))), v(1)
   CLOSE(88)
ENDIF
CALL cpu_time(cpu2)
ctime = cpu2 - cpu1 - overhead
IF(Medium)WRITE(Out_Unit,*)'   Time in sfunbas(Accumulate Matrix) = ', ctime

!-----------------------------------------------------------------------
! Calculate the derivative of the primitiv basis with respect to chi.
!-----------------------------------------------------------------------
CALL cpu_time(cpu1)
CALL kinchi   (nang, nbasis, mega, oscil, pn, parityfn, weight, philast, wptg, wpth,     &
               bfaci, dpndx, doscdz, dzdthd, dlnb, dxdchi, dthddchi, parityf2)
IF(rho<2.2d0)THEN
   OPEN(Unit=88, File='sfunbas_dbg.txt', Status='UNKNOWN', Position='APPEND')
   WRITE(88,'("DBG after kinchi: philast max|.|=",ES12.4," phipr max|.|=",ES12.4," nan count=",I5)') &
        MAXVAL(ABS(philast(1:nang,1:nbasis))), MAXVAL(ABS(phipr(1:nang,1:nbasis))), &
        COUNT(ISNAN(philast(1:nang,1:nbasis)))
   WRITE(88,'("DBG v after centrif: min/max=",2ES12.4)') MINVAL(v(1:nang)), MAXVAL(v(1:nang))
   CLOSE(88)
ENDIF
CALL cpu_time(cpu2)
ctime = cpu2 - cpu1 - overhead
IF(Medium)WRITE(Out_Unit,*)'   Time in sfunbas(kinchi) = ', ctime

!-----------------------------------------------------------------------
! Set v = 1/2*usys*rho**2*sin**2(ThetAPH)
!-----------------------------------------------------------------------
ScaleFac=ScaleFac/four
DO iang = 1, nang
   v(iang) = ScaleFac/(sthaph(iang))**2
ENDDO
!-----------------------------------------------------------------------
! Calculate matrix elements which have the partial with respect to
! chi.
!-----------------------------------------------------------------------
CALL cpu_time(cpu1)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)"Potential Matrix Elements (partial with respect to chi) in MatBas"
CALL potelmnt  (nang, nbasis, chanl, nangch, philast, v, s, rho, intwt, nbasis, nbasiss, nbasiss, jrot)
CALL cpu_time(cpu2)
ctime = cpu2 - cpu1 - overhead
IF(Medium)WRITE(Out_Unit,*)'   Time in sfunbas(potelmnt) = ', ctime

!-----------------------------------------------------------------------
! Accumulate Hamiltonian matrix elements in hamil
!-----------------------------------------------------------------------
CALL cpu_time(cpu1)
hamil=hamil+s
CALL cpu_time(cpu2)
ctime = cpu2 - cpu1 - overhead
IF(Medium)WRITE(Out_Unit,*)'   Time in sfunbas(Accumulate Matrix) = ', ctime

!-----------------------------------------------------------------------
! Solve the generalized eigenvalue problem.
! symmetric orthogonalization.  Russ Pack 88/9/22.
! diagonalize the overlap matrix.  Store results temporarily into
! vect and eigen.
!  -----------------------------------------------------------------
CALL cpu_time(cpu1)
!!!CALL rs(nbasiss,nbasiss,Ovrlp,eigen,1,vect,temp1,temp2,Info)
vect=Ovrlp
CALL DSYEV('V','U',nbasiss,vect,nbasiss,eigen,Ovrlp,nbasiss*nbasiss,info)
CALL cpu_time(cpu2)
ctime = cpu2 - cpu1 - overhead
IF(Medium)WRITE(Out_Unit,*)'   Time in sfunbas(DSYEV) = ', ctime

IF(Info/=0)THEN
   WRITE(Out_Unit,*)'Diagonalizing the Overlaps: Ovrlp'
   WRITE(Out_Unit,*)'Error returned from DSYEV: Info=',Info
   WRITE(Msg_Unit,*)'Diagonalizing the Overlaps: Ovrlp'
   WRITE(Msg_Unit,*)'Error returned from DSYEV: Info=',Info
   STOP 'SFunBas1'
ENDIF

IF(little)THEN
   WRITE(Out_Unit,*)' eigenvalues of the overlap matrix'
   DO k=1,min(5,max(1,nbasiss/10))
      WRITE(Out_Unit,'(1x,10es12.4)')(eigen(ibasis),ibasis=10*(k-1)+1,min(10*k,nbasiss))
   ENDDO
ENDIF
!  --------------------------------------------------------------------
!  throw out eigenvectors of overlap that have eigenvalues lt eigmin
!  form transform matrix s(nbasiss,neigen)=vect*(eigen**-1/2)
!  neigen is no of good basis fcns, mbad is no of lin. dependent ones.
!  --------------------------------------------------------------------
mbad=0   ! MBadStart
DO jbasis=1,nbasiss
   IF(eigen(jbasis)>eigmin)EXIT
   mbad=mbad+1
ENDDO
IF(rho<2.2d0)THEN
   OPEN(Unit=88, File='sfunbas_dbg.txt', Status='UNKNOWN', Position='APPEND')
   WRITE(88,'("DBG overlap lowest 12 eigen=",12ES12.3)') (eigen(jbasis),jbasis=1,12)
   CLOSE(88)
ENDIF
neigen=nbasiss-mbad
WRITE(*,'(4(A,I5),A,ES15.7)')"NBasiss=", NBasiss,"  MBad=",MBad, "  NEigen=", NEigen, "  Potntl_Min_Loc=", Potntl_Min_Loc(1), "  Potntl_Min=", Potntl_Min
WRITE(Out_Unit,'(4(A,I5),A,ES15.7)')"NBasiss=", NBasiss,"  MBad=",MBad, "  NEigen=", NEigen, "   Potntl_Min_Loc=", Potntl_Min_Loc(1), "  Potntl_Min=", Potntl_Min

DO jbasis=mbad+1,nbasiss
   eigen(jbasis)=one/sqrt(eigen(jbasis))
   jeigen=jbasis-mbad
   DO ibasis=1, nbasiss
      s(ibasis,jeigen)=vect(ibasis,jbasis)*eigen(jbasis)
   ENDDO
ENDDO

!  -------------------------------------------------------------
!  transform hamil(nbasiss,nbasiss) to good basis with rectangular
!  matrix s to make hnew(neigen,neigen).
!  use Ovrlp to temporarily store intermediate matrix.
!  -------------------------------------------------------------
CALL cpu_time(cpu1)
CALL dgemm('N', 'N', nbasiss, neigen, nbasiss,  1.d0, hamil, nbasiss, s,     nbasiss, 0.d0, Ovrlp, nbasiss)
CALL dgemm('T', 'N', neigen,  neigen, nbasiss,  1.d0, s,     nbasiss, Ovrlp, nbasiss, 0.d0, hnew,  neigen)
IF(rho<2.2d0)THEN
   OPEN(Unit=88, File='sfunbas_dbg.txt', Status='UNKNOWN', Position='APPEND')
   WRITE(88,'("DBG hnew max|.|=",ES12.4," diag min/max=",2ES12.4)') &
        MAXVAL(ABS(hnew(1:neigen*neigen))), &
        MINVAL([(hnew((ibasis-1)*neigen+ibasis),ibasis=1,neigen)]), &
        MAXVAL([(hnew((ibasis-1)*neigen+ibasis),ibasis=1,neigen)])
   WRITE(88,'("DBG hamil2 max|.|=",ES12.4," s-trans max|.|=",ES12.4)') &
        MAXVAL(ABS(hamil(1:nbasiss,1:nbasiss))), MAXVAL(ABS(s(1:nbasiss,1:neigen)))
   CLOSE(88)
ENDIF
CALL cpu_time(cpu2)
ctime = cpu2 - cpu1 - overhead
IF(Medium)WRITE(Out_Unit,*)'   Time in sfunbas(MatMul) = ', ctime

IF(rho<2.2d0)THEN
   OPEN(Unit=88, File='sfunbas_dbg.txt', Status='UNKNOWN', Position='APPEND')
   WRITE(88,'("DBG MxBasis=",I6," MxBasiss=",I6," nbasis=",I6," nbasiss=",I6," neigen=",I6)') &
        MxBasis, MxBasiss, nbasis, nbasiss, neigen
   WRITE(88,'("DBG theta min/max=",2ES13.5," |sin| min=",ES13.5," nzeros=",I5," nang=",I6)') &
        MINVAL(ThetAPH), MAXVAL(ThetAPH), MINVAL(ABS(sthaph)), &
        COUNT(ABS(sthaph)<1.d-10), nang
   CLOSE(88)
ENDIF

!  ---------------------------------------------------------------
!  diagonalize hnew and get its eigenvectors xvect and eigenvalues.
!  ---------------------------------------------------------------
CALL cpu_time(cpu1)
!!CALL rs(neigen,neigen,hnew,eigen,1,xvect,temp1,temp2,Info)
xvect=hnew
CALL DSYEV('V','U',neigen,xvect,neigen,eigen,hnew,neigen*neigen,info)
! PORT DEBUG (2026-10-01): trace the first rho to locate the -7 Ha problem.
IF(rho<2.2d0)THEN
   OPEN(Unit=88, File='sfunbas_dbg.txt', Status='UNKNOWN', Position='APPEND')
   WRITE(88,'("DBG sfunbas rho=",F7.3," potmin=",ES13.5," hamil diag min/max=",2ES13.5)') &
        rho, Potntl_Min, MINVAL([(hamil(ibasis,ibasis),ibasis=1,nbasiss)]), MAXVAL([(hamil(ibasis,ibasis),ibasis=1,nbasiss)])
   WRITE(88,'("DBG lowest 5 eigen=",5ES16.8)') (eigen(ibasis),ibasis=1,5)
   WRITE(88,'("DBG nbasiss=",I5," mbad=",I5," neigen=",I5)') nbasiss, mbad, neigen
   FLUSH(88)
   CLOSE(88)
ENDIF
IF(Potntl_Min>MinVal(Eigen))THEN
   WRITE(*,*)
   WRITE(*,*)"ERROR in SFunBas"
   WRITE(*,*)"Problem area: Rho=", rho
   WRITE(*,*)"Potntl_Min=",Potntl_Min,"  Potntl_Min_Loc=",Potntl_Min_Loc,"  EigMin=",MinVal(Eigen)
   WRITE(*,*)"Increase Eigmin:", EigMin
   WRITE(Out_Unit,*)
   WRITE(Out_Unit,*)"ERROR in SFunBas"
   WRITE(Out_Unit,*)"Problem area: Rho=", rho
   WRITE(Out_Unit,*)"Potntl_Min=",Potntl_Min,"  Potntl_Min_Loc=",Potntl_Min_Loc,"  EigMin=",MinVal(Eigen)   
   WRITE(Out_Unit,*)"Increase Eigmin:", EigMin
   STOP "SFunBas: Problem Area"
ENDIF
CALL cpu_time(cpu2)
ctime = cpu2 - cpu1 - overhead
IF(Medium)WRITE(Out_Unit,*)'   Time in sfunbas(DSYEV) = ', ctime

IF(Info/=0)THEN
   WRITE(Out_Unit,*)'Diagonalizing the Hamiltonian: hnew'
   WRITE(Out_Unit,*)'Error returned from DSYEV: Info=',Info
   WRITE(Msg_Unit,*)'Diagonalizing the Hamiltonian: hnew'
   WRITE(Msg_Unit,*)'Error returned from DSYEV: Info=',Info
   STOP 'SFunBas2'
ENDIF

!  ----------------------------------------------------------------
!  put coeffs of surface eigenfunctions relative to primitive basis
!  into vect(nbasiss,neigen)=s*xvect
!  ----------------------------------------------------------------
CALL cpu_time(cpu1)
CALL dgemm('N', 'N', nbasiss, neigen, neigen, 1.d0, s, nbasiss, xvect, neigen, 0.d0, vect, nbasiss)
CALL cpu_time(cpu2)
ctime = cpu2 - cpu1 - overhead
IF(Medium)WRITE(Out_Unit,*)'   Time in sfunbas(MatMul) = ', ctime


!  --------------------------------------------------------------
!  set down neigmin IF necessary.
!  --------------------------------------------------------------
IF(neigen<neigmin) neigmin=neigen
IF(ngood>neigmin) ngood=neigmin
!  -----------------------------------------------------------------
!  STOP IF neigmin gets too small.
!  -----------------------------------------------------------------
IF(neigmin<ngood)THEN
   WRITE(Msg_Unit,*)' neigmin too small. ','neigmin,ngood,neigen=',neigmin,ngood,neigen
   WRITE(Out_Unit,*)' neigmin too small. ','neigmin,ngood,neigen=',neigmin,ngood,neigen
   STOP 'neigmin'
ENDIF

!  --------------------------------------------------------------
!  WRITE vect(nbasiss,neigmin) on unit vec_unit ='vector'
!  presently leaves only the vect for the present, last rho value
!  done on vec_unit.
!  --------------------------------------------------------------
OPEN(Unit=vec_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Vector_ABM.bin', form='unformatted',status='unknown')
WRITE(vec_unit) nbasiss, neigmin, rho, mega
DO jbasis=1,neigmin
   WRITE(vec_unit)(vect(ibasis,jbasis), ibasis=1,nbasiss)
ENDDO
CLOSE(Unit=vec_unit)

IF(ithrho==IthStart_ABM)THEN
   OPEN(Unit=vfrst_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Vector_First.bin', form='unformatted', status='unknown')
   WRITE(vfrst_unit) nbasiss, neigmin, rho
   DO jbasis=1,neigmin
      WRITE(vfrst_unit)(vect(ibasis,jbasis), ibasis=1,nbasiss)
   ENDDO
   CLOSE(Unit=vfrst_unit)
ENDIF

!  --------------------------------------------------------------
!  CALL rdwrsurb to WRITE basis and its derivative w.o. chi.
!  WRITE surface fcn coeffs on iunit=sf33_unit or sf34_unit as
!  ithrho is even or odd.
!  --------------------------------------------------------------
CALL cpu_time(cpu1)
CALL rdwrsurb (2, iunit, mega, megamin, nang, nbasis, neigmin, phipr, philast, vect, 0, nbasis, nbasiss, nbasiss)
CALL cpu_time(cpu2)
ctime = cpu2 - cpu1 - overhead
IF(Medium)WRITE(Out_Unit,*)'   Time in sfunbas(rdwrsurf) = ', ctime

IF(ithrho==IthStart_ABM)THEN
   IF(mega==megamin)THEN
      OPEN(Unit=Pfrst_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Phiprim_First.bin',form='unformatted', status='unknown')
      REWIND Pfrst_Unit
   ENDIF
   WRITE(Pfrst_Unit) mega, nang, nbasis
   DO jbasis=1, nbasis
      WRITE(Pfrst_Unit) (phipr(iang,jbasis), iang=1, nang)
   ENDDO
   IF(mega==megamax)THEN
      CLOSE(Unit=Pfrst_Unit)
   ENDIF
ENDIF

!-----------------------------------------------------------------------
!  construct the primitive basis of the previous rho at the
!  quadrature points of the present rho.
!-----------------------------------------------------------------------
IF(loverlap.AND.ithrho/=1)THEN
   CALL cpu_time(cpu1)
   CALL basis (xptg, wptg, xpth, wpth, dpndx, doscdz,  ThetAPH, sthaph, cthaph, chiaph, pn, oscil,    &
               nang, nangch, parity, parityfn, parityf2, bfaci, dzdthd, dlnb, dxdchi, dxdtha,         &
               dthddtha, dthddchi, rho, mega, megamin, 2, weight, nfac, fnorm, chitau, iarran)
   CALL cpu_time(cpu2)
   ctime = cpu2 - cpu1 - overhead
   IF(Medium)WRITE(Out_Unit,*)'   Time in sfunbas(basis) = ', ctime

!
   CALL cpu_time(cpu1)
   CALL primitiv (nang, nbasis, mega, oscil, pn, parityfn, weight, philast, wptg, wpth, bfaci)
   CALL cpu_time(cpu2)
   ctime = cpu2 - cpu1 - overhead
   IF(Medium)WRITE(Out_Unit,*)'   Time in sfunbas(primitiv) = ', ctime

ENDIF
!  -------------------------------------------------------------------
!  get eigenvalues in electron volts.  Find average of
!  eigenvalues from nave1 to nave2.
!  -------------------------------------------------------------------
DO ibasis=1,neigmin
   eigenev(ibasis)=autoev*eigen(ibasis)
ENDDO

eigenav=zero
DO ibasis=1,ngood
   eigenav=eigenav+eigenev(ibasis)
ENDDO
eigenav=eigenav/(ngood+1)
!-----------------------------------------------------------------------
! WRITE the calculated eigenvalues in electron volts to a separate file.
!-----------------------------------------------------------------------
WRITE(sflvl_unit,'(/," $readx",5x,"distance in bohr")')
WRITE(sflvl_unit,'(1x,i5)')ione
WRITE(sflvl_unit,'(1x,f10.5)')rho
WRITE(sflvl_unit,'(" $ready",5x,"energy in ev.")')
WRITE(sflvl_unit,'(1x,i5)')neigmin
WRITE(sflvl_unit,'(1x,e14.7)')(eigenev(ibasis),ibasis=1,neigmin)
WRITE(sflvl_unit,'(" $end")')

WRITE(sflvl_unit2,'(1x,e14.7,10000(",",e14.7))')rho,(eigenev(ibasis),ibasis=1,neigmin)

!-----------------------------------------------------------------------
!  WRITE out additional information IF desired.
!  -----------------------------------------------------------

IF(Little)WRITE(Out_Unit,'(A)')'routine sfunbas'
WRITE(Out_Unit,*)
WRITE(Out_Unit,'(A,3I5)')'jtot, parity, mega = ',jtot,parity,mega
WRITE(Out_Unit,'(A,3I5,A,ES15.7)')'nang, nbasiss, neigen = ', nang,nbasiss,neigen, ' overlap eigmin = ', eigmin
WRITE(Out_Unit,'(A,3I5,ES15.7)')'ngood, nave1, nave2, eigenav=',ngood, nave1, nave2, eigenav
WRITE(Out_Unit,'(A)') ' calculated eigenvalues in eV:'
DO k=1,max(1,neigmin/10)
   WRITE(Out_Unit,'(1x,10d17.9)')(eigenev(ibasis),ibasis=10*(k-1)+1,min(10*k,neigmin))
ENDDO
WRITE(Msg_Unit,*)
WRITE(Msg_Unit,'(A,2ES15.7)')' rho, eigenav=',rho,eigenav
!WRITE(Msg_Unit,'(A,I5)')'(eigenev(ibasis),ibasis=1,min(ngood,neigmin)) with ngood=', ngood
!WRITE(Msg_Unit,'(5ES15.7)')(eigenev(ibasis),ibasis=1,min(ngood,neigmin))
WRITE(Msg_Unit,'(A,I5)')'(eigenev(ibasis),ibasis=1,min(ngood,neigmin,20)) with ngood=', ngood
WRITE(Msg_Unit,'(5ES15.7)')(eigenev(ibasis),ibasis=1,min(ngood,neigmin,20))
IF(medium.AND..NOT.full)THEN
   WRITE(Out_Unit,*)'Part of vect'
   DO ibasis = 1 , min(10,nbasiss)
      WRITE(Out_Unit,*)(vect(ibasis,jbasis),jbasis=1,min(10,nbasiss))
   ENDDO
ENDIF
IF(full)THEN
   WRITE(Out_Unit,*)'eigenvectors'
   CALL MxOut(vect, nbasiss, nbasiss)
   WRITE(Out_Unit,*)'primitive basis'
   CALL MxOut(phipr, nang, nbasis)
ENDIF

OPEN(Unit=weight_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Weights.bin',form='unformatted',status='unknown')
WRITE(weight_unit)ThetAPH
WRITE(weight_unit)chiaph
WRITE(weight_unit)weight
WRITE(weight_unit)eigen
WRITE(weight_unit)az, bz, cz
CLOSE(Unit=weight_unit,status='keep')

IF(ithrho==IthStart_ABM)THEN
   OPEN(Unit=wfrst_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Weights_First.bin',form='unformatted', status='unknown')
   WRITE(wfrst_unit)ThetAPH
   WRITE(wfrst_unit)chiaph
   WRITE(wfrst_unit)weight
   WRITE(wfrst_unit)eigen
   WRITE(wfrst_unit)az, bz, cz
   CLOSE(Unit=wfrst_unit,status='keep')
!   STOP "SfunBas TEMPMOD" !GregParker TEMPMOD
ENDIF

RETURN
ENDSUBROUTINE SFunBas
