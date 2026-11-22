SUBROUTINE ABM_APH
USE Numeric_Kinds_Module
USE Narran_Module
USE Parms_Module
USE TotJ_Module
USE Quantb_Module
USE FileUnits_Module
USE QCase_Module
USE Converge_Module
USE SfnDis_Module
USE Opts_Module
USE Gaussb_Module
USE Boundary_Module
USE APHChl_Module
USE Qall_Module, ONLY : Energy3_APH
IMPLICIT NONE

!            P U R P O S E   O F    S U B R O U T I N E
!     This program calculates the surface functions using
!     analytic basis functions. It also calculates overlap
!     matrix elements between the surfaces functions at the
!     previous rho with the surface functions at the current
!     rho. THEN potential energy interaction matrix elements
!     are calculated.
!
!           V A R I A B L E S
!     today     CHARACTER variable containing the current date.
!     hour      CHARACTER variable containing the current time.
!     rho       current hyperradius.
!     rholast   previous hyperradius.
!     rhonext   next hyperradius.
!     rhomin    minimum hyperradius.
!     rhomax    maximum hyperradius.
!     deltarho  fractional change in rho between sectors.
!     ithrho    ith-rho
!     naph      number of surface functions to be calculated.
!     jtot      total angular momentum.
!     mega      Projection of the total angular momentum on
!               the body-fixed z-axis.
!     megamin   minimum value of mega.
!     megamax   maximum value of mega.
!     parity    parity of the surface functions.
!     lsfunc    LOGICAL variable (true implies calculate
!               surface functions).
!     loverlap  LOGICAL variable (true implies calculate
!               overlap matrix elements).
!     lmatelm   LOGICAL variable (true implies calculate
!               coupling matrix elements).
!     amatrx   amatrx, bmatrx, and s are scratch matrices of
!              dimension mxbasis x mxbasis.
!     phinow   present primitive basis at present quadrature points.
!     philast  previous rho prim basis at  present quadrature points.
!     vect     contains the coeffs of the surf fcns in the prim. basis.

CHARACTER(LEN=10) hour
CHARACTER(LEN=8) today
CHARACTER(LEN=5) curzone
INTEGER nang, maxn, maxl, iunit, iunit1, dtvalues(8), nangch(narran+1), tea, i
REAL(Kind=WP_Kind) rho, rholast, rhonext
REAL(KIND=WP_Kind) cpu1, cpu2, overhead, tset, tsfun, tovr, tmat
! 2x2 (two electronic state) surface functions, enabled by environment
! variables APH3D_TWO_STATE / APH3D_V12_SCALE / APH3D_V22_SAME
LOGICAL two_state_run, v22_same
REAL(Kind=WP_Kind) v12_scale
CHARACTER(LEN=64) envstr
INTEGER envlen, envstat

!           A L L O C A T A B L E     A R R A Y S
INTEGER,            ALLOCATABLE:: nbasis(:),    neigmin(:), nbasiss(:)
REAL(Kind=WP_Kind), ALLOCATABLE:: philast(:,:), phinow(:,:)
REAL(Kind=WP_Kind), ALLOCATABLE:: ThetAPH(:),   chiaph(:),  s(:,:)
REAL(Kind=WP_Kind), ALLOCATABLE:: vect(:,:),    amatrx(:,:)
REAL(Kind=WP_Kind), ALLOCATABLE:: bmatrx(:,:),  sthaph(:),  weight(:)
REAL(Kind=WP_Kind), ALLOCATABLE:: eigen(:,:),   v(:),       cthaph(:)
REAL(Kind=WP_Kind), ALLOCATABLE:: cmatrix(:,:), dmatrix(:,:), chitau(:)

!            E X T E R N A L S
EXTERNAL popt, sfunbas, sfunbas_2x2, ovrbas, matbas, constant, ReadABM, openabm, setbasis

!-----------------------------------------------------------------------
! Determine printing options.
!-----------------------------------------------------------------------
LOGICAL little, medium, full
INTEGER ithsub, ithcall
DATA ithcall /0/, ithsub /0/
DATA little /.false./, medium /.false./, full /.false./
!DATA little /.True./, medium /.True./, full /.True./

!-----------------------------------------------------------------------
! Out_Unit is sfunout is the standard output file.
!-----------------------------------------------------------------------
WRITE(Out_Unit,*)'Called    ABM_APH'
CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/ABM_APH.txt', form='formatted', status='unknown')

!-------------------------------------------------------------
! Determine the time and todays date.
!-------------------------------------------------------------
CALL Date_And_Time(today, hour, curzone, dtvalues)
CALL Print_Date_And_Time(Today, Hour, CurZone, DtValues, Out_Unit)
!-----------------------------------------------------------------------
! Establish all conversion factors and physical constants
!-----------------------------------------------------------------------
CALL constant

!-----------------------------------------------------------------------
! Determine overhead time for CALL timing routine.
!-----------------------------------------------------------------------
CALL cpu_time(cpu1)
CALL cpu_time(cpu2)
overhead = cpu2 - cpu1

!-----------------------------------------------------------------------
! This routine reads in all the input DATA
!-----------------------------------------------------------------------
CALL ReadABM(nrho)


ALLOCATE(nbasis(0:mxmega), neigmin(0:mxmega), nbasiss(0:mxmega))
ALLOCATE(philast(mxang, mxbasis), phinow(mxang, mxbasis))
ALLOCATE(ThetAPH(mxang), chiaph(mxang), s(mxbasiss,mxbasiss))
ALLOCATE(vect(mxbasiss, mxbasiss), amatrx(mxbasiss,mxbasiss))
ALLOCATE(bmatrx(mxbasiss, mxbasiss), sthaph(mxang), weight(mxang))
ALLOCATE(eigen(mxbasiss, 0:mxmega), v(mxang), cthaph(mxang))
ALLOCATE(cmatrix(mxbasis,mxbasis),dmatrix(mxbasis,mxbasis),chitau(mxang))
ALLOCATE(xch(mxang), wch(mxang), xth(mxang), wth(mxang))
ALLOCATE(chanl(mxbasis, 0:mxmega), nvib(mxbasis), jrot(mxbasis,0:mxmega))

CALL popt ('main    ', little, medium, full, ithcall, ithsub)

!-----------------------------------------------------------------------
! Optional two-electronic-state (2x2) surface functions.  Off by default
! (the standard single-surface chain is completely unaffected).  Set
!   APH3D_TWO_STATE=1     to use SFunBas_2x2
!   APH3D_V12_SCALE=<x>   to scale the coupling block (0 => decoupled)
!   APH3D_V22_SAME=1      to use V11 for the second state (unit test)
!-----------------------------------------------------------------------
two_state_run = .FALSE.; v22_same = .FALSE.; v12_scale = 1.d0
CALL get_environment_variable('APH3D_TWO_STATE', envstr, envlen, envstat)
IF(envstat==0 .AND. envlen>0)THEN
   IF(envstr(1:1)=='1' .OR. envstr(1:1)=='T' .OR. envstr(1:1)=='t') two_state_run=.TRUE.
ENDIF
CALL get_environment_variable('APH3D_V12_SCALE', envstr, envlen, envstat)
IF(envstat==0 .AND. envlen>0) READ(envstr(1:envlen),*) v12_scale
CALL get_environment_variable('APH3D_V22_SAME', envstr, envlen, envstat)
IF(envstat==0 .AND. envlen>0)THEN
   IF(envstr(1:1)=='1' .OR. envstr(1:1)=='T' .OR. envstr(1:1)=='t') v22_same=.TRUE.
ENDIF
IF(two_state_run) WRITE(Out_Unit,'(A,L2,A,F8.4)') &
     'APH3D_TWO_STATE active: v22_same=', v22_same, '  v12_scale=', v12_scale

!-----------------------------------------------------------------------
! OPEN the input and output files.
!-----------------------------------------------------------------------
CALL openabm(jtot)

!  ---------------------------------------------------------------
!  set megamin
!  ---------------------------------------------------------------
megamin=0
IF((-1)**(jtot+parity)/=1) megamin=1
!-----------------------------------------------------------------------
! Set up the basis
!-----------------------------------------------------------------------
CALL cpu_time(cpu1)
NAPH=0
DO mega = megamin, megamax
   CALL setbasis (nang, nbasis, maxn, maxl, nangch, nbasiss)
   neigmin(mega)=nfreq(mega)
   NAPH=NAPH+neigmin(mega)
   WRITE(Out_Unit,*)'nfreq=',nfreq(mega)
   WRITE(Out_Unit,*)'neigmin=',neigmin
ENDDO
neigmin_sum=SUM(neigmin)

Rho=RhoMin
WRITE(Out_Unit,'(1x,A,5I5)')'MxAng,NAng,MxBasis,mxbasiss,MxMega=',MxAng,NAng,MxBasis,mxbasiss,MxMega
OPEN(Unit=SBFrst_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/SetBasis_First.bin', form='unformatted')
WRITE(SBFrst_Unit)mxang,mxbasis,mxbasiss,mxmega
WRITE(SBFrst_Unit)nang,nbasis,maxn,maxl,nangch,nbasiss,neigmin
WRITE(SBFrst_Unit)minvib, jmin, maxvib, jmax, chanl, nvib, jrot, jskip
WRITE(SBFrst_Unit)weau, ralpha, re, rx, noscil, nhermt, nglegn, nlegndre, balpha, calpha, dalpha, anharm, wexeau, intwt
WRITE(SBFrst_Unit)jtot, parity, mega, megamin, megamax, symmetry, jeven
WRITE(SBFrst_Unit)rho,scheme
CLOSE(Unit=SBFrst_Unit)

CALL cpu_time(cpu2)
tset = cpu2 - cpu1 - overhead
WRITE(Out_Unit,*)'Time used in setbasis = ',tset
IF(megamax>maxl) megamax = maxl

IF(lmatch)THEN
   nrho = 1
ENDIF

DO ithrho = IthStart_ABM, IthEnd_ABM
!-----------------------------------------------------------------------
! Initialize rho and several other variables.
!-----------------------------------------------------------------------
   !CALL rhocalc(rhomin, deltarho1, deltarho2, ithrho, rho, rholast, rhonext)
   CurrentSector=ithrho   !CurrentSector=IthRho+IthStart_ABM-1
   WRITE(Msg_Unit,*)
   WRITE(Msg_Unit,'("CurrentSector=",I5," at Rho=",ES22.15,A25)')CurrentSector,Basis_Dist(CurrentSector),Basis_Type(CurrentSector)
   WRITE(Out_Unit,*)
   WRITE(Out_Unit,'("CurrentSector=",I5," at Rho=",ES22.15,A25)')CurrentSector,Basis_Dist(CurrentSector),Basis_Type(CurrentSector)
   Rho=Basis_Dist(CurrentSector)
   !IF(Rho>=17.d0)STOP "Temporary Stop"
   IF(CurrentSector==1)THEN
      RhoLast=Basis_Dist(CurrentSector)
   ELSE
      RhoLast=Basis_Dist(CurrentSector-1)
   ENDIF
   RhoNext=Basis_Dist(CurrentSector+1)


   IF(qcase.AND.rho>=qswitch)THEN
      qcase = .false.
      DO mega = megamin, megamax
         CALL setbasis (nang, nbasis, maxn, maxl, nangch, nbasiss)
         WRITE(Out_Unit,*)'neigmin=',neigmin
      ENDDO
   ENDIF

!  --------------------------------------------------------------
!  set units on which surface functions are written and READ for use
!  in calculating matrix elements.
!  --------------------------------------------------------------
   IF((ithrho/2)*2==ithrho)THEN
      iunit=sf33_unit
      iunit1 = sf34_unit
   ELSE
      iunit=sf34_unit
      iunit1 = sf33_unit
   ENDIF

!-----------------------------------------------------------------------
! Loop over surface functions for each mega.
!-----------------------------------------------------------------------
   tea=0
   DO mega = megamin, megamax
!-----------------------------------------------------------------------
! Calculate the Surface functions at rho.
!-----------------------------------------------------------------------
      IF(lsfunc)    THEN
         CALL cpu_time(cpu1)
         IF(two_state_run)THEN
            CALL sfunbas_2x2(rho, nbasis(mega), jtot, parity, mega, nang, narran, chanl(1,mega),  &
                      nangch, philast, phinow, s, v, ThetAPH, chiaph, weight, vect, amatrx,       &
                      bmatrx, cmatrix, dmatrix, loverlap, ithrho,                                 &
                      neigmin(mega), eigen(1,mega), megamin, iunit, sthaph,                       &
                      cthaph, nbasiss(mega), jrot(1,mega), megamax, chitau, v12_scale, v22_same)
         ELSE
            CALL sfunbas(rho, nbasis(mega), jtot, parity, mega, nang, narran, chanl(1,mega),      &
                      nangch, philast, phinow, s, v, ThetAPH, chiaph, weight, vect, amatrx,    &
                      bmatrx, cmatrix, dmatrix, loverlap, ithrho,                              &
                      neigmin(mega), eigen(1,mega), megamin, iunit, sthaph,                    &
                      cthaph, nbasiss(mega), jrot(1,mega), megamax, chitau)
            DO i=1,neigmin(mega)
               tea=tea+1
               Energy3_APH(tea)=eigen(i,mega)
            ENDDO
         ENDIF
         CALL cpu_time(cpu2)
         tsfun = cpu2 - cpu1 - overhead
         IF(Medium)WRITE(Out_Unit,*)'Time in surface functions = ',tsfun
         IF(.NOT.two_state_run) CALL SFPlot_ABM_Main(Nbasis, Nbasiss, Neigmin, rho, nang, ithrho, mega)
      ENDIF
!-----------------------------------------------------------------------
! Calculate overlap integrals of the surface functions at rho with
! the surface functions at the previous rho (rholast).
!-----------------------------------------------------------------------
      IF(loverlap.AND.ithrho/=IthStart_ABM.AND..NOT.two_state_run)  THEN
          CALL cpu_time(cpu1)
          CALL ovrbas (rholast, rho, nang, nbasis(mega), chanl(1,mega),   &
&                      nangch, philast, phinow, bmatrx, vect, s,          &
&                      amatrx, intwt, ithrho, nbasiss(mega),              &
&                      neigmin(mega), iunit1, jrot(1,mega))
          CALL cpu_time(cpu2)
          tovr = cpu2 - cpu1 -overhead
          IF(Medium)WRITE(Out_Unit,*)'Time in ovrbas = ', tovr
       ENDIF
   Neigmin_Save(mega)=Neigmin(mega)
   ENDDO
   IF(SUM(Neigmin)<neigmin_sum)neigmin_sum=SUM(Neigmin)
   naph=neigmin_sum

!-----------------------------------------------------------------------
! Calculate matrix elements at several rho values in this sector.
!-----------------------------------------------------------------------
   IF(lmatelm.AND..NOT.lmatch.AND..NOT.two_state_run)   THEN
       CALL cpu_time(cpu1)
       CALL matbas (rho, nbasis, jtot, parity, megamin, megamax, nang, chanl, nangch, phinow, v,  &
&                   s, ThetAPH, chiaph, eigen, deltarho2, vect, amatrx, intwt, iunit, philast,    &
&                   bmatrx, sthaph, cthaph, neigmin, nbasiss, jrot)

       CALL cpu_time(cpu2)
       tmat = cpu2 - cpu1 - overhead
       IF(Medium)WRITE(Out_Unit,*)'Time in matbas = ', tmat
   ENDIF
   PhiLast=PhiNow
!-----------------------------------------------------------------------
! loop back to the beginning unless this was already lastrho.
!-----------------------------------------------------------------------
ENDDO

!-----------------------------------------------------------------------
! In the 2x2 check mode the energy-independent storage and the
! propagation are not meaningful (the 2N surface functions have no
! asymptotic (v,j,l) labels yet) - stop once the surface functions and
! the sfun2x2_dump.txt table have been produced.
!-----------------------------------------------------------------------
IF(two_state_run)THEN
   WRITE(Out_Unit,*)'Two-state (2x2) surface function check complete.'
   FLUSH(Out_Unit)
   STOP 'SFunBas_2x2 check complete'
ENDIF


!  --------------------------------------------------------------
!  on last rho WRITE one last line on 'Ovrlp' so that propagator
!  will use all sectors.
!  --------------------------------------------------------------
IF(loverlap.AND.ithrho/=1)THEN
   WRITE(Ovr_Unit)ithrho, rho, rhonext
   WRITE(Out_Unit,*)'ithrho,rho,rhonext',ithrho, rho, rhonext
ENDIF

WRITE(Out_Unit,*)' rholast, iunit=', rholast, iunit
WRITE(Out_Unit,*)' neigmin=',(neigmin(mega), mega=megamin, megamax)

DEALLOCATE(nbasis,  neigmin, nbasiss     )
DEALLOCATE(philast, phinow               )
DEALLOCATE(ThetAPH, chiaph,  s           )
DEALLOCATE(vect,    amatrx               )
DEALLOCATE(bmatrx,  sthaph,  weight      )
DEALLOCATE(eigen,   v,       cthaph      )
DEALLOCATE(cmatrix, dmatrix, chitau      )
DEALLOCATE(xch,     wch,     xth,     wth)
DEALLOCATE(chanl,   nvib,    jrot        )


CALL Date_And_Time(today, hour, curzone, dtvalues)
CALL Print_Date_And_Time(Today, Hour, CurZone, DtValues, Out_Unit)

DO IUnit=1,99
   IF((IUnit/=Msg_Unit).AND.(IUnit/=Out_Unit))CLOSE(IUnit)
ENDDO

CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/BasisForScatt.txt',Form='FORMATTED',ACCESS='APPEND')
WRITE(Out_Unit,*)'Completed ABM_APH'
WRITE(Out_Unit,*)

RETURN

ENDSUBROUTINE ABM_APH
