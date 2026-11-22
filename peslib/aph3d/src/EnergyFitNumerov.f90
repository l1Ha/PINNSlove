SUBROUTINE EnergyFitNumerov(NStates)
USE Numeric_Kinds_Module,    ONLY: WP_Kind, IW_Kind
USE Numbers_Module,          ONLY: Zero, Ten, WP_Epsilon
USE FileUnits_Module,        ONLY: Time_Unit
USE FileUnits_OneDim_Module, ONLY: Numerov_Unit, Numerovbin_Unit, Numerov_Oscil_Unit
USE OneDim_Parms_Module,     ONLY: nmax, MaxYterms, jmin, jmax, MinElec, MaxElec
USE QstateNumerov_Module,    ONLY: EigNumerov, LeftTP, RightTP, NthElec, Spin, Lambda, Nuvib, jrot, YIJ, Bound, QuasiB, Scatt, jval, nu, Elec
USE CriticalDist_Module,     ONLY: Vmax, Vinf
USE ElectronicState_Module,  ONLY: hstep, rvals, DiatomicPotential, Dipole, IthElec, KthElec
USE Convrsns_Module,         ONLY: autoev, cmm1toau
IMPLICIT NONE
CHARACTER(LEN=11), PARAMETER:: Method="Numerov   "
LOGICAL, PARAMETER :: IDbug=.false., OscillatorStrengths=.False.
INTEGER, PARAMETER :: lval=1
INTEGER IthState, KthState, j, mj, jp, mjp, ml
INTEGER NBound, NQuasi, NStates, NumNu, nYTerms, i, istate, nup, nutemp
REAL(KIND=WP_Kind) cpu1, cpu2
REAL(KIND=WP_Kind) Sum, coef, Cleb
REAL(KIND=WP_Kind), ALLOCATABLE:: Yfit(:,:), Psi(:,:)

Call CPU_TIME(cpu1)
CALL TimeInfoSub('Numerov','Begin EnergyFitDAF')
Nyterms=MaxYterms
ALLOCATE(YFit(jmin:jmax,0:NYterms-1))
ALLOCATE(NthElec(0:NStates))
ALLOCATE(Spin(0:NStates))
ALLOCATE(Lambda(0:NStates))
ALLOCATE(Nuvib(0:NStates))
ALLOCATE(jrot(0:NStates))
ALLOCATE(Bound(0:NStates))
ALLOCATE(QuasiB(0:NStates))
ALLOCATE(Scatt(0:NStates))
ALLOCATE(EigNumerov(0:NStates))
ALLOCATE(Psi(0:NMax,0:NStates))
REWIND(NumerovBin_Unit)
Yfit=Zero
IState=0
elec=0
DO WHILE(IState.lt.NStates)
   WRITE(Numerov_Unit,*)
   WRITE(Numerov_Unit,*)"Numerov States"
   WRITE(Numerov_Unit,*)"Index  n  S lm Nu  j  B QB Sc    EigNumerov"
   DO Nup=0,NStates
      READ(NumerovBin_Unit)Nutemp, jval, elec
      IF(Nutemp.lt.0)Exit
      nu=nutemp
      READ(NumerovBin_Unit)EigNumerov(IState)
      READ(NumerovBin_Unit)Psi(:,IState)
      NthElec(Istate)=elec
      Spin(Istate)   =0
      Lambda(Istate) =0
      Nuvib(Istate)  =Nutemp
      jrot(Istate)   =jval
      IF(EigNumerov(Istate)<vinf(elec))THEN
         Bound(Istate)  =.True.
         QuasiB(Istate) =.False.
         Scatt(Istate)  =.False.
      ELSEIF(EigNumerov(Istate).le.vmax(elec))THEN
         Bound(Istate)  =.False.
         QuasiB(Istate) =.True.
         Scatt(Istate)  =.True.
      ELSE
         Bound(Istate)  =.False.
         QuasiB(Istate) =.False.
         Scatt(Istate)  =.True.
      ENDIF
      WRITE(Numerov_Unit,'(I6,5I3,3(1x,L1),ES15.7)')Istate, NthElec(Istate), Spin(IState), Lambda(Istate), &
            Nuvib(Istate), jrot(Istate), Bound(Istate), QuasiB(Istate), Scatt(Istate), EigNumerov(Istate)
      IState=IState+1
   ENDDO
   WRITE(Numerov_Unit,*)
   WRITE(Numerov_Unit,*)"Numerov Energy Fit Parameters"
   WRITE(Numerov_Unit,*)" elec   nu jval  YFit"
   READ(NumerovBin_Unit)elec,jval,numnu,nbound,nquasi,nyterms
   READ(NumerovBin_Unit)(YFit(jval,i),i=0,NYterms-1)
   WRITE(Numerov_Unit,'(2x,3i5,7(ES15.7,"(Ha)"))')elec,nu,jval,(YFit(jval,i),I=0,NYTerms-1)
   WRITE(Numerov_Unit,'(2x,3i5,7(ES15.7,"(cm^-1)"))')elec,nu,jval,(YFit(jval,i)/cmm1toau,I=0,NYTerms-1)
ENDDO

IF(IDbug)THEN
   WRITE(Numerov_Unit,*)
   WRITE(Numerov_Unit,*)"None Zero Oscillator Strengths"
   WRITE(Numerov_Unit,'(2A9,A22,10A8)')'IthState','KthState','Oscillator Strength', 'IthElec', 'nu', 'j', 'mj', 'lval', 'ml', 'KthElec', 'nup', 'jp', 'mjp'
ENDIF

IF(OscillatorStrengths)THEN
WRITE(Numerov_Oscil_Unit)"None Zero Oscillator Strengths"
WRITE(Numerov_Oscil_Unit)'IthState','KthState','Oscillator Strength', 'IthElec', 'nu', 'j', 'mj', 'lval', 'ml', 'KthElec', 'nup', 'jp', 'mjp'
DO IthState=0,Nstates-1
   DO KthState=IthState,Nstates-1
      IthElec=NthElec(IthState)
      KthElec=NthElec(KthState)
      SUM=Zero
      DO I=0,nmax-1
         SUM=SUM+Psi(I,IthState)*Dipole(I,IthElec,KthElec)*Psi(I,KthState)
      ENDDO
      SUM=SUM*hstep
      IF(ABS(SUM).lt.Ten*WP_Epsilon)SUM=Zero
      j= jrot(IthState)
      nu=nuvib(IthState)
      jp=jrot(KthState)
      nup=nuvib(KthState)
      IF(ABS(j-lval)<jp.and.jp<j+lval.and.sum/=Zero)THEN
      DO mj=-j,j
         DO mjp=-jp,jp
            DO ml=-lval,lval
               coef=Cleb(j,lval,jp,mj,ml,mjp)
               IF(SUM/=Zero.and.coef/=Zero)THEN
                  IF(IDbug)THEN
                     WRITE(Numerov_Unit,'(2I9,ES22.14,10I8)')IthState, KthState, SUM*coef, IthElec, nu, j, mj, lval, ml, KthElec, nup, jp, mjp
                  ENDIF
                  WRITE(Numerov_Oscil_Unit)IthState, KthState, SUM*coef, IthElec, nu, j, mj, lval, ml, KthElec, nup, jp, mjp
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      ENDIF
   ENDDO
ENDDO
ENDIF

DEALLOCATE(YFit)
DEALLOCATE(Psi)

ALLOCATE(LeftTP (0:NStates))
ALLOCATE(RightTP(0:NStates))
DO elec=MinElec,MaxElec
    CALL FindTurningPoints(NStates, EigNumerov, Rvals, DiatomicPotential(:,elec), LeftTP, RightTP, Method, jrot, jval, NthElec)
ENDDO
elec=MinElec

DEALLOCATE(LeftTP)
DEALLOCATE(RightTP)
DEALLOCATE(NthElec)
DEALLOCATE(Spin)
DEALLOCATE(Lambda)
DEALLOCATE(Nuvib)
DEALLOCATE(jrot)
DEALLOCATE(Bound)
DEALLOCATE(QuasiB)
DEALLOCATE(Scatt)
DEALLOCATE(EigNumerov)

Call CPU_TIME(cpu2)
WRITE(Time_Unit,*)'Time in EnergyFitNumerov',cpu2-cpu1
CALL TimeInfoSub('Numerov','End EnergyFitNumerov')
WRITE(Time_Unit,*)
 
RETURN
END SUBROUTINE EnergyFitNumerov