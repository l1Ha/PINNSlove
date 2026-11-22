SUBROUTINE EnergyFitDAF(NStates)
USE Numeric_Kinds_Module,    ONLY: WP_Kind, IW_Kind
USE Numbers_Module,          ONLY: Zero, Ten, WP_Epsilon
USE FileUnits_Module,        ONLY: Time_Unit
USE FileUnits_OneDim_Module, ONLY: DAF_Unit, DAFbin_Unit, DAF_Oscil_Unit
USE OneDim_Parms_Module,     ONLY: nmax, MaxYterms, jmin, jmax, MinElec, MaxElec
USE QstateDAF_Module,        ONLY: EigDAF, LeftTP, RightTP, NthElec, Spin, Lambda, Nuvib, jrot, YIJ, Bound, QuasiB, Scatt, jval, nu, Elec
USE CriticalDist_Module,     ONLY: Vmax, Vinf
USE ElectronicState_Module,  ONLY: hstep, rvals, DiatomicPotential, Dipole, IthElec, KthElec
USE Convrsns_Module,         ONLY: autoev, cmm1toau
IMPLICIT NONE
CHARACTER(LEN=11), PARAMETER:: Method="Daf       "
LOGICAL, PARAMETER :: IDbug=.false., OscillatorStrengths=.False.
INTEGER NBound, NQuasi, NStates, NumNu, nYTerms, i, istate, nup, IthState, KthState
INTEGER j, mj, jp, mjp, ml
INTEGER, PARAMETER :: lval=1
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
ALLOCATE(EigDAF(0:NStates))
ALLOCATE(Psi(0:NMax,0:NStates))
REWIND(DAFBin_Unit)
Yfit=Zero
IState=0
DO WHILE(IState.lt.NStates)
   READ(DAFBin_Unit)elec,jval,numnu,nbound,nquasi,nYTerms
   WRITE(DAF_Unit,*)
   WRITE(DAF_Unit,*)"DAF States"
   WRITE(DAF_Unit,*)"Index  n  S lm Nu  j  B QB Sc    EigDAF"
   DO Nup=0,NumNu-1
      READ(DAFBin_Unit)Nu,EigDAF(IState)
      READ(DAFBin_Unit)Psi(:,IState)
      NthElec(Istate)=elec
      Spin(Istate)   =0
      Lambda(Istate) =0
      Nuvib(Istate)  =Nu
      jrot(Istate)   =jval
      IF(EigDAF(Istate)<vinf(elec))THEN
         Bound(Istate)  =.True.
         QuasiB(Istate) =.False.
         Scatt(Istate)  =.False.
      ELSEIF(EigDAF(Istate).le.vmax(elec))THEN
         Bound(Istate)  =.False.
         QuasiB(Istate) =.True.
         Scatt(Istate)  =.True.
      ELSE
         Bound(Istate)  =.False.
         QuasiB(Istate) =.False.
         Scatt(Istate)  =.True.
      ENDIF
      WRITE(DAF_Unit,'(I6,5I3,3(1x,L1),ES15.7)')Istate, NthElec(Istate), Spin(IState), Lambda(Istate), Nuvib(Istate), &
            jrot(Istate), Bound(Istate), QuasiB(Istate), Scatt(Istate), EigDAF(Istate)
      IState=IState+1
   ENDDO
   WRITE(DAF_Unit,*)
   WRITE(DAF_Unit,*)"DAF Energy Fit Parameters"
   WRITE(DAF_Unit,*)" elec   nu jval  YFit"
   READ(DAFBin_Unit)(YFit(jval,i),i=0,NYterms-1)
   WRITE(DAF_Unit,'(2x,3i5,7(ES15.7,"(Ha)"))')elec,nu,jval,(YFit(jval,i),I=0,NYTerms-1)
   WRITE(DAF_Unit,'(2x,3i5,7(ES15.7,"(cm^-1)"))')elec,nu,jval,(YFit(jval,i)/cmm1toau,I=0,NYTerms-1)
ENDDO

IF(IDbug)THEN
   WRITE(DAF_Unit,*)
   WRITE(DAF_Unit,*)"None Zero Oscillator Strengths"
   WRITE(DAF_Unit,'(2A9,A22,10A8)')'IthState','KthState','Oscillator Strength', 'IthElec', 'nu', 'j', 'mj', 'lval', 'ml', 'KthElec', 'nup', 'jp', 'mjp'
ENDIF

IF(OscillatorStrengths)THEN
WRITE(DAF_Oscil_Unit)"None Zero Oscillator Strengths"
WRITE(DAF_Oscil_Unit)'IthState','KthState','Oscillator Strength', 'IthElec', 'nu', 'j', 'mj', 'lval', 'ml', 'KthElec', 'nup', 'jp', 'mjp'
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
                     WRITE(DAF_Unit,'(2I9,ES22.14,10I8)')IthState, KthState, SUM*coef, IthElec, nu, j, mj, lval, ml, KthElec, nup, jp, mjp
                  ENDIF
                  WRITE(DAF_Oscil_Unit)IthState, KthState, SUM*coef, IthElec, nu, j, mj, lval, ml, KthElec, nup, jp, mjp
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      ENDIF
   ENDDO
ENDDO
ENDIF

ALLOCATE(LeftTP (0:NStates))
ALLOCATE(RightTP(0:NStates))
DO elec=MinElec,MaxElec
   CALL FindTurningPoints(NStates, EigDAF, Rvals, DiatomicPotential(:,elec), LeftTP, RightTP, Method, jrot, jval, NthElec)
ENDDO
elec=MinElec

DEALLOCATE(NthElec)
DEALLOCATE(Spin)
DEALLOCATE(Lambda)
DEALLOCATE(Nuvib)
DEALLOCATE(jrot)
DEALLOCATE(Bound)
DEALLOCATE(QuasiB)
DEALLOCATE(LeftTP)
DEALLOCATE(Scatt)
DEALLOCATE(RightTP)
DEALLOCATE(EigDAF)
DEALLOCATE(YFit)
DEALLOCATE(Psi)


Call CPU_TIME(cpu2)
WRITE(Time_Unit,*)'Time in EnergyFitDAF',cpu2-cpu1
CALL TimeInfoSub('Numerov','End EnergyFitDAF')
WRITE(Time_Unit,*)
 
RETURN
END SUBROUTINE EnergyFitDAF