SUBROUTINE BoundContinuum(NStates)
USE Numeric_Kinds_Module
USE Numbers_Module
USE FileUnits_Common_Module
USE FileUnits_OneDim_Module
USE OneDim_Parms_Module
USE QstateNumerov_Module
USE CriticalDist_Module
USE ElectronicState_Module
USE DiatomicPot_Module, ONLY: DiatomicPot
IMPLICIT NONE
LOGICAL, PARAMETER :: IDbug=.false.
INTEGER NStates, i, istate, IthState, Nup, numnu, nbound, nquasi
INTEGER j, mj, jp, mjp, ml, Kenergy, nutemp
INTEGER, PARAMETER :: lval=1
REAL(KIND=WP_Kind) SUM, Cleb, E, total, coef
REAL(KIND=WP_Kind), ALLOCATABLE:: Psi(:,:), PsiC(:)

OPEN(Unit=BC_Plot_Unit,file=TRIM(OutDIR)//"DiatomicPlots/"//TRIM(DiatomicPot)//"_BoundContinuum_Transistions.csv")
ALLOCATE(Psi(0:NMax,0:NStates))
ALLOCATE(PsiC(0:NMax))
ALLOCATE(EigNumerov(0:NStates))
ALLOCATE(NthElec(0:NStates))
ALLOCATE(Spin(0:NStates))
ALLOCATE(Lambda(0:NStates))
ALLOCATE(Nuvib(0:NStates))
ALLOCATE(jrot(0:NStates))
ALLOCATE(Bound(0:NStates))
ALLOCATE(QuasiB(0:NStates))
ALLOCATE(Scatt(0:NStates))
REWIND(NumerovBin_Unit)
IState=0
DO WHILE(IState.lt.NStates)
   IF(IDbug)THEN
      WRITE(Numerov_Unit,*)
      WRITE(Numerov_Unit,*)"Numerov States"
      WRITE(Numerov_Unit,*)"Index  n  S lm Nu  j  B QB Sc    EigNumerov"
   ENDIF
   DO Nup=0,NStates
      READ(NumerovBin_Unit)Nutemp, jval, elec
      IF(Nutemp<0)EXIT
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
      IF(IDbug)THEN
         WRITE(Numerov_Unit,'(I6,5I3,3(1x,L1),ES15.7)')Istate, NthElec(Istate), Spin(IState), Lambda(Istate), Nuvib(Istate), &
               jrot(Istate), Bound(Istate), QuasiB(Istate), Scatt(Istate), EigNumerov(Istate)
      ENDIF
      IState=IState+1
      IF(IState>NStates)Exit
   ENDDO
   IF(IState<NStates)THEN
      READ(NumerovBin_Unit)elec, jval, numnu, nbound, nquasi
      READ(NumerovBin_Unit)
      IF(IDbug)THEN
         WRITE(Numerov_Unit,'(2x,3i5,7ES15.7)')elec,nu,jval
      ENDIF
   ENDIF
ENDDO


IF(IDbug)THEN
   WRITE(Numerov_Unit,*)
   WRITE(Numerov_Unit,*)"Non Zero Oscillator Strengths"
   WRITE(Numerov_Unit,'(7A9,2A22)')'IthState','KEnergy', 'IthElec', 'nu', 'j', 'lval', 'jp', 'E', 'Oscillator Strength'
ENDIF
DO IthState=0,Nstates-1
   IthElec=NthElec(IthState)
   kthelec=NthElec(0)
   REWIND(ContinuumBin_Unit)
   DO Kenergy=0,NEnergy-1
      IF(Kenergy==0)THEN
         READ(ContinuumBin_Unit)
         READ(ContinuumBin_Unit)
      ENDIF
      READ(ContinuumBin_Unit)jp, E
      READ(ContinuumBin_Unit)psiC
      SUM=Zero
      DO I=0,nmax-1
         SUM=SUM+Psi(I,IthState)*Dipole(I,IthElec,KthElec)*PsiC(I)
      ENDDO
      SUM=SUM*hstep
      IF(ABS(SUM).lt.Ten*WP_Epsilon)SUM=Zero
      j= jrot(IthState)
      nu=nuvib(IthState)
      total=zero
      IF(SUM/=Zero)THEN
         DO mj=-j,j
            DO mjp=-jp,jp
               DO ml=-lval,lval
                  coef=Cleb(j,lval,jp,mj,ml,mjp)
                  total=total+sum*coef
               ENDDO
            ENDDO
         ENDDO
         IF(ABS(Total)<=Ten*WP_Epsilon)Total=zero
         IF(total.ne.Zero)THEN
            IF(IDbug)THEN
               WRITE(Numerov_Unit,'(7I9,2ES22.14)')IthState, KEnergy, IthElec, nu, j, lval, jp, E, total
            ENDIF
            WRITE(BC_Plot_Unit,'(ES22.14,",",ES22.14)')E,total
         ENDIF
      ENDIF
   ENDDO
ENDDO
CLOSE(BC_Plot_Unit)

DEALLOCATE(Psi)
DEALLOCATE(PsiC)
DEALLOCATE(EigNumerov)
DEALLOCATE(NthElec)
DEALLOCATE(Spin)
DEALLOCATE(Lambda)
DEALLOCATE(Nuvib)
DEALLOCATE(jrot)
DEALLOCATE(Bound)
DEALLOCATE(QuasiB)
DEALLOCATE(Scatt)
elec=MinElec
 
RETURN
END SUBROUTINE BoundContinuum