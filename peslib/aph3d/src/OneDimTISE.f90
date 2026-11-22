SUBROUTINE OneDimTISE(MaxVibFull)
USE Numeric_Kinds_Module
USE Numbers_Module
USE CriticalDist_Module
USE FileUnits_Module
USE FileUnits_OneDim_Module
USE Qstate_Module
USE OneDim_Parms_Module
USE ElectronicState_Module
USE DiatomicPot_Module
IMPLICIT NONE

! Variables
LOGICAL, PARAMETER:: idbug=.False., Morse_Test=.False., AllThree=.True., TimeInfo=.False.
INTEGER i, jmaximum, neig, MaxVibFull
REAL(KIND=WP_Kind), PARAMETER:: r0=Zero, tol=1.d-10
REAL(KIND=WP_Kind) r, mu, h
REAL(KIND=WP_Kind) cpu1, cpu2, overhead
REAL(KIND=WP_Kind) T_PotDiatomic, T_Morse_Pot2, T_Critical, T_Morse_Pot, T_BoundDAF, T_Numerov, T_FourierGrid
REAL(KIND=WP_Kind) T_Dipole, T_EnergyFitDAF, T_EnergyFitNumerov, T_EnergyFitFourierGrid, T_NumerovProp, T_Findvasy
REAL(KIND=WP_Kind) v(0:nmax), vp(0:nmax), vpp(0:nmax), m(0:nmax), mp(0:nmax), mpp(0:nmax)
REAL(KIND=WP_Kind) psi(0:nmax,0:nmax), eigdaf(0:nmax), eignumerov(0:nmax), eigenfourier(0:nmax)
REAL(KIND=WP_Kind) kinetic(0:nmax,0:nmax), toeplitz(0:nmax), dafx(0:nmax,0:ndaf), ab(0:42,0:nmax+42), bmat(0:nmax,0:nmax)

!Text files
OPEN(Unit=Out_Unit,file=TRIM(OutDIR)//"DiatomicOut/"//TRIM(DiatomicPot)//"_OneDimTISE.txt")
OPEN(Unit=Msg_Unit,file=TRIM(OutDIR)//"DiatomicOut/"//TRIM(DiatomicPot)//"_Msg.txt")
OPEN(Unit=Time_Unit,file=TRIM(OutDIR)//"DiatomicOut/"//TRIM(DiatomicPot)//"_Time.txt")
OPEN(Unit=Numerov_Unit,file=TRIM(OutDIR)//"DiatomicOut/"//TRIM(DiatomicPot)//"_Numerov_Results.txt")
OPEN(Unit=DAF_Unit,file=TRIM(OutDIR)//"DiatomicOut/"//TRIM(DiatomicPot)//"_DAF_Results.txt")
OPEN(Unit=Numerov_DVR_Unit,file=TRIM(OutDIR)//"DiatomicOut/"//TRIM(DiatomicPot)//"_Numerov_DVR.txt")
OPEN(Unit=DAF_DVR_Unit,file=TRIM(OutDIR)//"DiatomicOut/"//TRIM(DiatomicPot)//"_DAF_DVR.txt")
OPEN(Unit=FourierGrid_DVR_Unit,file=TRIM(OutDIR)//"DiatomicOut/"//TRIM(DiatomicPot)//"_FourierGrid_DVR.txt")
OPEN(Unit=Dbug_Unit,file=TRIM(OutDIR)//"DiatomicOut/"//TRIM(DiatomicPot)//"_Debug.txt")
OPEN(Unit=SPECTRA_Unit,file=TRIM(OutDIR)//"DiatomicOut/"//TRIM(DiatomicPot)//"_Spectra.txt")
OPEN(Unit=SYSPot_Unit,file=TRIM(OutDIR)//"DiatomicOut/"//TRIM(DiatomicPot)//"_Spectra.txt")
OPEN(Unit=FourierGrid_Unit,file=TRIM(OutDIR)//"DiatomicOut/"//TRIM(DiatomicPot)//"_FourierGrid_Results.txt")

!Plotting files
OPEN(Unit=Crit_Unit,file=TRIM(OutDIR)//"DiatomicPlots/"//TRIM(DiatomicPot)//"_CriticalDistances.csv")
OPEN(Unit=Pot_Unit,file=TRIM(OutDIR)//"DiatomicPlots/"//TRIM(DiatomicPot)//"_Potential.csv")
OPEN(Unit=PsiNumerovG_Unit,File=TRIM(OutDIR)//"DiatomicPlots/"//TRIM(DiatomicPot)//"_PsiNumerov_Grid.csv")
OPEN(Unit=PsiNumerovGS_Unit,File=TRIM(OutDIR)//"DiatomicPlots/"//TRIM(DiatomicPot)//"_PsiNumerov_Grid_Shifted.csv")
OPEN(Unit=PsiNumerovD_Unit,File=TRIM(OutDIR)//"DiatomicPlots/"//TRIM(DiatomicPot)//"_PsiNumerov_DVR.csv")
OPEN(Unit=PsiNumerovDS_Unit,File=TRIM(OutDIR)//"DiatomicPlots/"//TRIM(DiatomicPot)//"_PsiNumerov_DVR_Shifted.csv")
OPEN(Unit=PsiDAFG_Unit,File=TRIM(OutDIR)//"DiatomicPlots/"//TRIM(DiatomicPot)//"_PsiDAF_Grid.csv")
OPEN(Unit=PsiDAFGS_Unit,File=TRIM(OutDIR)//"DiatomicPlots/"//TRIM(DiatomicPot)//"_PsiDAF_Grid_Shifted.csv")
OPEN(Unit=PsiDAFD_Unit,File=TRIM(OutDIR)//"DiatomicPlots/"//TRIM(DiatomicPot)//"_PsiDAF_DVR.csv")
OPEN(Unit=PsiDAFDS_Unit,File=TRIM(OutDIR)//"DiatomicPlots/"//TRIM(DiatomicPot)//"_PsiDAF_DVR_Shifted.csv")
OPEN(Unit=EigPltNumerov_Unit,FILE=TRIM(OutDIR)//"DiatomicPlots/"//TRIM(DiatomicPot)//"_EigPltNumerov.csv")
OPEN(Unit=EigPltDAF_Unit,FILE=TRIM(OutDIR)//"DiatomicPlots/"//TRIM(DiatomicPot)//"_EigPltDAF.csv")
OPEN(Unit=EigPltFourierGrid_Unit,FILE=TRIM(OutDIR)//"DiatomicPlots/"//TRIM(DiatomicPot)//"_EigPltFourierGrid.csv")
OPEN(Unit=PsiFourierGridG_Unit,File=TRIM(OutDIR)//"DiatomicPlots/"//TRIM(DiatomicPot)//"_PsiFourierGrid_Grid.csv")
OPEN(Unit=PsiFourierGridGS_Unit,File=TRIM(OutDIR)//"DiatomicPlots/"//TRIM(DiatomicPot)//"_PsiFourierGrid_Shifted.csv")
OPEN(Unit=PsiFourierGridD_Unit,File=TRIM(OutDIR)//"DiatomicPlots/"//TRIM(DiatomicPot)//"_PsiFourierGrid_DVR.csv")
OPEN(Unit=PsiFourierGridDS_Unit,File=TRIM(OutDIR)//"DiatomicPlots/"//TRIM(DiatomicPot)//"_PsiFourierGrid_DVR_Shifted.csv")


!Binary files
OPEN(Unit=DAFBin_Unit,file=TRIM(OutDIR)//"BinOut/"//TRIM(DiatomicPot)//"_DAFBin.bin",form="unformatted")
OPEN(Unit=NumerovBin_Unit,file=TRIM(OutDIR)//"BinOut/"//TRIM(DiatomicPot)//"_NumerovBin.bin",form="unformatted")
OPEN(Unit=FourierGridBin_Unit,file=TRIM(OutDIR)//"BinOut/"//TRIM(DiatomicPot)//"_FourierGridBin.bin",form="unformatted")

IF(IARRAN==1)THEN
   T_Findvasy=0
   T_PotDiatomic=0
   T_Morse_Pot2=0
   T_Critical=0
   T_Morse_Pot=0
   T_BoundDAF=0
   T_Numerov=0
   T_FourierGrid=0
   T_Dipole=0
   T_EnergyFitDAF=0
   T_EnergyFitNumerov=0
   T_EnergyFitFourierGrid=0
   T_NumerovProp=0
ENDIF

Call CPU_TIME(cpu1)
!CALL TimeInfoSub('Main_Start','Start OneDimTISE')
CALL Findvasy(vasy, ri, hstep)
Call CPU_TIME(cpu2)
T_Findvasy=T_Findvasy+cpu2-cpu1

hstep=(ri-r0)/nmax
jval=0

!call testpot(jval,mu)

ALLOCATE(rvals(0:nmax))
IF(.not.ALLOCATED(DiatomicPotential))ALLOCATE(DiatomicPotential(0:nmax, MinElec:MaxElec))
DiatomicPotential=Zero

NStatesDAF(IARRAN)=0
nstatesNumerov(IARRAN)=0
nstatesFourierGrid(IARRAN)=0
OPEN(Unit=YFitPltDaf_Unit,file=TRIM(OutDIR)//"DiatomicPlots/"//TRIM(DiatomicPot)//"_Yfit_DAF.csv")
OPEN(Unit=YFitPltNumerov_Unit,file=TRIM(OutDIR)//"DiatomicPlots/"//TRIM(DiatomicPot)//"_Yfit_Numerov.csv")
OPEN(Unit=YFitPltFourier_Unit,file=TRIM(OutDIR)//"DiatomicPlots/"//TRIM(DiatomicPot)//"_Yfit_FourierGrid.csv")
IF(IARRAN==1)THEN
   OPEN(Unit=YFitBin_Unit,file=TRIM(OutDIR)//"BinOut/"//"Yfit.Bin", form='unformatted',action='write')
   WRITE(Pot_Unit,'(3(A7,","),4(A15,","))')"iarran", "elec", "jval", "r", "v(i)", "vp(i)", "vpp(i)"
ENDIF
DO elec=MinElec,MaxElec
   DO jval=jmin,jmax
      Call CPU_TIME(cpu1)
      WRITE(Screen_Unit,*)"elec=",elec," jval=",jval
      DO i=0,nmax
         r=i*hstep+r0
         rvals(i)=r
         CALL PotDiatomic(r, hstep, v(i), vp(i), vpp(i), jval, mu)
         WRITE(Pot_Unit,'(3(I7,","),4(ES15.7,","))')iarran, elec, jval, r, v(i), vp(i), vpp(i)
      ENDDO
      DiatomicPotential(:,elec)=v
      WRITE(Pot_Unit,'("&")')
      Call CPU_TIME(cpu2)
      T_PotDiatomic=T_PotDiatomic+cpu2-cpu1
      
      Call CPU_TIME(cpu1)
      IF(Morse_Test)THEN
         CALL Morse_Pot2(rvals, nmax, mu, v, vp, vpp)
         REWIND(UNIT=Pot_Unit)
         DO i=0,nmax
            r=i*hstep+r0
            rvals(i)=r
            WRITE(Pot_Unit,'(3(I7,","),4(ES15.7,","))')iarran, elec, jval, r, v(i), vp(i), vpp(i)
         ENDDO
         WRITE(Pot_Unit,'("&")')
      ENDIF 
      Call CPU_TIME(cpu2)
      T_Morse_Pot2=T_Morse_Pot2+cpu2-cpu1
      
      Call CPU_TIME(cpu1)
      WRITE(Screen_Unit,*)"CALLING Critical"
      CALL Critical(rvals, hstep, mu, v, vp, vpp, jval, elec)
      IF((re.eq.Zero).and.(rb.eq.Zero).and.(rc.eq.Zero))EXIT
      jmaximum=jval
      Call CPU_TIME(cpu2)
      T_Critical=T_Critical+cpu2-cpu1
      
      Call CPU_TIME(cpu1)
      IF(jval.eq.0)CALL Morse_Pot(rvals, nmax, mu, m, mp, mpp)
      Call CPU_TIME(cpu2)
      T_Morse_Pot=T_Morse_Pot+cpu2-cpu1

      Call CPU_TIME(cpu1)
      IF(IDBug)Write(dbug_unit,*)'nmax=',nmax
      WRITE(Screen_Unit,*)"CALLING BoundDAF"
      CALL BoundDAF(hstep, mu, v, eigdaf, psi, kinetic, toeplitz, dafx, ab, neig, MaxVibFull)    ! Use DAF Bound State code
      NStatesDAF(IARRAN)=NStatesDAF(IARRAN)+neig
      IF(NStatesDAF(IARRAN)>0)CALL DAFToDVR(NStatesDAF(IARRAN))
      Call CPU_TIME(cpu2)
      T_BoundDAF=T_BoundDAF+cpu2-cpu1

      IF(AllThree)THEN
         Call CPU_TIME(cpu1)
         WRITE(Screen_Unit,*)"CALLING NUMEROV"
         eignumerov=eigdaf
         ! PORT FIX (2026-10-01): Numerov_Matrix computes the same bound states
         ! but never writes BinOut/<mol>_NumerovBin.bin, which NumerovToDVR,
         ! EnergyFitNumerov and CompareMethods all read (=> "End of file").
         ! The legacy propagator routine writes that file in the expected
         ! format, so use it here.
         CALL Numerov(jval, hstep, mu, v, nmax, eignumerov, psi, elec, MaxYterms, neig, MaxVibFull) ! Use Numerov Propagator
         !CALL Numerov_Matrix(nmax+1, hstep, two*mu, rvals, V, eignumerov, psi, kinetic, bmat, neig, MaxVibFull)
         NStatesNumerov(IARRAN)=NStatesNumerov(IARRAN)+neig
         IF(NStatesNumerov(IARRAN)>0)CALL NumerovToDVR(NStatesNumerov(IARRAN))
         Call CPU_TIME(cpu2)
         T_Numerov=T_Numerov+cpu2-cpu1
      
         Call CPU_TIME(cpu1)
         WRITE(Screen_Unit,*)"CALLING FourierGrid"
         CALL FourierGrid(nmax/2, v, mu, eigenfourier, Psi, Kinetic, hstep, neig, MaxVibFull) ! Use FourierGrid Method
         NStatesFourierGrid(IARRAN)=NStatesFourierGrid(IARRAN)+neig
         IF(NStatesFourierGrid(IARRAN)>0)CALL FourierGridToDVR(NStatesFourierGrid(IARRAN))
         Call CPU_TIME(cpu2)
         T_FourierGrid=T_FourierGrid+cpu2-cpu1
      ENDIF
   ENDDO
   
   WRITE(Out_Unit,*)"Total Number of DAF States=", NStatesDAF(IARRAN)
   WRITE(DAF_Unit,*)
   WRITE(DAF_Unit,*)"Total Number of DAF States=", NStatesDAF(IARRAN)
   WRITE(DAF_Unit,*)"nstates, neig=",NStatesDAF(IARRAN),neig

   IF(AllThree)THEN
      WRITE(Out_Unit,*)    "Total Number of Numerov States=", NStatesNumerov(IARRAN)
      WRITE(Numerov_Unit,*)
      WRITE(Numerov_Unit,*)"Total Number of Numerov States=", NStatesNumerov(IARRAN)
      WRITE(Numerov_Unit,*)"nstates, neig=",NStatesNumerov(IARRAN),neig
  
      WRITE(Out_Unit,*)    "Total Number of FourierGrid States=", NStatesFourierGrid(IARRAN)
      WRITE(FourierGrid_Unit,*)
      WRITE(FourierGrid_Unit,*)"Total Number of FourierGrid States=", NStatesFourierGrid(IARRAN)
      WRITE(FourierGrid_Unit,*)"nstates, neig=",NStatesFourierGrid(IARRAN),neig
   ENDIF 
ENDDO
CLOSE(YFitPltDaf_Unit)
CLOSE(YFitPltNumerov_Unit)
CLOSE(YFitPltFourier_Unit)

Call CPU_TIME(cpu1)
WRITE(Screen_Unit,*)"CALLING DipoleFunction"
ALLOCATE(Dipole(0:nmax,MinElec:MaxElec,MinElec:MaxElec))
DO IthElec=MinElec,MaxElec
   DO KthElec=MinElec,MaxElec
      CALL DipoleFunction(rvals, Dipole(:,IthElec,KthElec), nmax, IthElec, KthElec)
   ENDDO
ENDDO
Call CPU_TIME(cpu2)
T_Dipole=T_Dipole+cpu2-cpu1

Call CPU_TIME(cpu1)
IF(FitDAFEnergies)THEN
   WRITE(Screen_Unit,*)"CALLING EnergyFitDAF: NStatesDAF(IARRAN)=", NStatesDAF(IARRAN)
   WRITE(*,*)"CALLING EnergyFitDAF: NStatesDAF(IARRAN)=", NStatesDAF(IARRAN)
   IF(IArran==1)OPEN(UNIT=DAF_Oscil_Unit,File=TRIM(OutDIR)//"BinOut/DAF_Oscil.bin",FORM="unformatted")
   WRITE(DAF_Oscil_Unit)"Diatomic Fragment=",DiatomicPot
   CALL EnergyFitDAF(NStatesDAF(IARRAN))
   IF(IArran==3)CLOSE(DAF_Oscil_Unit)
ENDIF
Call CPU_TIME(cpu2)
T_EnergyFitDAF=T_EnergyFitDAF+cpu2-cpu1

IF(AllThree)THEN
Call CPU_TIME(cpu1)
IF(FitNumerovEnergies)THEN
   WRITE(Screen_Unit,*)"CALLING EnergyFitNumerov: NStatesNumerov(IARRAN)=", NStatesNumerov(IARRAN)
   WRITE(*,*)"CALLING EnergyFitNumerov: NStatesNumerov(IARRAN)=", NStatesNumerov(IARRAN)
   IF(IArran==1)OPEN(UNIT=Numerov_Oscil_Unit,File=TRIM(OutDIR)//"BinOut/Numerov_Oscil.bin",FORM="unformatted")  
   WRITE(Numerov_Oscil_Unit)"Diatomic Fragment=",DiatomicPot
   CALL EnergyFitNumerov(NStatesNumerov(IARRAN))
   IF(IArran==3)CLOSE(Numerov_Oscil_Unit)
ENDIF
Call CPU_TIME(cpu2)
T_EnergyFitNumerov=T_EnergyFitNumerov+cpu2-cpu1

Call CPU_TIME(cpu1)
IF(FitFourierGridEnergies)THEN
   WRITE(Screen_Unit,*)"CALLING EnergyFitFourierGrid: NStatesFourierGrid(IARRAN)=", NStatesFourierGrid(IARRAN)
   WRITE(*,*)"CALLING EnergyFitFourierGrid: NStatesFourierGrid(IARRAN)=", NStatesFourierGrid(IARRAN)
   IF(IArran==1)OPEN(UNIT=FourierGrid_Oscil_Unit,File=TRIM(OutDIR)//"BinOut/FourierGrid_Oscil.bin",FORM="unformatted")  
   WRITE(FourierGrid_Oscil_Unit)"Diatomic Fragment=",DiatomicPot
   CALL EnergyFitFourierGrid(NStatesFourierGrid(IARRAN))
   IF(IArran==3)CLOSE(FourierGrid_Oscil_Unit)
ENDIF
Call CPU_TIME(cpu2)
T_EnergyFitFourierGrid=T_EnergyFitFourierGrid+cpu2-cpu1
ENDIF

continuum=.true.
Call CPU_TIME(cpu1)
IF(Continuum)THEN
   WRITE(*,*)"CALLING NumerovProp"
   WRITE(Screen_Unit,*)"CALLING NumerovProp"
   h=rvals(1)-rvals(0)
   DO elec=MinElec,MinElec
      DO jval=jmin,jmax
         CALL NumerovProp(nmax, mu, h, DiatomicPotential(:,elec), Emin, Emax, NEnergy, rvals, Psi)
      ENDDO
   ENDDO
   WRITE(Screen_Unit,*)"CALLING BoundContinuum"
   CALL BoundContinuum(NStatesDAF(IARRAN))
ENDIF
Call CPU_TIME(cpu2)
T_NumerovProp=T_NumerovProp+cpu2-cpu1


DEALLOCATE(Dipole)

FLUSH(YFitBin_Unit)
IF(IARRAN==3)THEN
   CLOSE(YFitBin_Unit)
   IF(TimeInfo)THEN
      WRITE(*,'(A,F10.5)')"T_Findvasy=",T_Findvasy
      WRITE(*,'(A,F10.5)')"T_PotDiatomic=",T_PotDiatomic
      WRITE(*,'(A,F10.5)')"T_Morse_Pot2=",T_Morse_Pot2
      WRITE(*,'(A,F10.5)')"T_Critical=",T_Critical
      WRITE(*,'(A,F10.5)')"T_Morse_Pot=",T_Morse_Pot
      WRITE(*,'(A,F10.5)')"T_BoundDAF=",T_BoundDAF
      WRITE(*,'(A,F10.5)')"T_Numerov=",T_Numerov
      WRITE(*,'(A,F10.5)')"T_FourierGrid=",T_FourierGrid
      WRITE(*,'(A,F10.5)')"T_Dipole=",T_Dipole
      WRITE(*,'(A,F10.5)')"T_EnergyFitDAF=",T_EnergyFitDAF
      WRITE(*,'(A,F10.5)')"T_EnergyFitNumerov=",T_EnergyFitNumerov
      WRITE(*,'(A,F10.5)')"T_EnergyFitFourierGrid=",T_EnergyFitFourierGrid
      WRITE(*,'(A,F10.5)')"T_NumerovProp=",T_NumerovProp
   
   
      WRITE(Out_Unit,'(A,F10.5)')"T_Findvasy=",T_Findvasy
      WRITE(Out_Unit,'(A,F10.5)')"T_PotDiatomic=",T_PotDiatomic
      WRITE(Out_Unit,'(A,F10.5)')"T_Morse_Pot2=",T_Morse_Pot2
      WRITE(Out_Unit,'(A,F10.5)')"T_Critical=",T_Critical
      WRITE(Out_Unit,'(A,F10.5)')"T_Morse_Pot=",T_Morse_Pot
      WRITE(Out_Unit,'(A,F10.5)')"T_BoundDAF=",T_BoundDAF
      WRITE(Out_Unit,'(A,F10.5)')"T_Numerov=",T_Numerov
      WRITE(Out_Unit,'(A,F10.5)')"T_FourierGrid=",T_FourierGrid
      WRITE(Out_Unit,'(A,F10.5)')"T_Dipole=",T_Dipole
      WRITE(Out_Unit,'(A,F10.5)')"T_EnergyFitDAF=",T_EnergyFitDAF
      WRITE(Out_Unit,'(A,F10.5)')"T_EnergyFitNumerov=",T_EnergyFitNumerov
      WRITE(Out_Unit,'(A,F10.5)')"T_EnergyFitFourierGrid=",T_EnergyFitFourierGrid
      WRITE(Out_Unit,'(A,F10.5)')"T_NumerovProp=",T_NumerovProp
   ENDIF
ENDIF

!CALL TimeInfoSub('Main_Stop','  End OneDimTISE')

CLOSE(Unit=Crit_Unit)
CLOSE(Unit=Pot_Unit)
CLOSE(Unit=Msg_Unit)
CLOSE(Unit=Time_Unit)
CLOSE(Unit=Dbug_Unit)

CLOSE(Unit=PsiNumerovG_Unit)
CLOSE(Unit=PsiNumerovGS_Unit)
CLOSE(Unit=PsiNumerovD_Unit)
CLOSE(Unit=PsiNumerovDS_Unit)
CLOSE(Unit=EigPltNumerov_Unit)
CLOSE(Unit=PsiDAFG_Unit)
CLOSE(Unit=PsiDAFGS_Unit)
CLOSE(Unit=PsiDAFD_Unit)
CLOSE(Unit=PsiDAFDS_Unit)
CLOSE(Unit=EigPltDAF_Unit)

DEALLOCATE(rvals)
!DEALLOCATE(DiatomicPotential)

!Binary files
CLOSE(Unit=DAFBin_Unit)
CLOSE(Unit=NumerovBin_Unit)
CLOSE(Unit=FourierGridBin_Unit)

! Body of OneDimTISE

RETURN
END SUBROUTINE OneDimTISE
