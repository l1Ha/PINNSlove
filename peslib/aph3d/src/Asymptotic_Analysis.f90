SUBROUTINE Asymptotic_Analysis(NProp)
! Author: Gregory A. Parker, Department of Physics and Astronomy, University of Oklahoma
USE Numeric_Kinds_Module
USE Numbers_Module
USE FileUnits_Asymptotic_Module
USE Energy_Module
USE EDeriv_Module
USE InputFile_Module
USE ReadK_Module
USE QAsy_Numbers_Module
USE TotalEng_Module
USE Masses_Module, ONLY: usys2
IMPLICIT NONE
INTEGER NOPEN, NProp, i
LOGICAL KNonSym                                                 ! True is K-Matrix is Non-Symmetric 
LOGICAL there
REAL(Kind=WP_Kind) EnergyNow
REAL(Kind=WP_Kind),ALLOCATABLE :: K_Matx(:,:)
REAL(Kind=WP_Kind),ALLOCATABLE :: KE_Matx(:,:) 
INTEGER inunit1, outunit1, inunit2, outunit2, DtValues(8)
CHARACTER(LEN=5)  CurZone
CHARACTER(LEN=8)  Today
CHARACTER(LEN=10) Hour
CHARACTER(LEN=5) ename
CHARACTER(LEN=19)K_File
CHARACTER(LEN=20)KE_File
CHARACTER(LEN=19), PARAMETER:: ProcName='Asymptotic_Analysis'
CHARACTER(LEN=6) Print_Flag
INTEGER kenergy

IF(.not.ALLOCATED(QNumbr))ALLOCATE(QNumbr(NQuant, NProp))                ! chanl, elec, nvib, jrot, lorb
IF(.not.ALLOCATED(xksq))ALLOCATE(xksq(NProp))              ! two times the reduced mass times the quantity of the total energy minus the int
IF(.not.ALLOCATED(AsymptEnergies))ALLOCATE(AsymptEnergies(NProp))          !
IF(.not.ALLOCATED(wvec))ALLOCATE(wvec(NProp))              !
ALLOCATE(EminLoc(NProp))           !
ALLOCATE(LastLoc(NProp),ElasLoc(NProp*NProp),InelasLoc(Nprop*NProp),ReactLoc(NProp*NProp))
IF(.not.ALLOCATED(QLabels))  ALLOCATE(QLabels(NQuant))
IF(.not.ALLOCATED(ConLabels))ALLOCATE(ConLabels(NConserved))
IF(.not.ALLOCATED(Conserved))ALLOCATE(Conserved(NConserved))
!NAtoms=3
!ALLOCATE(AtomicSymbol(NAtoms))  ! Atomic symbol
!ALLOCATE(Notes(NAtoms))         ! Notes
!ALLOCATE(AtomicNumber(NAtoms))  ! Atomic number
!ALLOCATE(MassNumber(NAtoms))    ! Mass number
!ALLOCATE(Abundance(NAtoms))     ! Isotopic abundance
!ALLOCATE(AtomicWeight(NAtoms))  ! Atomic Weight
!ALLOCATE(AtomicMass(NAtoms))    ! Atomic mass

CALL PoptAsy(ProcName, Print_Flag)
WRITE(Out_Unit,*)'Called    Asymptotic_Analysis'
CLOSE(Out_Unit)

OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/Asymptotic_Analysis.txt',Form='FORMATTED')
! K-Matrix and Energy derivative of K-Matrix

KNonSym=.False.  ! IF true generate and use non-symmetric K-Matrix

CALL Date_And_Time(today, hour, curzone, dtvalues)
CALL Print_Date_And_Time(Today, Hour, CurZone, DtValues, Out_Unit)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Entering:', ProcName

inunit1=20
outunit1=21
inunit2=22
outunit2=23

! FIND energy list
INQUIRE(file=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile, exist=there)
IF(THERE)THEN
   OPEN(inunit1,File=TRIM(InputDIR)//InputFile)
   READ(inunit1,TotEnergy)
   IF(eV_Input)THEN
      Efirst=Efirst_eV/autoeV
      DeltaEng=DeltaEng_eV/autoeV
   ENDIF
   CLOSE(inunit1)
ENDIF

OPEN(Unit=outunit1,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/Kout.txt')

OPEN(Unit=KSym_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Ksym.bin',Form='unformatted', Status='unknown')
OPEN(Unit=SReal_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/SReal.bin',Form='unformatted', Status='unknown')
OPEN(Unit=SImag_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/SImag.bin',Form='unformatted', Status='unknown')
OPEN(Unit=Probil_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Prob.bin',Form='unformatted', Status='unknown')
OPEN(Unit=Phase_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Phase.bin',Form='unformatted', Status='unknown')
OPEN(Unit=PhaseTran_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/PhaseTran.bin',Form='unformatted', Status='unknown')
OPEN(Unit=PhasePltUnit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'GraphicsOut/Phase.csv',Form='formatted', Status='unknown')


IF(EDeriv)THEN
   OPEN(Unit=EigLifePltUnit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'GraphicsOut/EigLife.csv',Form='formatted', Status='unknown')
   OPEN(Unit=outunit2,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/KEout.txt')
   OPEN(Unit=KESym_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/KEsym.bin',Form='unformatted', Status='unknown')
   OPEN(Unit=SEReal_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/SEReal.bin',Form='unformatted', Status='unknown')
   OPEN(Unit=SEImag_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/SEImag.bin',Form='unformatted', Status='unknown')
   OPEN(Unit=ProbilE_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/ProbE.bin',Form='unformatted', Status='unknown')
   OPEN(Unit=PhaseE_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/PhaseE.bin',Form='unformatted', Status='unknown')
   OPEN(Unit=PhaseTranE_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/PhaseTranE.bin',Form='unformatted', Status='unknown')
   OPEN(Unit=QReal_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/QReal.bin',Form='unformatted', Status='unknown')
   OPEN(Unit=QImag_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/QImag.bin',Form='unformatted', Status='unknown')
   OPEN(Unit=TDelayReal_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/TDelayReal.bin',Form='unformatted', Status='unknown')
   OPEN(Unit=TDelayImag_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/TDelayImag.bin',Form='unformatted', Status='unknown')
   OPEN(Unit=PhaseEPltUnit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'GraphicsOut/PhaseE.csv',Form='formatted', Status='unknown')
   OPEN(Unit=EigLifePltUnit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'GraphicsOut/EigLife.csv',Form='formatted', Status='unknown')
ENDIF


K_File='BinOut/K_Matrix.bin'
KE_File='BinOut/KE_Matrix.bin'
WRITE(Msg_Unit,'(2A,5x,A)')'InFile Read=',TRIM(K_File), TRIM(KE_File)
WRITE(Out_unit,'(2A,5x,A)')'InFile Read=',TRIM(K_File), TRIM(KE_File)
DO kenergy=1,nenergy
   IthEnergy=kenergy
   IF(kenergy<10)THEN
      WRITE(ename,'(i1)') kenergy
   ELSEIF(kenergy<100)THEN
      WRITE(ename,'(i2)') kenergy
   ELSEIF(kenergy<1000)THEN
      WRITE(ename,'(i3)') kenergy
   ELSEIF(kenergy<10000)THEN
      WRITE(ename,'(i4)') kenergy
   ELSEIF(kenergy<100000)THEN
      WRITE(ename,'(i5)') kenergy
   ELSE
      WRITE(Msg_Unit,*)'kenergy to large=',kenergy
      STOP 'main'
   ENDIF
   IF(kenergy==1)OPEN(Unit=inunit1,File=OutDIR(1:LEN(TRIM(OutDIR)))//K_File,form='unformatted',status='old') 
   IF(kenergy==1.AND.Ederiv)OPEN(Unit=inunit2,File=OutDIR(1:LEN(TRIM(OutDIR)))//KE_File,form='unformatted',status='old')
   CALL Read_K1(NOPEN,inunit1,outunit1)
   IF(Kenergy==nenergy)THEN
      Write(*,*)"nopen(Asymptotic_Analysis)=",nopen
   ENDIF
   ALLOCATE(K_Matx(nopen,nopen), KE_Matx(NOpen,NOpen))
   CALL Read_K2(NOPEN,k_matx,inunit1,outunit1)
   IF(EDeriv)THEN
      CALL Read_K1(NOPEN,inunit2,outunit2) 
      CALL Read_K2(NOPEN,ke_matx,inunit2,outunit2)
   ENDIF
   DO i=1,nopen
      xksq(i)=usys2*(Etot-AsymptEnergies(i))   ! PORT 2026-10-03: AsymptEnergies already relative (shifted in Kmatrx_Main)
   ENDDO
   CALL Print_QNumbr(NOpen, QNumbr, AsymptEnergies, xksq, wvec, NQuant,'Asymptotic_Analysis')
   energynow=efirst+deltaeng*(kenergy-1)
   IF(.not.KNonSym)THEN
      K_Matx=Half*(K_Matx+Transpose(K_Matx))
      IF(EDeriv)THEN
         KE_Matx=Half*(KE_Matx+Transpose(KE_Matx))
      ENDIF
   ENDIF
   CALL Asymptotic(NOpen, K_Matx, KE_Matx, EDeriv, KNonSym)        ! Perform asymptotic analysis
   DEALLOCATE(K_Matx, KE_Matx)
ENDDO

CALL Date_And_Time(today, hour, curzone, dtvalues)
CALL Print_Date_And_Time(Today, Hour, CurZone, DtValues, Out_Unit)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Leaving:', ProcName

CLOSE(inunit1)
CLOSE(inunit2)
CLOSE(outunit1)
CLOSE(outunit2)
CLOSE(Out_Unit)

CLOSE(KSym_Unit)
CLOSE(KESym_Unit)
CLOSE(SReal_Unit)
CLOSE(SImag_Unit)
CLOSE(SEReal_Unit)
CLOSE(SEImag_Unit)
CLOSE(Probil_Unit)
CLOSE(ProbilE_Unit)
CLOSE(Phase_Unit)
CLOSE(PhaseE_Unit)
CLOSE(PhaseTran_Unit)
CLOSE(PhaseTranE_Unit)
CLOSE(QReal_Unit)
CLOSE(QImag_Unit)
CLOSE(TDelayReal_Unit)
CLOSE(TDelayImag_Unit)
CLOSE(Unit=PhasePltUnit)
CLOSE(Unit=PhaseEPltUnit)
CLOSE(Unit=EigLifePltUnit)

!OPEN(unit=3037,File=TRIM(OutDIR)//'Output/'//'3037.txt',Form='formatted', Status='unknown')
OPEN(Unit=AP_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/Asymptotic_Plotting.txt')

CALL Asymptotic_Plotting('K_Matrix', 'K_Matrix', 'K', 'Delves', 'Delves')
CALL Asymptotic_Plotting('SReal', 'SReal', 'S_Real', 'Delves', 'Delves')
CALL Asymptotic_Plotting('SImag', 'SImag', 'S_Imag', 'Delves', 'Delves')
CALL Asymptotic_Plotting('Prob', 'Prob', 'Probability', 'Delves', 'Delves')
!CALL Asymptotic_Plotting('Phase', 'Phase', 'Phase', 'Index ', 'Delves')
CALL Asymptotic_Plotting('PhaseTran', 'PhaseTran', 'PhaseTran', 'Index ', 'Delves')
CALL Asymptotic_Plotting('SReal', 'SImag', 'S_Argand', 'Delves', 'Delves')

IF(EDeriv)THEN
   CALL Asymptotic_Plotting('KE_Matrix', 'KE_Matrix', 'KE', 'Delves', 'Delves')
   CALL Asymptotic_Plotting('SEReal', 'SEReal', 'SE_Real', 'Delves', 'Delves')
   CALL Asymptotic_Plotting('SEImag', 'SEImag', 'SE_Imag', 'Delves', 'Delves')
   CALL Asymptotic_Plotting('ProbE', 'ProbE', 'ProbabilityE', 'Delves', 'Delves')
   !CALL Asymptotic_Plotting('PhaseE', 'PhaseE', 'PhaseE', 'Index ', 'Delves')
   CALL Asymptotic_Plotting('PhaseTranE', 'PhaseTranE', 'PhaseTranE', 'Index ', 'Delves')
   CALL Asymptotic_Plotting('SEReal', 'SEImag', 'S_ArgandE', 'Delves', 'Delves')
   CALL Asymptotic_Plotting('QReal', 'QReal', 'Q_Real', 'Delves', 'Delves')
   CALL Asymptotic_Plotting('QImag', 'QImag', 'Q_Imag', 'Delves', 'Delves')   
   CALL Asymptotic_Plotting('QReal', 'QImag', 'Q_Argand', 'Delves', 'Delves')
   CALL Asymptotic_Plotting('TDelayReal', 'TDelayReal', 'TDelayReal_Real', 'Delves', 'Delves')
   CALL Asymptotic_Plotting('TDelayImag', 'TDelayImag', 'TDelayImag', 'Delves', 'Delves') 
   CALL Asymptotic_Plotting('TDelayReal', 'TDelayImag', 'TDelay_Argand', 'Delves', 'Delves')
ENDIF

!TmpModGregParker probably move elsewhere
!DEALLOCATE(QNumbr)   ! elec, chanl, nvib, jrot, lorb
!DEALLOCATE(xksq)     ! two times the reduced mass times the quantity of the total energy minus the int
!DEALLOCATE(AsymptEnergies) ! Internal energies 
!DEALLOCATE(wvec)     ! wavevector
!DEALLOCATE(EminLoc, LastLoc,ElasLoc, InelasLoc, ReactLoc)

CLOSE(AP_Unit)

CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/BasisForScatt.txt',Form='FORMATTED',ACCESS='APPEND')
WRITE(Out_Unit,*)'Completed Asymptotic_Analysis'
WRITE(Out_Unit,*)

RETURN
ENDSUBROUTINE Asymptotic_Analysis
