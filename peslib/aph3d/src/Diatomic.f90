SUBROUTINE Diatomic
USE DiatomicPot_Module
USE FileUnits_Module
USE FileUnits_OneDim_Module
USE OneDim_Parms_Module
USE PES_Module
USE NArran_Module
USE Convrsns_Module
USE Masses_Module
USE Quantb_module, ONLY: JminFull=>Jmin, JmaxFull=>Jmax, MinVibFull=>MinVib, MaxVibFull=>MaxVib
IMPLICIT NONE
INTEGER Ip1, Ip2
INTEGER :: DtValues(8)
CHARACTER(LEN=10) hour
CHARACTER(LEN=8) today
CHARACTER(LEN=5) curzone
CHARACTER(LEN=3) AS, BS

WRITE(Out_Unit,*)
WRITE(Out_unit,*)"At the beginning of Diatomic"
CALL Date_And_Time(today, hour, curzone, dtvalues)
CALL Print_Date_And_Time(Today, Hour, CurZone, DtValues, Out_Unit)

!PES_Name='H3 Porter-Karplus PES    '
FitDAFEnergies=.True.         ! IF TRUE the DAF eigenenergies will be fit
FitNumerovEnergies=.True.     ! IF TRUE the Numerov eigenenergies will be fit
FitFourierGridEnergies=.True. ! IF TRUE the FourierGrid eigenenergies will be fit
Continuum=.False.             ! IF TRUE Scattering eigenfunctions will be calculated
jmin=0                        ! Minimum rotational state (usually 0)
jmax=6                        ! Maximum rotational state
MinElec=0                     ! Minimum electronic state (usually 0)
MaxElec=0                     ! Maximum electronic state (usually 0 or 1)
NEnergy=101                   ! Number of continuum energies to calculate
Emin=1.d-15                   ! Minimum continuum energy (Hartree Units)
Emax=1.d-05                   ! Maximum continuum energy (Hartree Units)

CALL Constant
CALL ReadAllData

OPEN(Unit=Diatomic_Out_Unit,File=TRIM(Outdir)//"DiatomicOut/DiatomicOut.txt")
DO IARRAN=1,3
   Ith_Arrangement=IARRAN
   IF(IArran==1)THEN
       Ip1=2
       Ip2=3
   ELSEIF(IArran==2)THEN
       Ip1=3
       Ip2=1
   ELSE
       Ip1=1
       Ip2=2
   ENDIF
   jmin=JminFull(0,Iarran,0)
   jmax=JmaxFull(0,Iarran)
   MassA=mass(ip1)/amutoau
   MassB=mass(ip2)/amutoau
   ReducedMass=amutoau*MassA*MassB/(MassA+MassB)
   DiatomicPot=TRIM(AtomicSymbol(ip1))//TRIM(AtomicSymbol(ip2))
   WRITE(Msg_Unit,*)
   WRITE(Msg_Unit,*)DiatomicPot
   WRITE(Msg_Unit,*)"ip1=", ip1,"  ip2=",ip2
   WRITE(Msg_Unit,*)'Diatom=', TRIM(AtomicSymbol(ip1)),TRIM(AtomicSymbol(ip2))
   WRITE(Msg_Unit,*)'Diatomic Masses=', MassA, MassB
   WRITE(Msg_Unit,*)"Calling OneDimTISE"
   
   
   WRITE(Diatomic_Out_Unit,*)
   WRITE(Diatomic_Out_Unit,*)DiatomicPot
   WRITE(Diatomic_Out_Unit,*)"ip1=", ip1,"  ip2=",ip2
   WRITE(Diatomic_Out_Unit,*)'Diatom=', TRIM(AtomicSymbol(ip1)),TRIM(AtomicSymbol(ip2))
   WRITE(Diatomic_Out_Unit,*)'Diatomic Masses=', MassA, MassB
   WRITE(Diatomic_Out_Unit,*)"Calling OneDimTISE"
   CALL OneDimTISE(MaxVibFull(IARRAN))
ENDDO

! Compute the Dunham parameters for each arrangement channels
OPEN(Unit=YFitBin_Unit,file=TRIM(OutDIR)//"BinOut/"//"Yfit.Bin", form='unformatted',status='old',action='read')
OPEN(Unit=Yfitter_Unit,file=TRIM(OutDIR)//"DiatomicOut/"//"Yfit.txt", form='formatted',status='unknown',action='write') 
DO IARRAN=1,3
   Ith_Arrangement=IARRAN
   IF(IArran==1)THEN
       Ip1=2
       Ip2=3
   ELSEIF(IArran==2)THEN
       Ip1=3
       Ip2=1
   ELSE
       Ip1=1
       Ip2=2
   ENDIF
   jmin=JminFull(0,Iarran,0)
   jmax=JmaxFull(0,Iarran)
   MassA=mass(ip1)/amutoau
   MassB=mass(ip2)/amutoau
   ReducedMass=amutoau*MassA*MassB/(MassA+MassB)
   DiatomicPot=TRIM(AtomicSymbol(ip1))//TRIM(AtomicSymbol(ip2))
   CALL YFitter
ENDDO
CLOSE(YFitBin_Unit)
CLOSE(Yfitter_Unit)

CLOSE(UNIT=Diatomic_Out_Unit)

WRITE(Out_unit,*)"At the end of Diatomic"
CALL Date_And_Time(today, hour, curzone, dtvalues)
CALL Print_Date_And_Time(Today, Hour, CurZone, DtValues, Out_Unit)
WRITE(Msg_Unit,*)"End of Execution"
WRITE(Msg_Unit,*)"Completed Diatomic"
WRITE(Out_Unit,*)"End of Execution"
WRITE(Out_Unit,*)"Completed Diatomic"

RETURN
ENDSUBROUTINE Diatomic