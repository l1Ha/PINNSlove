SUBROUTINE ReadAllData
USE Numeric_Kinds_Module
USE FileUnits_Module
USE InputFile_Module
USE NameList_Module
IMPLICIT NONE
LOGICAL, Parameter :: Prt_Masses=.True.
INTEGER I, IERR, N_Irrep, j, nu, k, omega, Ichanl, IFlag, iarran
INTEGER, EXTERNAL :: CalcParm
CHARACTER(LEN=2) Irreps(6)
CHARACTER(LEN=3) Triatomic, SymGroup
CHARACTER(LEN=256) nl_iomsg

WRITE(Out_Unit,*)'Called    ReadAllData'
CLOSE(Out_Unit)

OPEN(Unit=Out_Unit,File=TRIM(OutDIR)//'Output/ReadAllData.txt',Form='FORMATTED')
CALL NumericInformation

! Input IFileName
OPEN(Unit=In_Unit,File=TRIM(InputDIR)//'Input/AllBasis.txt',form='formatted',status='old',action='read')
   REWIND(In_Unit)
   READ(In_Unit,  NML=IFileName, IOSTAT=IERR)
   IF (IERR/=0)THEN
      WRITE(Msg_Unit,*)' Missing a necessary namelist: IFileName'
      STOP ' Missing a necessary namelist: IFileName'
   ENDIF
   WRITE(Out_Unit,NML=IFilename)
CLOSE(In_Unit)
WRITE(Out_Unit,'("InputFile=",A)')InputFile(1:LEN(InputFile))

! Atomic Masses
INQUIRE (File=InputDIR(1:LEN(TRIM(InputDIR)))//TRIM(InputFile), exist = there)
IF(there)THEN
   Triatomic="ABC"
   SymGroup="MxS"
   N_Irrep=1
   Irreps=["A"," "," "," "," "," "]
   REWIND In_Unit
   OPEN(Unit=In_Unit,File=TRIM(InputDIR)//TRIM(InputFile),status='old',action='read')
      READ(In_Unit,  NML=AtomicMasses, IOSTAT=IERR)                            ! Atomic Masses
      IF(IERR==-1)THEN
         IERR=0
         REWIND In_Unit
         READ(In_Unit, NML=Atoms, IOSTAT=IERR)
         IF(Prt_Masses)THEN
            WRITE(Out_Unit,NML=Atoms)
            WRITE(Msg_Unit,NML=Atoms)
            WRITE(Out_Unit,*)
         ENDIF
         DO I=1,3
            CALL AtomicWeights(AtomicSymbol(I), MassNumber(I), AtomicNumber(I), Mass(I), Abundance(I), AtomicWeight(I), Notes(I))
            IF(Prt_Masses)THEN
            WRITE(Out_Unit, '(A,I3)')'Atomic Number = ', AtomicNumber(I)
            WRITE(Out_Unit, '(A,A3)')'Atomic Symbol = ', TRIM(AtomicSymbol(I))
            WRITE(Out_Unit, '(A,I3)')'Mass Number = '  , MassNumber(I)
            WRITE(Out_Unit, '(A,1PES23.15)')'Relative Atomic Mass = '  , Mass(I)
            WRITE(Out_Unit, '(A,1PES23.15)')'Isotopic Composition = '  , Abundance(I)
            WRITE(Out_Unit, '(A,1PES23.15)')'Standard Atomic Weight = ', AtomicWeight(I)
            WRITE(Out_Unit, '(A,A)')'Notes = ', TRIM(Notes(I))
            WRITE(Out_Unit,*)
            ENDIF
         ENDDO
         IF(AtomicSymbol(2)/=AtomicSymbol(3).or.MassNumber(2)/=Massnumber(3))THEN
            Triatomic="ABC"
            IF(SymGroup=="MxS")THEN
               SymGroup="C2"
               N_Irrep=2
               Irreps=["A ","B ","  ","  ","  ","  "]
            ENDIF
         ELSEIF(AtomicSymbol(1)==AtomicSymbol(2).and.MassNumber(1)==MassNumber(2))THEN
            Triatomic="A3 "            
            IF(SymGroup=="MxS")THEN
               SymGroup="C6v"
               N_Irrep=6
               Irreps=["A1","A2","B1","B2","E1","E2"]
            ENDIF
         ELSE
            Triatomic="AB2"
            IF(SymGroup=="MxS")THEN
               SymGroup="C2v"
               N_Irrep=4
               Irreps=["A1","A2","B1","B2","  ","  "]
            ENDIF
         ENDIF
         WRITE(Out_Unit,'(A10,A3)')"Triatomic=", Triatomic
         WRITE(Out_Unit,'(A10,A3)')"SymGroup= ",SymGroup 
         WRITE(Out_Unit,'(A10,I2)')"N_Irrep=  ",N_Irrep 
         WRITE(Out_Unit,'(A10,6A3)')"Irreps=   ",(Irreps(i),i=1,N_Irrep)
         
         WRITE(Msg_Unit,'(A10,A3)')"Triatomic=", Triatomic
         WRITE(Msg_Unit,'(A10,A3)')"SymGroup= ",SymGroup 
         WRITE(Msg_Unit,'(A10,I2)')"N_Irrep=  ",N_Irrep 
         WRITE(Msg_Unit,'(A10,6A3)')"Irreps=   ",(Irreps(i),i=1,N_Irrep)
      ELSEIF (IERR/=0)THEN
         WRITE(Msg_Unit,*)' Missing a necessary namelist: Atoms or AtomicMasses'
         STOP ' Missing a necessary namelist: Atoms or AtomicMasses'
      ENDIF
      WRITE(Out_Unit,NML=AtomicMasses)
   CLOSE(In_Unit)
ELSE
   WRITE(Out_Unit,  * ) 'Data file must exist: ', InputDIR(1:LEN(TRIM(InputDIR)))//InputFile
   WRITE(Msg_Unit,  * ) 'Data file must exist: ', InputDIR(1:LEN(TRIM(InputDIR)))//InputFile
   STOP "ReadAllData"
ENDIF

! Read in Name for the Potential
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile)
   READ(In_Unit,  NML=PotentialNames, IOSTAT=IERR)                         ! Potential Names
   WRITE(Out_Unit,NML=PotentialNames)
   IF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR no PotentialNames'
      WRITE(Msg_Unit,*)'ERROR no PotentialNames'
      STOP 'No Potential'
   ENDIF
CLOSE(In_Unit)

! Read in FlowOptions
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile)
   READ(In_Unit,  NML=FlowOptions, IOSTAT=IERR)                          ! Flow Options
   WRITE(Out_Unit,NML=FlowOptions)
   IF(IERR/=0)THEN
      WRITE(Out_Unit,*)'Warning using default FlowOptions'
      WRITE(Msg_Unit,*)'Warning using default FlowOptions'
   ENDIF
CLOSE(In_Unit)

! Get Physical Constants
CALL Constant

! Calculate Mass Factors
CALL MassFac(.True.)

! Sector Definitions
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//TRIM(InputFile),form='formatted',status='old')
   INQUIRE(file=InputDIR(1:LEN(TRIM(InputDIR)))//TRIM(InputFile), exist=there)
   READ(In_Unit,  NML=SectorDef, IOSTAT=IERR)
   WRITE(Out_Unit,NML=SectorDef)
   IF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR no Namelist SectorDef'
      WRITE(Msg_Unit,*)'ERROR no Namelist SectorDef'
      STOP 'ERROR no Namelist SectorDef'
   ENDIF
   RhoMin=RhoStart
   RhoMax=RStartJacobi
   Nrho=NSectors
CLOSE(In_Unit)

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
   READ(In_Unit,  NML=sfundist, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'SFunDist')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist SFunDist'
      WRITE(Msg_Unit,*)'ERROR with Namelist SFunDist'
      STOP 'ReadAllData: SFunDist'
   ENDIF
   IFlag=CalcParm(NRho, RhoMin, DeltaRho1, DeltaRho2, RhoMax, Kase_CalcParm)
   IF(IFlag/=0)THEN
      WRITE(Out_Unit,'("Error Boundaries: IFlag=",L2)')IFlag
      WRITE(Msg_Unit,'("Error Boundaries: IFlag=",L2)')IFlag
      STOP 'ReadAllData: IFlag/=0'
   ENDIF
   WRITE(Out_Unit,NML=sfundist)
CLOSE(In_Unit)

! Contour Plotting Data
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
   READ(In_Unit,  NML=contour, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'Contour')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Contour11'
      WRITE(Msg_Unit,*)'ERROR with Namelist Contour12'
   ENDIF
   WRITE(Out_Unit,NML=contour)
CLOSE(In_Unit)

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
READ(In_Unit,  NML=quantum, IOSTAT=IERR)
IF(IERR/=0)THEN
   WRITE(Out_Unit,*)'ERROR with Namelist Quantum'
   WRITE(Msg_Unit,*)'ERROR with Namelist Quantum'
   STOP 'ReadAllData: Quantum'
ENDIF
CLOSE(In_Unit)
DO iarran=1,3
   DO nu=MinVib(iarran),MaxVib(iarran)
      jmax(nu,iarran)=jmax(0,iarran)
   ENDDO
ENDDO
   nfreq=0
   DO Omega=0,JTot
      DO k=1,3
         DO nu=minvib(k),maxvib(k)
            DO j=jmin(nu,k,Omega),jmax(nu,k)
               nfreq(Omega)=nfreq(Omega)+1
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   numinval=MINVAL(MinVib)
   numaxval=MAXVAL(MaxVib)
   jminval=MINVAL(jmin)
   jmaxval=MAXVAL(jmax)
   WRITE(Out_Unit,NML=quantum)

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
   READ(In_Unit,  NML=momentum, IOSTAT=IERR)
   IF(IERR/=0)THEN
      WRITE(Out_Unit,*)'Warning no momentum used'
      WRITE(Msg_Unit,*)'Warning no momentum used'
   ENDIF
   WRITE(Out_Unit,NML=momentum)
CLOSE(In_Unit)

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
   READ(In_Unit, NML=convrg, IOSTAT=IERR)
   nsfunc=SUM(nfreq)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'Convrg')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Convrg'
      WRITE(Msg_Unit,*)'ERROR with Namelist Convrg'
      STOP 'ReadAllData: Convrg'
   ENDIF
   IF(ngood<=1.or.ngood>nfreq(0))ngood=max(1,nfreq(0)/2)
   IF(nave1<=1)nave1=1
   IF(nave2<=1.or.nave2>ngood)nave2=max(1,ngood/2)
   WRITE(Out_Unit,NML=convrg)
CLOSE(In_Unit)

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
   READ(In_Unit,  NML=dipole, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'Dipole')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Dipole'
      WRITE(Msg_Unit,*)'ERROR with Namelist Dipole'
   ENDIF
   WRITE(Out_Unit,NML=dipole)
CLOSE(In_Unit)

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
   READ(In_Unit,  NML=etachanl, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'EtaChanl')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist EtaChanl'
      WRITE(Msg_Unit,*)'ERROR with Namelist EtaChanl'
      STOP 'ReadAllData: EtaChanl'
   ENDIF
   WRITE(Out_Unit,NML=etachanl)
CLOSE(In_Unit)

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
   READ(In_Unit,  NML=fbr, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'FBR')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist FBR'
      WRITE(Msg_Unit,*)'ERROR with Namelist FBR'
      STOP 'ReadAllData: FBR'
   ENDIF
   WRITE(Out_Unit,NML=fbr)
CLOSE(In_Unit)

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
   READ(In_Unit,  NML=gauss, IOSTAT=IERR)
   DO k=1,Narran
      nlegndre(k)=jmax(0,k)
   ENDDO
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'Gauss')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Gauss'
      WRITE(Msg_Unit,*)'ERROR with Namelist Gauss'
      STOP 'ReadAllData: Gauss'
   ENDIF
   WRITE(Out_Unit,NML=gauss)
CLOSE(In_Unit)

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
   READ(In_Unit,  NML=IntChanl, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'IntChanl')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist IntChanl'
      WRITE(Msg_Unit,*)'ERROR with Namelist IntChanl'
      STOP 'ReadAllData: IntChanl'
   ENDIF
   WRITE(Out_Unit,NML=IntChanl)
CLOSE(In_Unit)

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
   READ(In_Unit,  NML=integopt, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'IntegOpt')
      int_opt='ABM_Integrals'
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist IntegOpt'
      WRITE(Msg_Unit,*)'ERROR with Namelist IntegOpt'
      STOP 'ReadAllData: IntegOpt'
   ENDIF
   WRITE(Out_Unit,NML=integopt)
CLOSE(In_Unit)

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
   READ(In_Unit,  NML=logder, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'Logder')
      EDeriv=.True.
      nsteps=20
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Logder'
      WRITE(Msg_Unit,*)'ERROR with Namelist Logder'
      STOP 'ReadAllData: Logder'
   ENDIF
   WRITE(Out_Unit,NML=logder)
CLOSE(In_Unit)

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
   READ(In_Unit,  NML=mshdat, IOSTAT=IERR)
   nsfunc=SUM(NFreq)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'MshDat')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist MshDat'
      WRITE(Msg_Unit,*)'ERROR with Namelist MshDat'
      STOP 'ReadAllData: MshDat'
   ENDIF
   WRITE(Out_Unit,NML=mshdat)
CLOSE(In_Unit)

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
   READ(In_Unit,  NML=options, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'Options')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Options'
      WRITE(Msg_Unit,*)'ERROR with Namelist Options'
      STOP 'ReadAllData: Options'
   ENDIF
   WRITE(Out_Unit,NML=options)
CLOSE(In_Unit)

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
   !READ(In_Unit,  NML=plot, IOSTAT=IERR)
   !WRITE(Out_Unit,NML=Plot)
CLOSE(In_Unit)

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
   READ(In_Unit,  NML=potzero, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'PotZero')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist PotZero'
      WRITE(Msg_Unit,*)'ERROR with Namelist PotZero'
      STOP 'ReadAllData: PotZero'
   ENDIF
   WRITE(Out_Unit,NML=PotZero)
CLOSE(In_Unit)

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
   READ(In_Unit,  NML=qcasep, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'QCasep')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist QCasep'
      WRITE(Msg_Unit,*)'ERROR with Namelist QCasep'
      STOP 'ReadAllData: QCasep'
   ENDIF
   WRITE(Out_Unit,NML=qcasep)
CLOSE(In_Unit)

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
   READ(In_Unit,  NML=qdivp, IOSTAT=IERR, IOMSG=nl_iomsg)
   IF(IERR/=0 .AND. IERR/=-1) WRITE(Msg_Unit,*)'DEBUG QDivp IERR=',IERR,' MSG=',TRIM(nl_iomsg)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'QDivp')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist QDivp'
      WRITE(Msg_Unit,*)'ERROR with Namelist QDivp'
      STOP 'ReadAllData: QDivp'
   ENDIF
   WRITE(Out_Unit,NML=qdivp)
CLOSE(In_Unit)

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
   READ(In_Unit,  NML=regins, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'Regins')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Regins'
      WRITE(Msg_Unit,*)'ERROR with Namelist Regins'
      STOP 'ReadAllData: Regins'
   ENDIF
   WRITE(Out_Unit,NML=regins)
CLOSE(In_Unit)

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
   READ(In_Unit,  NML=rstart, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'RStart')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist RStart'
      WRITE(Msg_Unit,*)'ERROR with Namelist RStart'
      STOP 'ReadAllData: RStart'
   ENDIF
   WRITE(Out_Unit,NML=rstart)
CLOSE(In_Unit)

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
   READ(In_Unit,  NML=savnam, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'SavNam')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist SavNam'
      WRITE(Msg_Unit,*)'ERROR with Namelist SavNam'
      STOP 'ReadAllData: SavNam'
   ENDIF
   WRITE(Out_Unit,NML=savnam)
CLOSE(In_Unit)

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
   READ(In_Unit,  NML=spectro, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'Spectro')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Spectro'
      WRITE(Msg_Unit,*)'ERROR with Namelist Spectro'
      STOP 'ReadAllData: Spectro'
   ENDIF
   WRITE(Out_Unit,NML=spectro)
CLOSE(In_Unit)

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
   READ(In_Unit,  NML=storage, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'Storage')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Storage'
      WRITE(Msg_Unit,*)'ERROR with Namelist Storage'
      STOP 'ReadAllData: Storage'
   ENDIF
   WRITE(Out_Unit,NML=storage)
CLOSE(In_Unit)

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
   READ(In_Unit,  NML=stst, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'STST')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist STST'
      WRITE(Msg_Unit,*)'ERROR with Namelist STST'
      STOP 'ReadAllData: STST'
   ENDIF
   WRITE(Out_Unit,NML=stst)
CLOSE(In_Unit)

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
   READ(In_Unit,  NML=system, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'System')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist System'
      WRITE(Msg_Unit,*)'ERROR with Namelist System'
      STOP 'ReadAllData: System'
   ENDIF
   WRITE(Out_Unit,NML=system)
CLOSE(In_Unit)

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
   READ(In_Unit,  NML=sysu, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'sysu')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist sysu'
      WRITE(Msg_Unit,*)'ERROR with Namelist sysu'
      STOP 'ReadAllData: sysu'
   ENDIF
   WRITE(Out_Unit,NML=sysu)
CLOSE(In_Unit)

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
   READ(In_Unit,  NML=testing, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'Testing')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Testing'
      WRITE(Msg_Unit,*)'ERROR with Namelist Testing'
      STOP 'ReadAllData: Testing'
   ENDIF
   WRITE(Out_Unit,NML=testing)
CLOSE(In_Unit)

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
   READ(In_Unit,  NML=Title_Labels, IOSTAT=IERR)
   IF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Title_Labels'
      WRITE(Msg_Unit,*)'ERROR with Namelist Title_Labels'
      STOP 'ReadAllData: Title_Labels'
   ENDIF
   WRITE(Out_Unit,NML=Title_Labels)
CLOSE(In_Unit)

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
   READ(In_Unit,  NML=TotEnergy, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'TotEnergy')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist TotEnergy'
      WRITE(Msg_Unit,*)'ERROR with Namelist TotEnergy'
      STOP 'ReadAllData: TotEnergy'
   ENDIF
   IF(eV_Input)THEN
      Efirst=Efirst_eV/autoeV
      DeltaEng=DeltaEng_eV/autoeV
   ENDIF
   WRITE(Out_Unit,NML=TotEnergy)
CLOSE(In_Unit)

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,form='formatted',status='old')
   READ(In_Unit,  NML=vivs, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'Vivs')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Vivs'
      WRITE(Msg_Unit,*)'ERROR with Namelist Vivs'
      STOP 'ReadAllData: Vivs'
   ENDIF
   WRITE(Out_Unit,NML=VIVS)
CLOSE(In_Unit)

!--------------------------
!READ(In_Unit,NML=basis_n, IOSTAT=IERR)
!READ(In_Unit,NML=dafthetaopt, IOSTAT=IERR)
!READ(In_Unit,NML=debug, IOSTAT=IERR)
!READ(In_Unit,NML=dipole2, IOSTAT=IERR)
!READ(In_Unit,NML=gauss, IOSTAT=IERR)
!READ(In_Unit,NML=kmatrix, IOSTAT=IERR)
!READ(In_Unit,NML=pmaovr, IOSTAT=IERR)
!READ(In_Unit,NML=potpar, IOSTAT=IERR)
!READ(In_Unit,NML=sfplot_ABM, IOSTAT=IERR)
!READ(In_Unit,NML=umat_trans, IOSTAT=IERR)
!------------------------------

CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/QuantumStates.txt',Form='FORMATTED')
DO ichanl=1,Narran
  WRITE(*,*)
  WRITE(*,*)"Ichanl=", IChanl, " Jtot=", Jtot, "  jmax(Ichanl)=", jmax(0,Ichanl)
  WRITE(Out_Unit,*)
  WRITE(Out_Unit,*)"Ichanl=", IChanl, "  jmax(Ichanl)=", jmax(0,Ichanl)
  CALL QuantumStates(Jtot,jmax(0,ichanl))
ENDDO

CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/BasisForScatt.txt',Form='FORMATTED',ACCESS='APPEND')
WRITE(Out_Unit,*)'Completed ReadAllData'
WRITE(Out_Unit,*)

RETURN
ENDSUBROUTINE ReadAllData
