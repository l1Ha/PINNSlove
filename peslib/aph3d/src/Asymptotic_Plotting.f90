SUBROUTINE Asymptotic_Plotting(AFile,BFile,CFile, RowFrame, ColumnFrame)
! Author: Gregory A. Parker
! Department of Physics and Astronomy
! University of Oklahoma
USE Numeric_Kinds_Module
USE Narran_Module
USE FileUnits_Asymptotic_Module
USE Energy_Module
USE EDeriv_Module
USE ReadK_Module
USE InputFile_Module
USE Convrsns_Module
USE QAsy_Numbers_Module, ONLY : NQuant, QNumbr
IMPLICIT NONE
INTEGER, PARAMETER:: NPRT=300
INTEGER M, k, DtValues(8), Initial_Channel, Final_Channel, ElasticStates, ElasticMax, jtime/0/
INTEGER NOPEN, NOpen_Max, NopenMax, Num_Elastic, Num_React, Num_Inelastic, index, LastQn
INTEGER, ALLOCATABLE :: QI_Elastic(:,:), QF_Elastic(:,:)
INTEGER, ALLOCATABLE :: QI_React(:,:), QF_React(:,:)
INTEGER, ALLOCATABLE :: QI_Inelastic(:,:), QF_Inelastic(:,:)  
LOGICAL There, PrntQuant, Elastic, q
REAL(Kind=WP_Kind) EnergyNow
REAL(Kind=WP_Kind),ALLOCATABLE :: AMatx(:,:), ALIST(:)
REAL(Kind=WP_Kind),ALLOCATABLE :: BMatx(:,:), BLIST(:)
INTEGER i, f, TypePlot, isave
CHARACTER(LEN=5)  CurZone
CHARACTER(LEN=6) RowFrame, ColumnFrame
CHARACTER(LEN=8)  Today
CHARACTER(LEN=10) Hour
CHARACTER(LEN=*) AFile
CHARACTER(LEN=*) BFile
CHARACTER(LEN=*) CFile
CHARACTER(LEN=19), PARAMETER:: ProcName='Asymptotic_Plotting'
CHARACTER(LEN=6) Print_Flag
CHARACTER(LEN=2)  NQuantStr
CHARACTER(LEN=5)  NQuantStr5
CHARACTER(LEN=30) NQuantFormatStr
CHARACTER(LEN=32) NQuantFormatStr32
CHARACTER(LEN=63) NQuantFormatStr60
CHARACTER(LEN=63) NQuantFormatStr62, QType, QTyped
INTEGER kenergy

CALL PoptAsy(ProcName, Print_Flag)
IF(AFile==BFile)THEN
   TypePlot=1
ELSE
   TypePlot=2
ENDIF

WRITE(AP_Unit,*)
WRITE(AP_Unit,*)"Currently in: ",TRIM(ProcName)

IF(AFile=='Phase'.or.AFile=='PhaseE')THEN
   TypePlot=0
ENDIF

IF(CFile=='Probability')THEN
   PrntQuant=.True.
   OPEN(Unit=Prob_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/Prob.txt',Form='formatted', Status='unknown')
ELSE
   PrntQuant=.False.
ENDIF

CALL Date_And_Time(today, hour, curzone, dtvalues)
CALL Print_Date_And_Time(Today, Hour, CurZone, DtValues, Out_Unit)
WRITE(AP_Unit,*)
WRITE(AP_Unit,*)
WRITE(AP_Unit,*)
WRITE(AP_Unit,*)'Entering: ', ProcName


INQUIRE(file=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile, exist=there)
IF(THERE)THEN
   OPEN(In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile)
   READ(In_Unit,NML=TotEnergy)
   IF(eV_Input)THEN
      Efirst=Efirst_eV/autoeV
      DeltaEng=DeltaEng_eV/autoeV
   ENDIF
   WRITE(AP_Unit,NML=TotEnergy)
CLOSE(In_Unit)
ENDIF

NopenMax=NPRT
ALLOCATE(ALIST(NopenMax*NopenMax), BLIST(NopenMax*NopenMax))
ALLOCATE(QI_Elastic(NQuant,NopenMax*NopenMax), QF_Elastic(NQuant,NopenMax*NopenMax) )
ALLOCATE(QI_React(NQuant,NopenMax*NopenMax), QF_React(NQuant,NopenMax*NopenMax) )
ALLOCATE(QI_Inelastic(NQuant,NopenMax*NopenMax), QF_Inelastic(NQuant,NopenMax*NopenMax) )

DO Initial_Channel=1,Narran
   DO Final_Channel=1,Narran
      ElasticMax=1
      IF(Initial_Channel==Final_Channel)ElasticMax=2
      DO ElasticStates=1,ElasticMax
         IF(Initial_Channel==Final_Channel)THEN
            IF(ElasticStates==1)THEN
               QType="_Elastic"
            ELSE
               QType="_Inelastic"
            ENDIF
         ELSE
            QType="_Reactive"
         ENDIF
         WRITE(NQuantStr5,'(A1,I1,A2,I1)') "_",Initial_Channel,"to",Final_Channel
         QTyped=TRIM(QType)//NQuantStr5
         WRITE(AP_Unit,'(A)') TRIM(QTyped)
         
   OPEN(Unit=In_UnitA,File=TRIM(OutDIR)//'BinOut/'//AFile//'.bin',Form='unformatted', Status='old')
   IF(AFile/=BFile)THEN
      OPEN(Unit=In_UnitB,File=TRIM(OutDIR)//'BinOut/'//BFile//'.bin',Form='unformatted', Status='old')
   ENDIF
   IF(TRIM(RowFrame)=='Delves'.and.TRIM(ColumnFrame)=='Delves')THEN
      OPEN(Unit=Plt_Unit,File=TRIM(OutDIR)//'GraphicsOut/'//CFile//TRIM(QTyped)//'.csv',Form='formatted', Status='unknown')
   ELSE
      OPEN(Unit=Plt_Unit,File=TRIM(OutDIR)//'GraphicsOut/'//CFile//'.csv',Form='formatted', Status='unknown')
   ENDIF
   WRITE(AP_Unit,'(6(A,2x))')'AFile=', AFile,'BFile=', BFile, 'CFile=', CFile
   
   ! Create and Print headings
   DO kenergy=1,nenergy
      energynow=efirst+deltaeng*(kenergy-1)
      CALL Read_K1(NOPEN,In_UnitA,AP_Unit)
      flush(msg_unit)
      IF(AFile/=BFile)CALL Read_K1(NOPEN,In_UnitB,AP_Unit)
      ALLOCATE(AMatx(nopen,nopen), BMatx(NOpen,NOpen))
      !IF(kenergy==nenergy)ALLOCATE(ALIST(nopen*nopen), BLIST(NOpen*NOpen))
      CALL Read_K2(NOPEN,AMatx,In_UnitA,AP_Unit)
      IF(AFile/=BFile)CALL Read_K2(NOPEN,BMatx,In_UnitB,AP_Unit)
      CALL SortEnergy(AsymptEnergies,LastLoc,NOpen)
      M=MIN(NOPEN,NPRT)
      Nopen_Max=MIN(Nopen,NPRT)
      WRITE(NQuantStr,'(I1)') NQuant
      WRITE(NQuantFormatStr62,*)'("Final  ",",","Etot(au)",",","Etot(eV)",",",190001('//NQuantStr//'i3,","))'
      WRITE(NQuantFormatStr60,*)'("Initial",",","Etot(au)",",","Etot(eV)",",",190001('//NQuantStr//'i3,","))'
      IF(kenergy==nenergy)THEN
         !ALLOCATE(QI_Elastic(NQuant,NOpen*NOpen), QF_Elastic(NQuant,NOpen*NOpen) )
         !ALLOCATE(QI_React(NQuant,NOpen*NOpen), QF_React(NQuant,NOpen*NOpen) )
         !ALLOCATE(QI_Inelastic(NQuant,NOpen*NOpen), QF_Inelastic(NQuant,NOpen*NOpen) )
         IF(TRIM(RowFrame)=='Delves'.and.TRIM(ColumnFrame)=='Delves')THEN
            QI_Elastic=0
            QF_Elastic=0
            QI_React=0
            QF_React=0
            QI_Inelastic=0
            QF_Inelastic=0
            Num_Elastic=0
            Num_React=0
            Num_Inelastic=0
            DO f=1,NOpen_Max
               DO i=1,NOpen_Max
                  Elastic=.True.
                  IF(TRIM(QType)/="_Elastic")Elastic=.False.
                  IF(i/=f)Elastic=.False.
                  IF(Initial_Channel/=QNumbr(1,LastLoc(i)))Elastic=.false.
                  IF(Elastic)THEN
                     Num_Elastic=Num_Elastic+1
                     ElasLoc(Num_Elastic)=i+(f-1)*NOpen_Max
                     QI_Elastic(:,Num_Elastic)=QNumbr(:,LastLoc(i))
                     QF_Elastic(:,Num_Elastic)=QNumbr(:,LastLoc(f))
                  ELSEIF(QNumbr(1,LastLoc(i))/=QNumbr(1,LastLoc(f)).and.QNumbr(1,LastLoc(f))==Initial_Channel.and.QNumbr(1,LastLoc(i))==Final_Channel)THEN !TMPMOD
                     Num_React=Num_React+1
                     ReactLoc(Num_React)=i+(f-1)*NOpen_Max
                     QI_React(:,Num_React)=QNumbr(:,LastLoc(i))
                     QF_React(:,Num_React)=QNumbr(:,LastLoc(f))
                  ELSE
                     Elastic=.True.
                     DO k=1,NQuant
                        IF(QNumbr(k,LastLoc(i))/=QNumbr(k,LastLoc(f)))Elastic=.False.
                     ENDDO
                     IF(.not.Elastic.and.QNumbr(1,LastLoc(i))==Initial_Channel.and.QNumbr(1,LastLoc(f))==Initial_Channel)THEN
                        Num_Inelastic=Num_Inelastic+1
                        InelasLoc(Num_InElastic)=i+(f-1)*NOpen_Max
                        QI_Inelastic(:,Num_Inelastic)=QNumbr(:,LastLoc(i))
                        QF_Inelastic(:,Num_Inelastic)=QNumbr(:,LastLoc(f))
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
            IF(TRIM(QType)=="_Elastic")THEN
               IF(TypePlot==1)THEN
                  WRITE(Plt_Unit,NQuantFormatStr62)((QF_Elastic(k,i),k=1,NQuant),i=1,Num_Elastic)
                  WRITE(Plt_Unit,NQuantFormatStr60)((QI_Elastic(k,f),k=1,NQuant),f=1,Num_Elastic)
               ELSEIF(TypePlot==2)THEN
                  WRITE(Plt_Unit,NQuantFormatStr62)((QF_Elastic(k,i),k=1,NQuant),(QI_Elastic(k,i),k=1,NQuant),i=1,Num_Elastic)
                  WRITE(Plt_Unit,NQuantFormatStr60)((QI_Elastic(k,f),k=1,NQuant),(QF_Elastic(k,f),k=1,NQuant),f=1,Num_Elastic)
               ENDIF
            ELSEIF(TRIM(QType)=="_Reactive")THEN
               IF(TypePlot==1)THEN
                  WRITE(Plt_Unit,NQuantFormatStr62)((QF_React(k,i),k=1,NQuant),i=1,Num_React)
                  WRITE(Plt_Unit,NQuantFormatStr60)((QI_React(k,f),k=1,NQuant),f=1,Num_React)
               ELSEIF(TypePlot==2)THEN
                  WRITE(Plt_Unit,NQuantFormatStr62)((QF_React(k,i),k=1,NQuant),(QI_React(k,i),k=1,NQuant),i=1,Num_React)
                  WRITE(Plt_Unit,NQuantFormatStr60)((QI_React(k,f),k=1,NQuant),(QF_React(k,f),k=1,NQuant),f=1,Num_React)
               ENDIF
            ELSE
               IF(TypePlot==1)THEN
                  WRITE(Plt_Unit,NQuantFormatStr62)((QF_Inelastic(k,i),k=1,NQuant),i=1,Num_Inelastic)
                  WRITE(Plt_Unit,NQuantFormatStr60)((QI_Inelastic(k,f),k=1,NQuant),f=1,Num_Inelastic)
               ELSEIF(TypePlot==2)THEN
                  WRITE(Plt_Unit,NQuantFormatStr62)((QF_Inelastic(k,i),k=1,NQuant),(QI_Inelastic(k,i),k=1,NQuant),i=1,Num_Inelastic)
                  WRITE(Plt_Unit,NQuantFormatStr60)((QI_Inelastic(k,f),k=1,NQuant),(QF_Inelastic(k,f),k=1,NQuant),f=1,Num_Inelastic)
               ENDIF
            ENDIF
         ELSEIF(TRIM(RowFrame)=='APH'.and.TRIM(ColumnFrame)=='APH')THEN
            WRITE(Plt_Unit,'("Etot",",",190001(i2,A2,5i2,","))')((i,i=1,NOpen_Max),f=1,NOpen_Max)
            WRITE(Plt_Unit,'("Etot",",",190001(i2,A2,5i2,","))')((f,i=1,NOpen_Max),f=1,NOpen_Max)
         ELSEIF(TRIM(RowFrame)=='Index'.and.TRIM(ColumnFrame)=='Delves')THEN
	         WRITE(NQuantStr,'(I1)') NQuant
            WRITE(Plt_Unit,'(3(","),190001(i5,","))') ((i,i=1,NOpen_Max),f=1,NOpen_Max)
            WRITE(Plt_Unit,NQuantFormatStr62) (((QNumbr(k,f), k=1,NQuant),f=1,NOpen_Max),i=1,NOpen_Max)
         ENDIF
      ENDIF
      DEALLOCATE(AMatx,BMatx)
   ENDDO

   REWIND In_UnitA
   REWIND In_UnitB
   ! Print files for plotting
   DO kenergy=1,nenergy
      energynow=efirst+deltaeng*(kenergy-1)
      CALL Read_K1(NOPEN,In_UnitA,AP_Unit)
      IF(AFile/=BFile)CALL Read_K1(NOPEN,In_UnitB,AP_Unit)
      flush(msg_unit)
      ALLOCATE(AMatx(nopen,nopen), BMatx(NOpen,NOpen))
      amatx=0.d0
      bmatx=0.d0
      CALL Read_K2(NOPEN,AMatx,In_UnitA,AP_Unit)
      
      IF(AFile/=BFile)CALL Read_K2(NOPEN,BMatx,In_UnitB,AP_Unit)
      CALL SortEnergy(AsymptEnergies,EminLoc,NOpen)
      
      M=MIN(NOPEN,NPRT)
      index=0
      IF(Initial_Channel==1.and.Final_Channel==1.and.CFile=='Probability'.and.ElasticStates==1)THEN
      IF(Kenergy==1)THEN
         open(unit=681,file=TRIM(OutDir)//"Output/quant.csv")
         open(unit=682,file=TRIM(OutDir)//"Output/proba.csv")
      ENDIF
!----
!----
!----
      WRITE(681,'(1x,",",190001(I3,",",5I3,","))')(EminLoc(f),(QNumbr(k,(f)),k=1,NQuant),f=1,NOpen)
      
      WRITE(682,*)
      WRITE(682,*)"Kenergy=",Kenergy
      WRITE(682,*)"kenergy=",kenergy      
      WRITE(682,'(1x,",",",",190001(I3,","))') (EminLoc(f),f=1,NOpen)
      WRITE(682,'(1x,",",",",190001(5I3,","))') ((QNumbr(k,EminLoc(f)),k=1,NQuant),f=1,NOpen)
      DO f=1,NOpen
         WRITE(682,'(I3,",",5I3,",",190001(es19.12,","))')EminLoc(f),(QNumbr(k,EminLoc(f)),k=1,NQuant),(AMatx(EminLoc(i),EminLoc(f)),i=1,NOpen) !Min(85,Nopen))
      ENDDO
      ENDIF
      IF(TRIM(RowFrame)=='Delves'.and.TRIM(ColumnFrame)=='Delves')THEN
         IF(TypePlot==0)THEN
            WRITE(Plt_Unit,'(1x,I6,",",190001(es19.12,","))')kenergy,EnergyNow,EnergyNow*autoev,(AMatx(i,1),i=1,M)
         ELSEIF(TypePlot==1)THEN
            index=0
            AList=0.d0
            BList=0.d0
            DO i=1,M
               DO f=1,M
                  Elastic=.True.
                  IF(TRIM(QType)/="_Elastic")Elastic=.False.
                  IF(i/=f)Elastic=.False.
                  IF(Initial_Channel/=QNumbr(1,EminLoc(i)))Elastic=.false.
                  IF(Elastic.and.TRIM(QType)=="_Elastic")THEN
                     IF(index<Num_Elastic)THEN
                        index=index+1
                        AList(index)=AMatx(EminLoc(i),EminLoc(i))
                        IF(TypePlot==2)BList(index)=BMatx(EminLoc(i),EminLoc(i))
                     ELSE
                        WRITE(*,'(i5,4(A,2x))')kenergy,TRIM(QType), TRIM(AFile), TRIM(BFile), TRIM(CFile)
                        WRITE(*,*)M, index, Num_Elastic, i, EminLoc(i), f, EminLoc(f)
                        WRITE(*,*)QNumbr(:,EminLoc(i))
                        WRITE(*,*)QNumbr(:,EminLoc(f))
                        STOP "In Asymptotic_Plotting: Elastic"
                     ENDIF
                  ELSEIF(TRIM(QType)=="_Reactive".and.Initial_Channel==QNumbr(1,EminLoc(i)).and.Final_Channel==QNumbr(1,EminLoc(f)))THEN
                     IF(index<Num_React)THEN
                        IF(Index==0)isave=EminLoc(i)
                        IF(Isave==EminLoc(i))THEN
                        index=index+1
                        AList(index)=AMatx(EminLoc(i),EminLoc(f))
                        IF(TypePlot==2)BList(index)=BMatx(EminLoc(i),EminLoc(f))
                        ENDIF
                     ELSE
                        WRITE(*,'(i5,4(A,2x))')kenergy, TRIM(QType), TRIM(AFile), TRIM(BFile), TRIM(CFile)
                        WRITE(*,*)M, index, Num_React, i, EminLoc(i), f, EminLoc(f)
                        WRITE(*,*)QNumbr(:,EminLoc(i))
                        WRITE(*,*)QNumbr(:,EminLoc(f))
                        STOP "In Asymptotic_Plotting: Reactive"
                     ENDIF
                  ELSEIF(QNumbr(1,EminLoc(i))==QNumbr(1,EminLoc(f)).and.TRIM(QType)=="_Inelastic".and.QNumbr(1,EminLoc(i))==Initial_Channel)THEN
                     Elastic=.True.
                     DO k=1,NQuant
                        IF(QNumbr(k,EminLoc(i))/=QNumbr(k,EminLoc(f)))Elastic=.False.
                     ENDDO
                     IF(.not.Elastic)THEN
                        IF(index<Num_Inelastic)THEN
                        IF(Index==0)isave=EminLoc(i)
                        IF(Isave==EminLoc(i))THEN
                           index=index+1
                           AList(index)=AMatx(EminLoc(i),EminLoc(f))
                           IF(TypePlot==2)BList(index)=BMatx(EminLoc(i),EminLoc(f))
                           ENDIF
                        ELSE
                           WRITE(*,'(i5,4(A,2x))')kenergy, TRIM(QType), TRIM(AFile), TRIM(BFile), TRIM(CFile)
                           WRITE(*,*)M, index, Num_Inelastic, i, EminLoc(i), f, EminLoc(f)
                           WRITE(*,*)QNumbr(:,EminLoc(i))
                           WRITE(*,*)QNumbr(:,EminLoc(f))
                           STOP "In Asymptotic_Plotting: Inelastic"
                        ENDIF
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
            IF(index>0)WRITE(Plt_Unit,'(1x,I6,",",190001(es19.12,","))')kenergy,EnergyNow,EnergyNow*autoev,(AList(i),i=1,index)
         ELSEIF(TypePlot==2)THEN
            index=0
            AList=0.d0
            BList=0.d0
            DO f=1,M
               DO i=1,M
                  Elastic=.True.
                  IF(TRIM(QType)/="_Elastic")Elastic=.False.
                  IF(i/=f)Elastic=.False.
                  IF(Initial_Channel/=QNumbr(1,EminLoc(i)))Elastic=.false.
                  IF(Elastic.and.TRIM(QType)=="_Elastic")THEN
                     IF(index<Num_Elastic)THEN
                        IF(Index==0)isave=EminLoc(i)
                        IF(Isave==EminLoc(i))THEN
                        index=index+1
                        AList(index)=AMatx(EminLoc(i),EminLoc(f))
                        IF(TypePlot==2)BList(index)=BMatx(EminLoc(i),EminLoc(f))
                        ENDIF
                     ELSE
                        WRITE(*,'(i5,4(A,2x))')kenergy,TRIM(QType), TRIM(AFile), TRIM(BFile), TRIM(CFile)
                        WRITE(*,*)M, index, Num_Elastic, i, EminLoc(i), f, EminLoc(f)
                        WRITE(*,*)QNumbr(:,EminLoc(i))
                        WRITE(*,*)QNumbr(:,EminLoc(f))
                        STOP "In Asymptotic_Plotting: Elastic loc2"
                     ENDIF
                  !ELSEIF(QNumbr(1,EminLoc(i))/=QNumbr(1,EminLoc(f)).and.TRIM(QType)=="_Reactive")THEN
                  ELSEIF(TRIM(QType)=="_Reactive".and.Initial_Channel==QNumbr(1,EminLoc(i)).and.Final_Channel==QNumbr(1,EminLoc(f)))THEN
                     IF(index<Num_React)THEN
                        IF(Index==0)isave=EminLoc(i)
                        IF(Isave==EminLoc(i))THEN
                        index=index+1
                        AList(index)=AMatx(EminLoc(i),EminLoc(f))
                        IF(TypePlot==2)BList(index)=BMatx(EminLoc(i),EminLoc(f))
                        ENDIF
                     ELSE
                        WRITE(*,'(i5,4(A,2x))')kenergy, TRIM(QType), TRIM(AFile), TRIM(BFile), TRIM(CFile)
                        WRITE(*,*)M, index, Num_React, i, EminLoc(i), f, EminLoc(f)
                        WRITE(*,*)QNumbr(:,EminLoc(i))
                        WRITE(*,*)QNumbr(:,EminLoc(f))
                        STOP "In Asymptotic_Plotting: Reactive loc2"
                     ENDIF
                  ELSEIF(QNumbr(1,EminLoc(i))==QNumbr(1,EminLoc(f)).and.TRIM(QType)=="_Inelastic".and.QNumbr(1,EminLoc(i))==Initial_Channel)THEN
                     Elastic=.True.
                     DO k=1,NQuant
                        IF(QNumbr(k,EminLoc(i))/=QNumbr(k,EminLoc(f)))Elastic=.False.
                     ENDDO
                     IF(.not.Elastic)THEN
                        IF(index<Num_Inelastic)THEN
                        IF(Index==0)isave=EminLoc(i)
                        IF(Isave==EminLoc(i))THEN
                           index=index+1
                           AList(index)=AMatx(EminLoc(i),EminLoc(f))
                           IF(TypePlot==2)BList(index)=BMatx(EminLoc(i),EminLoc(f))
                           ENDIF
                        ELSE
                           WRITE(*,'(i5,4(A,2x))')kenergy, TRIM(QType), TRIM(AFile), TRIM(BFile), TRIM(CFile)
                           WRITE(*,*)M, index, Num_Inelastic, i, EminLoc(i), f, EminLoc(f)
                           WRITE(*,*)QNumbr(:,EminLoc(i))
                           WRITE(*,*)QNumbr(:,EminLoc(f))
                           STOP "In Asymptotic_Plotting: Inelastic loc2"
                        ENDIF
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
            IF(index>0)WRITE(Plt_Unit,'(1x,I6,",",190001(es19.12,","))')kenergy,EnergyNow,EnergyNow*autoev,(AList(i),BList(i),i=1,index)
         ENDIF
      ELSE
         WRITE(Plt_Unit,'(1x,I6,",",190001(es19.12,","))')kenergy,EnergyNow,EnergyNow*autoev,((AMatx(EminLoc(i),EminLoc(f)),i=1,NOpen),f=1,Nopen)
      ENDIF
      IF(PrntQuant.and.((kenergy==1).or.(kenergy==nenergy)))THEN
         WRITE(Prob_Unit,*)
         WRITE(Prob_Unit,'(1x,A,I5,A,ES15.7,A,I5,A,I5)') 'Kenergy=', Kenergy, ' Etot=', Etot, ' Nopen=', Nopen
         WRITE(Prob_Unit,'(4x,A,1x,A,1x,A,2x,A,2x,A,2x,A,4x,A,4x,A,4x,A)')'i','chanl', 'elect', 'nvib', 'jrot', 'lorb', 'AsymptEnergies(Ha)', 'AsymptEnergies(eV)', 'AsymptEnergies(cm^-1)'
         WRITE(NQuantStr,'(I1)') NQuant
         WRITE(NQuantFormatStr,*) '(1x,I4,'//NQuantStr//'I6,3ES22.7,I6)'
         DO i=1,nopen
            f=EminLoc(i)
            WRITE(Prob_Unit,NQuantFormatStr) i, (QNumbr(k,f),k=1,NQuant), AsymptEnergies(f), AsymptEnergies(f)*autoev, AsymptEnergies(f)/cmm1toau, f
         ENDDO
      ENDIF
      DEALLOCATE(AMatx, BMatx)
   ENDDO

   CALL Date_And_Time(today, hour, curzone, dtvalues)
   CALL Print_Date_And_Time(Today, Hour, CurZone, DtValues, Out_Unit)
   WRITE(AP_Unit,*)'Leaving:', ProcName

   CLOSE(In_UnitA)
   IF(AFile/=BFile)CLOSE(In_UnitB)
   CLOSE(Plt_Unit)

   !DEALLOCATE(ALIST, BLIST)
   !DEALLOCATE(QI_Elastic, QF_Elastic )
   !DEALLOCATE(QI_React, QF_React )
   !DEALLOCATE(QI_Inelastic, QF_Inelastic )
   ENDDO
ENDDO
ENDDO

DEALLOCATE(ALIST, BLIST)
DEALLOCATE(QI_Elastic, QF_Elastic )
DEALLOCATE(QI_React, QF_React )
DEALLOCATE(QI_Inelastic, QF_Inelastic )

IF(CFile=='Probability')THEN
   PrntQuant=.True.
   OPEN(Prob_Unit)
ENDIF

RETURN
ENDSUBROUTINE Asymptotic_Plotting
