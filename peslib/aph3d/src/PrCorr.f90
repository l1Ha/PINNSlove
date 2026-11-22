SUBROUTINE PrCorr(neig,nrho,track,elevels,rhovals,curve)
USE Numeric_Kinds_Module
USE Numbers_Module
USE FileUnits_Module
USE PES_MODULE
USE QAsy_Numbers_Module, ONLY : Nquant, Qnumbr, Qlabels, AsymptEnergies
USE Convrsns_Module
IMPLICIT NONE
INTEGER neig, nrho, track(neig,NRho), i, j, jrho, k1, k2, jj, MinInt, MaxInt, jquant
REAL(Kind=WP_Kind) curve(neig,NRho), Xmin, Xmax
REAL(Kind=WP_Kind) elevels(neig,NRho), rhovals(NRho)
CHARACTER(LEN=80) Title, SubTitle, XLabel, YLabel, MTVFile

WRITE(Out_Unit,*)'Printing plot DATA'
OPEN(Unit=WgPairsDiabatic_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'GraphicsOut/WgPairsDiabatic_Unit.rbw')
OPEN(Unit=WgReadyDiabatic_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'GraphicsOut/wgreadyDiabatic.rbw')
OPEN(Unit=DiabaticLevels_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'GraphicsOut/DiabaticLevels.graph')
OPEN(Unit=WgTableDiabatic_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'GraphicsOut/wgtableDiabatic.rbw')
OPEN(Unit=WgDiatomDiabatic_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'GraphicsOut/wgdiatomDiabatic.txt')
OPEN(Unit=DiabaticLevels_CSV_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'GraphicsOut/DiabaticLevels.csv')
OPEN(Unit=DiabaticLevelsInfo_CSV_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'GraphicsOut/DiabaticLevelsInfo.csv')
DO j=1,neig
   curve(j,1)=elevels(j,1)
   DO jrho=2,nrho
      k1=0
      k2=0
      DO jj=1,neig
         IF(track(jj,jrho-1)==j)k1=jj
         IF(track(jj,jrho)==j)k2=jj
      ENDDO
      IF(k1/=0.AND.k2/=0)THEN
         curve(j,jrho)=elevels(k2,jrho)
!------------------------------------------------------------
! WRITE DATA in Bob's pair FORMAT
!------------------------------------------------------------
         WRITE(WgPairsDiabatic_Unit,*)' $inp newcrv=.true.,ndata=2,'
         WRITE(WgPairsDiabatic_Unit,4)rhovals(jrho-1),rhovals(jrho)
         WRITE(WgPairsDiabatic_Unit,5)elevels(k1,jrho-1),elevels(k2,jrho)
         WRITE(WgPairsDiabatic_Unit,*)' $END'
      ENDIF
   ENDDO
ENDDO

OPEN(Unit=Track_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/curve-map.txt',status='unknown')
DO j=1,neig
    WRITE(Track_Unit,'("jeig=",I5," neig=",I5," nrho=",I5)')j,neig,nrho
    WRITE(Track_Unit,'("Track(jeig,jrho),jrho=1,nrho)")')
    WRITE(Track_Unit,'(20I5)')(track(j,jrho),jrho=1,nrho)
ENDDO
CLOSE(Unit=Track_Unit)

!------------------------------------------------------------
! WRITE DATA for xvgr
!------------------------------------------------------------
DO j=1,neig
   WRITE(DiabaticLevels_Unit,7)(rhovals(jrho),curve(j,jrho),jrho=1,nrho)
   IF(j/=neig)WRITE(DiabaticLevels_Unit,*)'&'
ENDDO

DO jrho=1,nrho
   WRITE(DiabaticLevels_CSV_Unit,'(1x,e14.7,100000(",",e14.7))')rhovals(jrho),(curve(j,jrho),j=1,neig)
ENDDO
WRITE(DiabaticLevelsInfo_CSV_Unit,'(1x,e14.7,100000(",",I8))')   rhovals(nrho),( Track(j,nrho),j=1,neig)
DO jquant=1,nquant !tmpmod gregparker
  WRITE(DiabaticLevelsInfo_CSV_Unit,'(1x,A14,100000(",",I8))')   TRIM(Qlabels(jquant)),(Qnumbr(jquant,Track(j,nrho)),j=1,neig)
ENDDO
WRITE(DiabaticLevelsInfo_CSV_Unit,'(1x,A14,100000(",",e14.7))')'Asymptenergies',(AsymptEnergies(Track(j,nrho))*autoev,j=1,neig)
WRITE(DiabaticLevelsInfo_CSV_Unit,'(1000(5I3,","))')((QNumbr(j,jj),j=1,nquant),jj=1,19) !QNumbr not allocated..Greg

!-------------------------------------------------------------
! WRITE DATA in Bob's readx ready FORMAT
!-------------------------------------------------------------
WRITE(WgReadyDiabatic_Unit,*)' $readx'
WRITE(WgReadyDiabatic_Unit,*)nrho
WRITE(WgReadyDiabatic_Unit,6)(rhovals(jrho),jrho=1,nrho)
DO j=1,neig
   WRITE(WgReadyDiabatic_Unit,*)' $ready'
   WRITE(WgReadyDiabatic_Unit,*)nrho
   WRITE(WgReadyDiabatic_Unit,6)(curve(j,jrho),jrho=1,nrho)
   WRITE(WgReadyDiabatic_Unit,*)'0$CONTINUE'
ENDDO


!---------------------------------------------------------------
! WRITE DATA in Bob's Table FORMAT
!---------------------------------------------------------------
WRITE(WgTableDiabatic_Unit,*)'0$readtb'
WRITE(WgTableDiabatic_Unit,*)neig,nrho,'  -1','  -1'
DO jrho=1,nrho
   WRITE(WgTableDiabatic_Unit,7)rhovals(jrho)
   WRITE(WgTableDiabatic_Unit,6)(curve(j,jrho),j=1,neig)
ENDDO
WRITE(WgTableDiabatic_Unit,*)'0$REWIND'

!-----------------------------------------------------------------
! WRITE DATA for Bob's Diatom code
!-----------------------------------------------------------------
WRITE(WgDiatomDiabatic_Unit,*)
WRITE(WgDiatomDiabatic_Unit,*)'      vshift=0.0'
WRITE(WgDiatomDiabatic_Unit,*)'vlargeR=',rhovals(nrho)
WRITE(WgDiatomDiabatic_Unit,*)
WRITE(WgDiatomDiabatic_Unit,*)'smallR=exp'
WRITE(WgDiatomDiabatic_Unit,*)'largeR=exp'
WRITE(WgDiatomDiabatic_Unit,*)
WRITE(WgDiatomDiabatic_Unit,*)'unit_distance=bohr'
WRITE(WgDiatomDiabatic_Unit,*)'unit_energy=ev'
WRITE(WgDiatomDiabatic_Unit,*)
WRITE(WgDiatomDiabatic_Unit,*)'fit=spline DATA=xypairs'
WRITE(WgDiatomDiabatic_Unit,*)
WRITE(WgDiatomDiabatic_Unit,*)' ! Distance (bohr)     Energy(eV)'
WRITE(WgDiatomDiabatic_Unit,*)'   --------------      ---------'
WRITE(WgDiatomDiabatic_Unit,*)
WRITE(WgDiatomDiabatic_Unit,*)'0$REWIND'
WRITE(WgDiatomDiabatic_Unit,7)(rhovals(jrho),curve(1,jrho),jrho=1,nrho)

MinInt=Rhovals(1)
MaxInt=Rhovals(NRho)
XMin=One*MinInt
XMax=One*MaxInt
Title=PES_Name
SubTitle='Diabatic Curves'
XLabel='Rho'
YLabel='Energy(eV)'
MTVFile='GraphicsOut/DiabaticLevels.mtv'
OPEN(Unit=DiaAdia_Bin_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/DiabaticLevels.bin',form='unformatted')
WRITE(DiaAdia_Bin_Unit)Subtitle
WRITE(DiaAdia_Bin_Unit)Nrho
WRITE(DiaAdia_Bin_Unit)(RhoVals(i),i=1,Nrho)
WRITE(DiaAdia_Bin_Unit)Neig
WRITE(DiaAdia_Bin_Unit)((Curve(j,i),i=1,nrho),j=1,neig)
CLOSE(DiaAdia_Bin_Unit)
CALL Graf_2d(NRho, Neig, Rhovals, curve, Xmin, Xmax, Title, SubTitle, XLabel, YLabel, MTVFile)

SubTitle='Adiabatic Curves'
MTVFile='GraphicsOut/AdiabaticLevels.mtv'
OPEN(Unit=DiaAdia_Bin_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/AdiabaticLevels.bin',form='unformatted')
WRITE(DiaAdia_Bin_Unit)Subtitle
WRITE(DiaAdia_Bin_Unit)Nrho
WRITE(DiaAdia_Bin_Unit)(RhoVals(i),i=1,Nrho)
WRITE(DiaAdia_Bin_Unit)Neig
WRITE(DiaAdia_Bin_Unit)((Elevels(j,i),i=1,nrho),j=1,neig)
CLOSE(DiaAdia_Bin_Unit)
CALL Graf_2d(NRho, Neig, Rhovals, elevels, Xmin, Xmax, Title, SubTitle, XLabel, YLabel, MTVFile)

!----------------------------------------------------------------
! CLOSE units
!----------------------------------------------------------------
CLOSE(Unit=WgPairsDiabatic_Unit)
CLOSE(Unit=WgReadyDiabatic_Unit)
CLOSE(Unit=DiabaticLevels_Unit)
CLOSE(Unit=WgTableDiabatic_Unit)

WRITE(Out_Unit,*)'Completed writing plot DATA'

4  FORMAT(1x,'   xdata(1)=',2e15.7)
5  FORMAT(1x,'   ydata(1)=',2e15.7)
6  FORMAT(1x,5e15.7)
7  FORMAT(1x,2e15.7)

RETURN
ENDSUBROUTINE PrCorr
