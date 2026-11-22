SUBROUTINE CompareMethods
USE Numeric_Kinds_Module
USE Numbers_Module
USE FileUnits_Module
USE FileUnits_OneDim_Module
USE OneDim_Parms_Module
USE DiatomicPot_MODULE
USE Masses_Module
USE Convrsns_Module
USE Region_Module
USE Oops_Module
USE Quantb_module, ONLY: JminFull=>Jmin, JmaxFull=>Jmax, MinVib, MaxVib   !, MinVib=>MinVibFull, MaxVib=>MaxVibFull
IMPLICIT NONE
LOGICAL, PARAMETER :: TmpJacobi=.TRUE.
INTEGER(KIND=IW_Kind) NBound, NQuasi, NStates, NumNu, nYTerms, i, istate, nup, jval, ip1, ip2, ip3, iread, NuTemp, jTemp, elecTemp
INTEGER(KIND=IW_Kind) l3, m3, k3, n3, Index, nt, nchanacc, ntmax, karran, nstat(3), nnuj(3), nnuj_nhermt, nhermt_iarran, MinNStates, nvibmax, nrotmax
INTEGER(KIND=IW_Kind), ALLOCATABLE:: j_Daf(:),    j_Numerov(:),    j_FourierGrid(:)
INTEGER(KIND=IW_Kind), ALLOCATABLE:: elec_Daf(:), elec_Numerov(:), elec_FourierGrid(:)
INTEGER(KIND=IW_Kind), ALLOCATABLE:: Nu_Daf(:),   Nu_Numerov(:),   Nu_FourierGrid(:)
INTEGER(KIND=IW_Kind), ALLOCATABLE:: j_Harmonic(:,:), nu_Harmonic(:,:)
REAL(KIND=WP_Kind), ALLOCATABLE:: YfitDaf(:,:),         EigDaf(:),         PsiDaf(:,:)
REAL(KIND=WP_Kind), ALLOCATABLE:: YfitNumerov(:,:),     EigNumerov(:),     PsiNumerov(:,:)
REAL(KIND=WP_Kind), ALLOCATABLE:: YfitFourierGrid(:,:), EigFourierGrid(:), PsiFourierGrid(:,:)
REAL(KIND=WP_Kind), ALLOCATABLE:: PsiHarmonic(:,:), EigHarmonic(:,:,:)
REAL(KIND=WP_Kind) iread2, pdiff, x3, XsqHarmonic, XkHarmonic, junk

OPEN(Unit=9945,File=TRIM(Outdir)//"DiatomicOut/CompareMethods(Ha).txt")
OPEN(Unit=9946,File=TRIM(Outdir)//"DiatomicOut/CompareMethods(eV).txt")

DO IARRAN=1,3
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
   DiatomicPot=TRIM(AtomicSymbol(ip1))//TRIM(AtomicSymbol(ip2))

! The DAF method
Nyterms=MaxYterms
NStates=NStatesDAF(IARRAN)
ALLOCATE(YFitDaf(jmin:jmax,0:NYterms-1))
ALLOCATE(EigDaf(0:NStates-1))
ALLOCATE(PsiDaf(0:NMax,0:NStates-1))
ALLOCATE(j_Daf(0:NStates-1))
ALLOCATE(Nu_Daf(0:NStates-1))
ALLOCATE(Elec_Daf(0:NStates-1))
elec_Daf=0
j_Daf=0
Nu_Daf=0
EigDaf=0.d0
YFitDaf=0.d0
PsiDaf=0.d0
OPEN(Unit=DAFBin_Unit,file=TRIM(OutDIR)//"BinOut/"//TRIM(DiatomicPot)//"_DAFBin.bin",form="unformatted")
REWIND(DAFBin_Unit)
YfitDaf=Zero
IState=0
NStates=NStatesDAF(IARRAN)
DO WHILE(IState<NStates)
   READ(DAFBin_Unit)elec_Daf(IState), j_Daf(IState), numnu, nbound, nquasi, nYTerms
   jval=j_Daf(IState)
   DO Nup=0,NumNu-1
      READ(DAFBin_Unit)Nu_Daf(IState), EigDaf(IState)
      READ(DAFBin_Unit)PsiDaf(:,IState)
      IState=IState+1
      IF(IState>=NStates)EXIT
      elec_Daf(IState)=elec_Daf(IState-1)
      j_Daf(IState)=j_Daf(IState-1)
   ENDDO
   READ(DAFBin_Unit)(YFitDaf(jval,i),i=0,NYterms-1)
ENDDO
CLOSE(Unit=DAFBin_Unit)

! The Numerov method
NStates=NStatesNumerov(IARRAN)
Nyterms=MaxYterms
ALLOCATE(YFitNumerov(jmin:jmax,0:NYterms-1))
ALLOCATE(EigNumerov(0:NStates-1))
ALLOCATE(PsiNumerov(0:NMax,0:NStates-1))
ALLOCATE(j_Numerov(0:NStates-1))
ALLOCATE(nu_Numerov(0:NStates-1))
ALLOCATE(elec_Numerov(0:NStates-1))
elec_Numerov=0
j_Numerov=0
Nu_Numerov=0
EigNumerov=0.d0
YFitNumerov=0.d0
PsiNumerov=0.d0
OPEN(Unit=NumerovBin_Unit,file=TRIM(OutDIR)//"BinOut/"//TRIM(DiatomicPot)//"_NumerovBin.bin",form="unformatted")
REWIND(NumerovBin_Unit)
YfitNumerov=Zero
IState=0
DO IState=0,NStates-1
   READ(NumerovBin_Unit)NuTemp, jTemp, elecTemp
   j_Numerov(IState)=jTemp
   Nu_Numerov(IState)=NuTemp
   IF(NuTemp<0)THEN
      READ(NumerovBin_Unit)elec_Numerov(IState), j_Numerov(IState), numnu, nbound, nquasi, nYterms
      jval=j_Numerov(IState)
      READ(NumerovBin_Unit)(YFitNumerov(jval,i),i=0,nYterms-1)
      READ(NumerovBin_Unit)nu_Numerov(IState), j_Numerov(IState), elec_Numerov(IState)
   ELSEIF(IState==NStates-1)THEN
      READ(NumerovBin_Unit)EigNumerov(IState)
      READ(NumerovBin_Unit)PsiNumerov(:,IState)
      READ(NumerovBin_Unit)NuTemp
      EXIT
   ENDIF
   jval=j_Numerov(IState)
   READ(NumerovBin_Unit)EigNumerov(IState)
   READ(NumerovBin_Unit)PsiNumerov(:,IState)
ENDDO
IState=NStates-1
READ(NumerovBin_Unit)elec_Numerov(IState), j_Numerov(IState), numnu, nbound, nquasi, nYterms !tmpmod gregparker is this correct?
jval=j_Numerov(IState)
READ(NumerovBin_Unit)(YFitNumerov(jval,i),i=0,nYterms-1)
CLOSE(Unit=NumerovBin_Unit)

! The FourierGrid method
NStates=NStatesFourierGrid(IARRAN)
Nyterms=MaxYterms
ALLOCATE(YFitFourierGrid(jmin:jmax,0:NYterms-1))
ALLOCATE(EigFourierGrid(0:NStates-1))
ALLOCATE(PsiFourierGrid(0:NMax,0:NStates-1))
ALLOCATE(j_FourierGrid(0:NStates-1))
ALLOCATE(nu_FourierGrid(0:NStates-1))
ALLOCATE(elec_FourierGrid(0:NStates-1))
elec_FourierGrid=0
j_FourierGrid=0
Nu_FourierGrid=0
EigFourierGrid=0.d0
YFitFourierGrid=0.d0
PsiFourierGrid=0.d0
OPEN(Unit=FourierGridBin_Unit,file=TRIM(OutDIR)//"BinOut/"//TRIM(DiatomicPot)//"_FourierGridBin.bin",form="unformatted")
REWIND(FourierGridBin_Unit)
YfitFourierGrid=Zero
IState=0
DO WHILE(IState<NStates)
   READ(FourierGridBin_Unit)elec_FourierGrid(IState), j_FourierGrid(IState), numnu, nbound, nquasi, nYTerms
   jval=j_FourierGrid(IState)
   DO Nup=0,NumNu-1
      READ(FourierGridBin_Unit)Nu_FourierGrid(IState), EigFourierGrid(IState)
      READ(FourierGridBin_Unit)PsiFourierGrid(:,IState)
      IState=IState+1
      IF(IState>=NStates)EXIT
      elec_FourierGrid(IState)=elec_FourierGrid(IState-1)
      j_FourierGrid(IState)=j_FourierGrid(IState-1)
   ENDDO
   READ(FourierGridBin_Unit)(YFitFourierGrid(jval,i),i=0,NYterms-1)
ENDDO
CLOSE(Unit=FourierGridBin_Unit)

IF(TmpJacobi)THEN
! Harmonic oscillator method used in Jacobi
IF(IARRAN==1)THEN
   iregion='jacobi'
   ioops=.false.
   WRITE(*,*)"iregion=", iregion,"  ioops=", ioops, "Loc2"
   OPEN(Unit=bas_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Jac_Basis.bin', form = 'unformatted', status = 'Old')
   REWIND(Bas_Unit)
   DO karran=1,3
      READ(bas_unit)nstat(karran)
      READ(bas_unit)
      READ(bas_unit)
      READ(bas_unit)
      READ(bas_unit)
   ENDDO
   ntmax=SUM(nstat)
   REWIND(Bas_Unit)
ENDIF
READ(bas_unit)nt, nnuj(IARRAN), nchanacc, nnuj_nhermt, nhermt_iarran, ioops, iregion
!WRITE(*,'(6A15)')"IARRAN","nt", "nnuj(IARRAN)", "nchanacc", "nnuj_nhermt", "nhermt_iarran"
!WRITE(*,'(6I15)')IARRAN, nt, nnuj(IARRAN), nchanacc, nnuj_nhermt, nhermt_iarran
WRITE(9945,*)
WRITE(9945,*)
WRITE(9945,*)
WRITE(9945,*)DiatomicPot, ioops, iregion
WRITE(9945,'(6A15)')"IARRAN","nt", "nnuj(IARRAN)", "nchanacc", "nnuj_nhermt", "nhermt_iarran"
WRITE(9945,'(6I15)')IARRAN, nt, nnuj(IARRAN), nchanacc, nnuj_nhermt, nhermt_iarran
WRITE(9946,*)
WRITE(9946,*)
WRITE(9946,*)
WRITE(9946,*)DiatomicPot
WRITE(9946,'(6A15)')"IARRAN","nt", "nnuj(IARRAN)", "nchanacc", "nnuj_nhermt", "nhermt_iarran"
WRITE(9946,'(6I15)')IARRAN, nt, nnuj(IARRAN), nchanacc, nnuj_nhermt, nhermt_iarran
IF(IARRAN==1)THEN 
   ALLOCATE(j_Harmonic(0:ntmax-1,1:narran))
   ALLOCATE(nu_Harmonic(0:ntmax-1,1:narran))
   ALLOCATE(PsiHarmonic(0:NMax,0:ntmax-1))
   nvibmax=MAXVAL(MaxVib)
   nrotmax=MAXVAL(jmaxFull)
   ALLOCATE(EigHarmonic(0:nvibmax,0:nrotmax,1:Narran))   !tmpmodGregParker fix this problem
   j_Harmonic=0
   nu_Harmonic=0
   PsiHarmonic=0.d0
   EigHarmonic=0.d0
ENDIF
READ(bas_unit)
READ(bas_unit)(Nu_Harmonic(i,IARRAN), j_Harmonic(i,IARRAN), l3, m3, k3, x3, n3, i=0,nt-1)
IF(IARRAN==1)THEN
   Index=0
ELSEIF(IARRAN==2)THEN
   Index=nnuj(1)
ELSEIF(IARRAN==3)THEN
   Index=nnuj(1)+nnuj(2)
ENDIF
   READ(bas_unit)!(PsiHarmonic(Index+i), i=1,nnuj(IARRAN))
   READ(bas_unit)(XsqHarmonic,XkHarmonic,EigHarmonic(Nu_Harmonic(i,IARRAN),j_Harmonic(i,IARRAN),IARRAN), i=0,nt-1)
ENDIF
   

! Output the comparison chart in Ha
MinNStates=MIN(NStatesDAF(IARRAN),NStatesNumerov(IARRAN),NStatesFourierGrid(IARRAN))
DO IState=0,MinNStates-1
   IF(IState==0)WRITE(9945,'(4A10,5A23)')"IState", "Elec_Daf", "Nu_Daf", "j_Daf", "EigDaf(Ha)", "EigNumerov(Ha)", "EigFourierGrid(Ha)","EigHarmonic(Ha)","Percent Diff"
   IF(Nu_Daf(IState)<=MaxVib(IARRAN))THEN
      IF(elec_Daf(IState)/=elec_Numerov(IState).or.elec_Numerov(IState)/=elec_FourierGrid(IState))THEN
         WRITE(9945,'(I5,A,4I5)')IState, "  Elec=", elec_Daf(IState), elec_Numerov(IState),   elec_FourierGrid(IState)
      ENDIF
      IF(j_Daf(IState)/=j_Numerov(IState).or.j_Numerov(IState)/=j_FourierGrid(IState))THEN
         WRITE(9945,'(I5,A,4I5)')IState, "  jrot=", j_Daf(IState),    j_Numerov(IState),      j_FourierGrid(IState)
      ENDIF
      IF(Nu_Daf(IState)/=nu_Numerov(IState).or.nu_Numerov(IState)/=nu_FourierGrid(IState))THEN
         WRITE(9945,'(I5,A,4I5)')IState, "  mvib=", Nu_Daf(IState),   Nu_Numerov(IState),     Nu_FourierGrid(IState)
      ENDIF
      pdiff=100.d0*MAX(ABS(EigDaf(IState)-EigNumerov(IState)),ABS(EigDaf(IState)-EigFourierGrid(IState)),ABS(EigDaf(IState)-EigHarmonic(Nu_Daf(IState),j_Daf(IState),IARRAN)))/EigDaf(IState)
      WRITE(9945,'(4I10,4ES23.15,F23.7)')IState, Elec_Daf(IState), Nu_Daf(IState), j_Daf(IState), EigDaf(IState), EigNumerov(IState), EigFourierGrid(IState), EigHarmonic(Nu_Daf(IState),j_Daf(IState),IARRAN),pdiff
   ENDIF
ENDDO

! Dunham energy fit parameters in Ha
WRITE(9945,*)"Dunham energy fit parameters in Ha"
DO jval=jmin,jmax
   WRITE(9945,*)
   WRITE(9945,'(A,I5)')"jval=",jval
   DO i=0,NYterms-1
      WRITE(9945,'(2I5,"  YFit=",3ES23.15)')i,jval,YFitDAF(jval,i), YFitNumerov(jval,i), YFitFourierGrid(jval,i)
   ENDDO
ENDDO


! Output the comparison chart in eV
DO IState=0,MinNStates-1
   IF(IState==0)WRITE(9946,'(4A10,5A23)')"IState", "Elec_Daf", "Nu_Daf", "j_Daf", "EigDaf(eV)", "EigNumerov(eV)", "EigFourierGrid(eV)","EigHarmonic(eV)","Percent Diff"
   IF(Nu_Daf(IState)<=MaxVib(IARRAN))THEN
      IF(elec_Daf(IState)/=elec_Numerov(IState).or.elec_Numerov(IState)/=elec_FourierGrid(IState))THEN
         WRITE(9946,'(I5,A,4I5)')IState, "  Elec=", elec_Daf(IState), elec_Numerov(IState),   elec_FourierGrid(IState)
      ENDIF
      IF(j_Daf(IState)/=j_Numerov(IState).or.j_Numerov(IState)/=j_FourierGrid(IState))THEN
         WRITE(9946,'(I5,A,4I5)')IState, "  jrot=", j_Daf(IState),    j_Numerov(IState),      j_FourierGrid(IState)
      ENDIF
      IF(Nu_Daf(IState)/=nu_Numerov(IState).or.nu_Numerov(IState)/=nu_FourierGrid(IState))THEN
         WRITE(9946,'(I5,A,4I5)')IState, "  mvib=", Nu_Daf(IState),   Nu_Numerov(IState),     Nu_FourierGrid(IState)
      ENDIF
      pdiff=100.d0*MAX(ABS(EigDaf(IState)-EigNumerov(IState)),ABS(EigDaf(IState)-EigFourierGrid(IState)),ABS(EigDaf(IState)-EigHarmonic(Nu_Daf(IState),j_Daf(IState),IARRAN)))/EigDaf(IState)
      WRITE(9946,'(4I10,4ES23.15,F23.7)')IState, Elec_Daf(IState), Nu_Daf(IState), j_Daf(IState), EigDaf(IState)*autoev, EigNumerov(IState)*autoev, EigFourierGrid(IState)*autoev, EigHarmonic(Nu_Daf(IState),j_Daf(IState),IARRAN)*autoev,pdiff
   ENDIF
ENDDO

! Dunham energy fit parameters in eV
WRITE(9945,*)"Dunham energy fit parameters in eV"
DO jval=jmin,jmax
   WRITE(9946,*)
   WRITE(9946,'(A,I5)')"jval=",jval
   DO i=0,NYterms-1
      WRITE(9946,'(2I5,"  YFit=",3ES23.15)')i,jval,YFitDAF(jval,i)*autoev, YFitNumerov(jval,i)*autoev, YFitFourierGrid(jval,i)*autoev
   ENDDO
ENDDO

DEALLOCATE(YFitDaf)
DEALLOCATE(EigDaf)
DEALLOCATE(PsiDaf)
DEALLOCATE(j_Daf)
DEALLOCATE(Nu_Daf)
DEALLOCATE(Elec_Daf)

DEALLOCATE(YFitNumerov)
DEALLOCATE(EigNumerov)
DEALLOCATE(PsiNumerov)
DEALLOCATE(j_Numerov)
DEALLOCATE(nu_Numerov)
DEALLOCATE(elec_Numerov)

DEALLOCATE(YFitFourierGrid)
DEALLOCATE(EigFourierGrid)
DEALLOCATE(PsiFourierGrid)
DEALLOCATE(j_FourierGrid)
DEALLOCATE(nu_FourierGrid)
DEALLOCATE(elec_FourierGrid)

IF(TmpJacobi.and.IARRAN==3)THEN
DEALLOCATE(j_Harmonic)
DEALLOCATE(nu_Harmonic)
DEALLOCATE(PsiHarmonic)
DEALLOCATE(EigHarmonic)
ENDIF


ENDDO

CLOSE(Unit=9945)
CLOSE(Unit=9946)

RETURN
ENDSUBROUTINE CompareMethods