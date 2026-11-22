SUBROUTINE YFitter
USE Numeric_Kinds_Module
USE Numbers_Module
!USE CriticalDist_Module
USE FileUnits_Module
USE FileUnits_OneDim_Module
USE Qstate_Module
USE OneDim_Parms_Module, ONLY: MaxYterms, IARRAN, jmin, jmax, MinElec, MaxElec
USE ElectronicState_Module
USE DiatomicPot_Module, ONLY: DiatomicPot
USE Quantb_module, ONLY: jMinVal, jMaxVal, NuMinVal, NuMaxVal, MinVib, MaxVib
USE Convrsns_Module
USE Spectro_Module
USE Numbers_Module, ONLY: Zero
IMPLICIT NONE
INTEGER(KIND=IW_Kind) IMethod, NYTerms, numnu, Nbound, NQuasi, i, j, info, kterm, nyterms2
INTEGER(KIND=IW_Kind), PARAMETER :: lwork=1000
CHARACTER(LEN=11) Method
CHARACTER(LEN=25) PotID
! PORT FIX (2026-10-01): YFitter is called once per arrangement (IARRAN=1..3)
! but allocates its work arrays only on the IARRAN==1 call; without SAVE they
! were deallocated on return and the IARRAN=2/3 calls wrote through dangling
! descriptors (SIGSEGV at line 65).
REAL(KIND=WP_Kind), ALLOCATABLE, SAVE :: YFitDAF(:,:,:), YFitNumerov(:,:,:), YFitFourier(:,:,:)
REAL(KIND=WP_Kind), ALLOCATABLE, SAVE :: work(:), q(:,:), amat(:,:), yFit(:), yFitj(:), Eng(:)
REAL(KIND=WP_Kind), ALLOCATABLE, SAVE :: Yjk_Daf(:,:), Yjk_Numerov(:,:), Yjk_Fourier(:,:), Yjk(:,:,:)
REAL(KIND=WP_Kind), ALLOCATABLE, SAVE :: Spectro_Eng_Daf(:,:,:,:), Spectro_Eng_Numerov(:,:,:,:), Spectro_Eng_Fourier(:,:,:,:), Spectro_Eng(:,:,:,:)

nYterms=MIN(MaxYterms,NuMaxVal)

IF(IARRAN==1)THEN
   ALLOCATE(Yfit(MaxYterms))
   ALLOCATE(yFitj(0:jMaxVal))
   ALLOCATE(Amat(0:jMaxVal,MaxYterms))
   REWIND(YfitBin_Unit)
   ALLOCATE(YFitDaf(0:jMaxVal,0:MaxYterms-1,1:3))
   ALLOCATE(YFitNumerov(0:jMaxVal,0:MaxYterms-1,1:3))
   ALLOCATE(YFitFourier(0:jMaxVal,0:MaxYterms-1,1:3))
   ALLOCATE(Work(lwork))
   ALLOCATE(Yjk(0:MaxYterms-1,0:MaxYterms-1, NArran))
   ALLOCATE(Yjk_Daf(0:MaxYterms-1,0:MaxYterms-1))
   ALLOCATE(Yjk_Numerov(0:MaxYterms-1,0:MaxYterms-1))
   ALLOCATE(Yjk_Fourier(0:MaxYterms-1,0:MaxYterms-1))
   ALLOCATE(Spectro_Eng(MinElec:MaxElec,0:NuMaxVal,0:jMaxVal,1:NArran))
   ALLOCATE(Spectro_Eng_Daf(MinElec:MaxElec,0:NuMaxVal,0:jMaxVal,1:NArran))
   ALLOCATE(Spectro_Eng_Numerov(MinElec:MaxElec,0:NuMaxVal,0:jMaxVal,1:NArran))
   ALLOCATE(Spectro_Eng_Fourier(MinElec:MaxElec,0:NuMaxVal,0:jMaxVal,1:NArran))
   ALLOCATE(Eng(0:NuMaxVal))
   Yjk=Zero
   Yjk_Daf=Zero
   Yjk_Numerov=Zero
   Yjk_fourier=Zero
   YFitDaf=Zero
   YFitNumerov=Zero
   YFitFourier=Zero
ENDIF

WRITE(YFitter_Unit,*)
WRITE(YFitter_Unit,'(2A,I5)')DiatomicPot, "  IArran=", Iarran
DO j=jmin,jmax
   WRITE(YFitter_Unit,*)
   WRITE(YFitter_Unit,*)"j=",j
   DO IMethod=1,3
     READ(YfitBin_Unit)elec,jval,numnu,nbound,nquasi,nYTerms2
     READ(YfitBin_Unit)Method, PotID, Iarran, elec, jval, (Yfit(i+1),i=0,nYterms2-1)
     IF(IMethod==1.and.TRIM(Method)=="DAF")THEN
        WRITE(YFitter_Unit,'(A11,5x,A5,3I5,5ES15.7)')TRIM(Method),TRIM(PotID), Iarran, elec, jval, (Yfit(i+1),i=0,nYterms-1)
        DO i=0,nYterms-1
           YFitDaf(jval, i, IArran)=yfit(I+1)
        ENDDO
      ENDIF
      IF(IMethod==2.and.TRIM(Method)=="Numerov")THEN
         WRITE(YFitter_Unit,'(A11,5x,A5,3I5,5ES15.7)')TRIM(Method),TRIM(PotID), Iarran, elec, jval, (Yfit(i+1),i=0,nYterms-1)
         DO i=0,nYterms-1
            YFitNumerov(jval, i, IArran)=yfit(I+1)
         ENDDO
      ENDIF
      IF(IMethod==3.and.TRIM(Method)=="FourierGrid")THEN
         WRITE(YFitter_Unit,'(A11,5x,A5,3I5,5ES15.7)')TRIM(Method),TRIM(PotID), Iarran, elec, jval, (Yfit(i+1),i=0,nYterms-1)
         DO i=0,nYterms-1
            YFitFourier(jval, i, IArran)=yfit(I+1)
         ENDDO
      ENDIF
   ENDDO
ENDDO

WRITE(YFitter_Unit,'(A)')"DAF Dunham Y Parameters"
WRITE(YFitter_Unit,'(2A,I5)')DiatomicPot, "  IArran=", Iarran
DO j=jmin,jmax
   WRITE(YFitter_Unit,'(1x,I5,2x,5ES23.15)')j,(YfitDAF(j,i,IArran)/cmm1toau,i=0,nYterms-1)
ENDDO

WRITE(YFitter_Unit,*)
WRITE(YFitter_Unit,'(A)')"Numerov Dunham Y Parameters"
WRITE(YFitter_Unit,'(2A,I5)')DiatomicPot, "  IArran=", Iarran
DO j=jmin,jmax
   WRITE(YFitter_Unit,'(1x,I5,2x,5ES23.15)')j,(YfitNumerov(j,i,IArran)/cmm1toau,i=0,nYterms-1)
ENDDO

WRITE(YFitter_Unit,*)
WRITE(YFitter_Unit,'(A)')"FourierGrid Dunham Y Parameters"
WRITE(YFitter_Unit,'(2A,I5)')DiatomicPot, "  IArran=", Iarran
DO j=jmin,jmax
   WRITE(YFitter_Unit,'(1x,I5,2x,5ES23.15)')j,(YfitFourier(j,i,IArran)/cmm1toau,i=0,nYterms-1)
ENDDO

WRITE(YFitter_Unit,*)
WRITE(YFitter_Unit,'(2A,I5)')DiatomicPot, "  IArran=", Iarran
WRITE(YFitter_Unit,'(A)')"DAF Dunham Y(j,k) Parameters"
WRITE(YFitter_Unit,'(4x,5(A16,2x,I5))')("j=",j,j=0,NYterms-1)
yfitj=Zero
DO i=0,NYterms-1
   DO j=jmin,jmax
      yfitj(j)=YFitDaf(j,i,IArran)
   ENDDO
   DO kterm=0,nYTerms-1
      DO j=0,jMaxVal
         amat(j,kterm+1)=(j*(j+1))**kterm
      ENDDO
   ENDDO
   nyterms2=MIN(jMaxVal+1,nyterms)
   CALL dgels('N',jmax+1,nYterms2,1,amat,jMaxVal+1,Yfitj,jMaxVal+1,work,lwork,info)
   WRITE(YFitter_Unit,'(A,1x,I5,2x,5ES23.15)')"k=",i,(Yfitj(kterm)/cmm1toau,kterm=0,nYterms2-1)
   DO kterm=0,nYterms2-1
      Yjk_Daf(i,kterm)=yfitj(kterm)
   ENDDO
ENDDO

WRITE(YFitter_Unit,*)
WRITE(YFitter_Unit,'(A)')"Numerov Dunham Y(j,k) Parameters"
WRITE(YFitter_Unit,'(4x,5(A16,2x,I5))')("j=",j,j=0,NYterms-1)
yfitj=Zero
DO i=0,NYterms-1
   DO j=jmin,jmax
      yfitj(j)=YFitNumerov(j,i,IArran)
   ENDDO
   DO kterm=0,nYTerms-1
      DO j=0,jMaxVal
         amat(j,kterm+1)=(j*(j+1))**kterm
      ENDDO
   ENDDO
   nyterms2=MIN(jMaxVal+1,nyterms)
   CALL dgels('N',jmax+1,nYterms2,1,amat,jMaxVal+1,Yfitj,jMaxVal+1,work,lwork,info)
   WRITE(YFitter_Unit,'(A,1x,I5,2x,5ES23.15)')"k=",i,(Yfitj(kterm)/cmm1toau,kterm=0,nYterms2-1)
   DO kterm=0,nYterms2-1
      Yjk_Numerov(i,kterm)=yfitj(kterm)
   ENDDO
ENDDO

WRITE(YFitter_Unit,*)
WRITE(YFitter_Unit,'(A)')"FourierGrid Dunham Y(j,k) Parameters"
WRITE(YFitter_Unit,'(4x,5(A16,2x,I5))')("j=",j,j=0,NYterms-1)
yfitj=Zero
DO i=0,NYterms-1
   DO j=jmin,jmax
      yfitj(j)=YFitFourier(j,i,IArran)
   ENDDO
   DO kterm=0,nYTerms-1
      DO j=0,jMaxVal
         amat(j,kterm+1)=(j*(j+1))**kterm
         IF(kterm==0)amat(j,kterm+1)=One
      ENDDO
   ENDDO
   nyterms2=MIN(jMaxVal+1,nyterms)
   CALL dgels('N',jmax+1,nYterms2,1,amat,jMaxVal+1,Yfitj,jMaxVal+1,work,lwork,info)
   WRITE(YFitter_Unit,'(A,1x,I5,2x,5ES23.15)')"k=",i,(Yfitj(kterm)/cmm1toau,kterm=0,nYterms2-1)
   DO kterm=0,nYterms2-1
      Yjk_Fourier(i,kterm)=yfitj(kterm)
   ENDDO
ENDDO

!Generated spectra from the DAF Dunham Yjk parameters
WRITE(YFitter_Unit,*)
WRITE(YFitter_Unit,*)"Energies using Yjk_Daf Expansion"
WRITE(YFitter_Unit,'(4x,100(A16,2x,I5))')("nu=",nu,nu=0,NuMaxVal)
DO j=jmin,jmax
   DO nu=NuMinVal,NuMaxVal
      Eng(nu)=Zero
      DO i=0,nYterms-1
         DO Kterm=0,nYterms-1
            Eng(nu)=Eng(nu)+Yjk_Daf(kterm,i)*((j*(j+1))**i)*((nu+half)**kterm)
         ENDDO
      ENDDO
      Spectro_Eng_Daf(elec,nu,j,IArran)=Eng(nu)
   ENDDO
   WRITE(YFitter_Unit,'(A,1x,I5,2x,100ES23.15)')"j=",j,(Eng(nu)/cmm1toau,nu=NuMinVal,NuMaxVal)
ENDDO

!Generated spectra from the Numerov Dunham Yjk parameters
WRITE(YFitter_Unit,*)
WRITE(YFitter_Unit,*)"Energies using Yjk_Numerov Expansion"
WRITE(YFitter_Unit,'(4x,100(A16,2x,I5))')("nu=",nu,nu=0,NuMaxVal)
DO j=jmin,jmax
   DO nu=NuMinVal,NuMaxVal
      Eng(nu)=Zero
      DO i=0,nYterms-1
         DO Kterm=0,nYterms-1
            Eng(nu)=Eng(nu)+Yjk_Numerov(kterm,i)*((j*(j+1))**i)*((nu+half)**kterm)
         ENDDO
      ENDDO
      Spectro_Eng_Numerov(elec,nu,j,IArran)=Eng(nu)
   ENDDO
   WRITE(YFitter_Unit,'(A,1x,I5,2x,100ES23.15)')"j=",j,(Eng(nu)/cmm1toau,nu=NuMinVal,NuMaxVal)
ENDDO

!Generated spectra from the Fourier Dunham Yjk parameters
WRITE(YFitter_Unit,*)
WRITE(YFitter_Unit,*)"Energies using Yjk_Fourier Expansion"
WRITE(YFitter_Unit,'(4x,100(A16,2x,I5))')("nu=",nu,nu=0,NuMaxVal)
DO j=jmin,jmax
   DO nu=NuMinVal,NuMaxVal
      Eng(nu)=Zero
      DO i=0,nYterms-1
         DO Kterm=0,nYterms-1
            Eng(nu)=Eng(nu)+Yjk_Fourier(kterm,i)*((j*(j+1))**i)*((nu+half)**kterm)
         ENDDO
      ENDDO
      Spectro_Eng_Fourier(elec,nu,j,IArran)=Eng(nu)
   ENDDO
   WRITE(YFitter_Unit,'(A,1x,I5,2x,100ES23.15)')"j=",j,(Eng(nu)/cmm1toau,nu=NuMinVal,NuMaxVal)
ENDDO

Yjk(:,:,IArran)=Yjk_Daf
Spectro_Eng=Spectro_Eng_DAF

!Use input data to set default
we(IArran)=wecm(IArran)*cmm1toau
wexe(IArran)=wexecm(IArran)*cmm1toau
weye(IArran)=weyecm(IArran)*cmm1toau
be(IArran)=becm(IArran)*cmm1toau
de(IArran)=decm(IArran)*cmm1toau
alfe(IArran)=alfecm(IArran)*cmm1toau
!These constants are usually not read in so these may be zero
wzero(IArran)=wzerocm(IArran)*cmm1toau
bete(IArran)=betecm(IArran)*cmm1toau
game(IArran)=gamecm(IArran)*cmm1toau
he(IArran)=hecm(IArran)*cmm1toau

IF(MaxVib(IArran)>=30000)THEN
wzero(IArran)=Yjk(0,0,IArran)
we(IArran)=Yjk(1,0,IArran)
wexe(IArran)=-Yjk(2,0,IArran)
weye(IArran)=Yjk(3,0,IArran)
be(IArran)=Yjk(0,1,IArran)
de(IArran)=Yjk(0,2,IArran)
he(IArran)=Yjk(0,3,IArran)
alfe(IArran)=-Yjk(1,1,IArran)
bete(IArran)=Yjk(1,2,IArran)
game(IArran)=Yjk(2,1,IArran)
!now in inverse centemeters
wzerocm(IArran)=Yjk(0,0,IArran)/cmm1toau
wecm(IArran)=Yjk(1,0,IArran)/cmm1toau
wexecm(IArran)=-Yjk(2,0,IArran)/cmm1toau
weyecm(IArran)=Yjk(3,0,IArran)/cmm1toau
becm(IArran)=Yjk(0,1,IArran)/cmm1toau
decm(IArran)=Yjk(0,2,IArran)/cmm1toau
hecm(IArran)=Yjk(0,3,IArran)/cmm1toau
alfecm(IArran)=-Yjk(1,1,IArran)/cmm1toau
betecm(IArran)=Yjk(1,2,IArran)/cmm1toau
gamecm(IArran)=Yjk(2,1,IArran)/cmm1toau
ELSEIF(MaxVib(IArran)>=2)THEN
wzero(IArran)=Yjk(0,0,IArran)
we(IArran)=Yjk(1,0,IArran)
wexe(IArran)=-Yjk(2,0,IArran)
weye(IArran)=Yjk(3,0,IArran)
be(IArran)=Yjk(0,1,IArran)
de(IArran)=Yjk(0,2,IArran)
he(IArran)=Yjk(0,3,IArran)
alfe(IArran)=-Yjk(1,1,IArran)
bete(IArran)=Yjk(1,2,IArran)
game(IArran)=Yjk(2,1,IArran)
!now in inverse centemeters
wzerocm(IArran)=Yjk(0,0,IArran)/cmm1toau
wecm(IArran)=Yjk(1,0,IArran)/cmm1toau
wexecm(IArran)=-Yjk(2,0,IArran)/cmm1toau
becm(IArran)=Yjk(0,1,IArran)/cmm1toau
decm(IArran)=Yjk(0,2,IArran)/cmm1toau
hecm(IArran)=Yjk(0,3,IArran)/cmm1toau
alfecm(IArran)=-Yjk(1,1,IArran)/cmm1toau
betecm(IArran)=Yjk(1,2,IArran)/cmm1toau
gamecm(IArran)=Yjk(2,1,IArran)/cmm1toau
ELSEIF(MaxVib(IArran)>=1)THEN
wzero(IArran)=Yjk(0,0,IArran)
we(IArran)=Yjk(1,0,IArran)
be(IArran)=Yjk(0,1,IArran)
de(IArran)=Yjk(0,2,IArran)
he(IArran)=Yjk(0,3,IArran)
alfe(IArran)=-Yjk(1,1,IArran)
bete(IArran)=Yjk(1,2,IArran)
!now in inverse centemeters
wzerocm(IArran)=Yjk(0,0,IArran)/cmm1toau
wecm(IArran)=Yjk(1,0,IArran)/cmm1toau
becm(IArran)=Yjk(0,1,IArran)/cmm1toau
decm(IArran)=Yjk(0,2,IArran)/cmm1toau
hecm(IArran)=Yjk(0,3,IArran)/cmm1toau
alfecm(IArran)=-Yjk(1,1,IArran)/cmm1toau
betecm(IArran)=Yjk(1,2,IArran)/cmm1toau
ELSEIF(MaxVib(IArran)>=0)THEN
wzero(IArran)=Yjk(0,0,IArran)
IF(Yjk(0,1,IArran)/=Zero)be(IArran)=Yjk(0,1,IArran)
IF(Yjk(0,2,IArran)/=Zero)de(IArran)=Yjk(0,2,IArran)
IF(Yjk(0,3,IArran)/=Zero)he(IArran)=Yjk(0,3,IArran)
!now in inverse centemeters
wzerocm(IArran)=Yjk(0,0,IArran)/cmm1toau
IF(Yjk(0,1,IArran)/=Zero)becm(IArran)=Yjk(0,1,IArran)/cmm1toau
IF(Yjk(0,2,IArran)/=Zero)decm(IArran)=Yjk(0,2,IArran)/cmm1toau
IF(Yjk(0,3,IArran)/=Zero)hecm(IArran)=Yjk(0,3,IArran)/cmm1toau
ENDIF


WRITE(YFitter_Unit,*)
WRITE(YFitter_Unit,*)"Spectroscopic constants for ", TRIM(DiatomicPot), "  IArran=", Iarran
WRITE(YFitter_Unit,'(A15,ES23.15,A17,ES23.15)')"wzero(IArran)=", wzero(IArran), "wzerocm(IArran)=", wzerocm(IArran)
WRITE(YFitter_Unit,'(A15,ES23.15,A17,ES23.15)')"we(IArran)=",    we(IArran),    "wecm(IArran)=",    wecm(IArran)
WRITE(YFitter_Unit,'(A15,ES23.15,A17,ES23.15)')"wexe(IArran)=",  wexe(IArran),  "wexecm(IArran)=",  wexecm(IArran)
WRITE(YFitter_Unit,'(A15,ES23.15,A17,ES23.15)')"weye(IArran)=",  weye(IArran),  "weyecm(IArran)=",  weyecm(IArran)
WRITE(YFitter_Unit,'(A15,ES23.15,A17,ES23.15)')"be(IArran)=",    be(IArran),    "becm(IArran)=",    becm(IArran)
WRITE(YFitter_Unit,'(A15,ES23.15,A17,ES23.15)')"de(IArran)=",    de(IArran),    "decm(IArran)=",    decm(IArran)
WRITE(YFitter_Unit,'(A15,ES23.15,A17,ES23.15)')"he(IArran)=",    he(IArran),    "hecm(IArran)=",    hecm(IArran)
WRITE(YFitter_Unit,'(A15,ES23.15,A17,ES23.15)')"alfe(IArran)=",  alfe(IArran),  "alfecm(IArran)=",  alfecm(IArran)
WRITE(YFitter_Unit,'(A15,ES23.15,A17,ES23.15)')"bete(IArran)=",  bete(IArran),  "betecm(IArran)=",  betecm(IArran)
WRITE(YFitter_Unit,'(A15,ES23.15,A17,ES23.15)')"game(IArran)=",  game(IArran),  "gamecm(IArran)=",  gamecm(IArran)


IF(IARRAN==3)THEN
   DEALLOCATE(Yfit)
   DEALLOCATE(Amat)
   DEALLOCATE(Work)
ENDIF

RETURN
ENDSUBROUTINE YFitter