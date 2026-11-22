SUBROUTINE Delves_Basis
USE FileUnits_Module
USE Narran_Module
USE Parms_Module
USE NumNuj_Module
USE Regins_Module
USE chltot_Module
USE VFunc_Module
USE region_Module
USE Storage_Module
USE RhoValue_Module
USE GaussQuady_Module
USE qnumsp_Module
USE CHerm1_Module
USE CHerm2_Module
USE InputFile_Module
USE Oops_Module
USE Energy_Module
USE EngPro_module
USE TotalEng_Module
IMPLICIT NONE
LOGICAL :: ijacobi
CHARACTER(LEN=10) hour
CHARACTER(LEN=8) today
CHARACTER(LEN=5) curzone
INTEGER dtvalues(8)
INTEGER :: njacobi, IERR
REAL(Kind=WP_Kind), ALLOCATABLE :: tstore (:,:)
REAL(Kind=WP_Kind), ALLOCATABLE :: vecold (:,:)
REAL(Kind=WP_Kind), ALLOCATABLE :: xsq (:)
REAL(Kind=WP_Kind), ALLOCATABLE :: xk (:)
EXTERNAL title, constant, readaph, JacBasis

!----------------------------------------------------------------------
! Get primary input DATA.
!----------------------------------------------------------------------
irdindep = .false.
iwrindep = .true.
WRITE(Out_Unit,*)'Called    Delves_Basis'
CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'DiatomicOut/Jacobi_Basis.txt',Form='FORMATTED')
OPEN(Unit=bas_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Jac_Basis.bin', form = 'unformatted', status = 'unknown')
WRITE(Out_Unit,*)'Entering Delves_Basis'
!-------------------------------------------------------------
! Determine the time and todays date.
!-------------------------------------------------------------
CALL Date_And_Time(today, hour, curzone, dtvalues)
!-------------------------------------------------------------
!   WRITE the time and date to unit Out_Unit.
!-------------------------------------------------------------
CALL Print_Date_And_Time(Today, Hour, CurZone, DtValues, Out_Unit)

!-----------------------------------------------------------------------
! Establish all conversion factors and physical constants
!-----------------------------------------------------------------------

CALL constant
!-----------------------------------------------------------------------
! READ in all input DATA
!-----------------------------------------------------------------------

CALL readaph
IF(.NOT.allocated(xpth))ALLOCATE(xpth(maxhermt,narran),wpth(maxhermt,narran))
!--------------------------------------------------------------------
! Obtain the asymptotic Jacobi basis functions.
!--------------------------------------------------------------------
kenergy = 1
Emin=efirst
ejump=deltaeng
Etot = efirst + (kenergy - 1) * deltaeng
iregion = 'jacobi '
ijacobi = .true.
ioops=.false.
WRITE(Out_Unit,  * ) ' G E N E R A T E    J A C O B I    B A S I S'
WRITE(*,*)"iregion=", iregion,"  ioops=", ioops, "Loc3"
IF(.NOT.ALLOCATED(numax))ALLOCATE(numax(nvibrot+500),numin(nvibrot+500))
IF(.NOT.ALLOCATED(chlast))ALLOCATE(chlast(maxosc*maxosc,0:mxl,narran))
IF(.NOT.ALLOCATED(chnow))ALLOCATE(chnow(maxosc*maxosc,0:mxl,narran))
IF(.NOT.ALLOCATED(chinuj))ALLOCATE(chinuj(nvbrthrt))
IF(.NOT.ALLOCATED(tstore))ALLOCATE(tstore(nvibrot,nvibrot*150))
IF(.NOT.ALLOCATED(vecold))ALLOCATE(vecold(nvibrot,nvibrot*150))
IF(.NOT.ALLOCATED(xsq))ALLOCATE(xsq(nvibrot*13+500))
IF(.NOT.ALLOCATED(xk))ALLOCATE(xk(nvibrot*13+500))

OPEN(Unit=Oliver_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/oliver.txt',status='unknown')
CALL JacBasis (nchanl, vecold, tstore, xsq, xk, finish, chinuj, nnuj, njacobi)
CLOSE(Out_Unit)
CLOSE(bas_unit)

WRITE(Out_Unit,*)
WRITE(Out_Unit,  * ) ' G E N E R A T E    D E L V E S    B A S I S'
iregion = 'delves '
ijacobi = .false.
ioops=.true.
WRITE(Out_Unit,  * ) ' D E L V E S    B A S I S'
WRITE(*,*)"iregion=", iregion,"  ioops=", ioops, "Loc3"
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'DiatomicOut/Delves_Basis.txt',Form='FORMATTED')
OPEN(Unit=bas_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Del_Basis.bin', form = 'unformatted', status = 'unknown')
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile, status = 'old')
READ(In_Unit,NML=rhovalue,IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'RhoValue')
      RhoVal=Finish
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist RhoValues'
      WRITE(Msg_Unit,*)'ERROR with Namelist RhoValues'
      STOP 'ReadABM: RhoValues'
   ENDIF
CLOSE(unit=In_Unit)
WRITE(Out_Unit,NML=rhovalue)
iregion = 'delves'
ioops=.true.
CALL upsiln (vecold, tstore, xsq, xk, rhoval, chinuj)

WRITE(Out_Unit,*)'Leaving Delves_Basis'
CALL Date_And_Time(today, hour, curzone, dtvalues)
CALL Print_Date_And_Time(Today, Hour, CurZone, DtValues, Out_Unit)
!DEALLOCATE(chinuj)
!DEALLOCATE(chlast)
!DEALLOCATE(chnow)
!IF(Allocated(Numax))DEALLOCATE(numax, numin)
DEALLOCATE(tstore)
DEALLOCATE(vecold)
DEALLOCATE(xsq)
DEALLOCATE(xk)
!DEALLOCATE(xpth,wpth)
CLOSE(SBFrst_Unit)


CLOSE(Out_Unit)
CLOSE(bas_unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/BasisForScatt.txt',Form='FORMATTED',ACCESS='APPEND')
WRITE(Out_Unit,*)'Leaving Delves_Basis'
WRITE(Out_Unit,*)

RETURN
ENDSUBROUTINE Delves_Basis
