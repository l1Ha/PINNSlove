SUBROUTINE Jacobi_Basis
USE Numeric_Kinds_Module
USE FIleUnits_Module
USE Narran_Module
USE Parms_Module
USE TotalEng_Module
USE engpro_Module
USE NumNuj_Module
USE Regins_Module
USE chltot_Module
USE VFunc_Module
USE region_Module
USE Storage_Module
USE GaussQuady_Module
USE qnumsp_Module
USE CHerm1_Module
USE CHerm2_Module
USE Oops_Module
IMPLICIT NONE
LOGICAL :: ijacobi
CHARACTER(LEN=10) hour
CHARACTER(LEN=8) today
CHARACTER(LEN=5) curzone
INTEGER dtvalues(8)
INTEGER :: njacobi
REAL(Kind=WP_Kind), ALLOCATABLE :: tstore (:,:)
REAL(Kind=WP_Kind), ALLOCATABLE :: vecold (:,:)
REAL(Kind=WP_Kind), ALLOCATABLE :: xsq (:)
REAL(Kind=WP_Kind), ALLOCATABLE :: xk (:)
EXTERNAL title, constant, readaph, JacBasis

!----------------------------------------------------------------------
! Get primary input DATA.
!----------------------------------------------------------------------

CALL ReadAllData
WRITE(Out_Unit,*)'Called    Jacobi_Basis'
CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'DiatomicOut/Jacobi_Basis.txt',Form='FORMATTED')
irdindep = .false.
iwrindep = .true.
WRITE(*,*)"iregion=", iregion,"  ioops=", ioops, "Loc4"
OPEN(Unit=bas_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Jac_Basis.bin', form = 'unformatted', status = 'unknown')
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
IF(.not.ALLOCATED(xpth))ALLOCATE(xpth(maxhermt,narran),wpth(maxhermt,narran))
!--------------------------------------------------------------------
! Obtain the asymptotic Jacobi basis functions.
!--------------------------------------------------------------------
kenergy = 1
Etot = emin + (kenergy - 1) * ejump
iregion = 'jacobi '
ijacobi = .true.
ioops=.false.
WRITE(Out_Unit,  * ) ' G E N E R A T E    J A C O B I    B A S I S'
IF(.not.allocated(chinuj))ALLOCATE(chinuj(nvbrthrt))
IF(.not.ALLOCATED(Numax))ALLOCATE(numax(nvibrot),numin(nvibrot))
IF(.not.ALLOCATED(chnow))THEN
   ALLOCATE(chlast(maxosc*maxosc,0:mxl,narran))
   ALLOCATE(chnow(maxosc*maxosc,0:mxl,narran))
ENDIF
ALLOCATE(tstore(nvibrot,nvibrot*150))
ALLOCATE(vecold(nvibrot,nvibrot*150))
ALLOCATE(xsq(nvibrot*13))
ALLOCATE(xk(nvibrot*13))
CALL JacBasis (nchanl, vecold, tstore, xsq, xk, finish, chinuj, nnuj, njacobi)

WRITE(Out_Unit,*)'Leaving Jacobi_Basis'
CALL Date_And_Time(today, hour, curzone, dtvalues)
CALL Print_Date_And_Time(Today, Hour, CurZone, DtValues, Out_Unit)
!DEALLOCATE(chinuj)
!DEALLOCATE(chlast)
!DEALLOCATE(chnow)
!DEALLOCATE(numax, numin)
DEALLOCATE(tstore)
DEALLOCATE(vecold)
DEALLOCATE(xsq)
DEALLOCATE(xk)
!DEALLOCATE(xpth,wpth)

CLOSE(Out_Unit)
CLOSE(bas_unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/BasisForScatt.txt',Form='FORMATTED',ACCESS='APPEND')
WRITE(Out_Unit,*)'Leaving Jacobi_Basis'
WRITE(Out_Unit,*)

RETURN
ENDSUBROUTINE jacobi_basis
