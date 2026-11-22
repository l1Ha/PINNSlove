SUBROUTINE aph3d (w, vecnew, tstore, vecold, xsq, xk, ndimen)
!
!
!            P U R P O S E   O F    S U B R O U T I N E
!
!            I N P U T    A R G U M E N T S
!
!            O U T P U T    A R G U M E N T S
!
!
!this routine is called by:
!         aph_to_delves
!this routine calls
!         popt,aphget,upsiln,aphdelbf
!-----------------------------------------------------------------------

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>


USE fileunits_Module
USE Numeric_Kinds_Module
USE Narran_Module
USE TotalEng_Module
USE Title1_Module
USE region_Module
USE Regins_Module
USE Region_Module
USE Oops_Module
USE Storage_Module
USE convrsns_Module
USE aphchl_Module
USE VFunc_Module
USE NumNuj_Module
USE chltot_Module
USE engpro_Module
!USE Gaussb_Module
USE GaussQuady_Module
USE qnumsp_Module
USE cherm1_Module
USE cherm2_Module
USE Parms_Module
USE Nidm_Module
USE APHchl_Module, ONLY: Neigmin_Sum
USE InputFile_Module

IMPLICIT NONE
CHARACTER(LEN=10) hour
CHARACTER(LEN=8) today
CHARACTER(LEN=5) curzone
LOGICAL INext, IAPH
LOGICAL  little, medium, full
INTEGER dtvalues(8)
INTEGER ithcall, ithsub, ndimen, njacobi
REAL(Kind=WP_Kind) version, eredev
REAL(Kind=WP_Kind) w(ndimen**2)
REAL(Kind=WP_Kind) vecnew(ndimen**2)
REAL(Kind=WP_Kind) tstore(ndimen**2)
REAL(Kind=WP_Kind) vecold(ndimen**2)
REAL(Kind=WP_Kind) xsq(ndimen)
REAL(Kind=WP_Kind) xk(ndimen)
! PORT (2026-10-10): $basis_n is read here so that naph equals the requested
! propagation space even when the ABM surface-function stage is skipped.
INTEGER, PARAMETER :: maxl_nm = 2000
INTEGER, ALLOCATABLE :: lam_n(:)
INTEGER lamfst, lamlst
INTEGER :: ierr_lam
NAMELIST/basis_n/ lam_n, lamfst, lamlst
EXTERNAL  popt, aphget,upsiln,aphdelbf
DATA iaph/.false./
!------------------------------------------------------------------
! Determine printing options.
!------------------------------------------------------------------
DATA little/.false./,medium/.false./,full/.false./, ithcall/0/,ithsub/0/,version/3.3d0/

CALL popt('aph3d   ',little,medium,full,ithcall,ithsub)
!-------------------------------------------------------------------------
!ALLOCATE memory
!-------------------------------------------------------------------------
!ALLOCATE(xpth(maxhermt,narran),wpth(maxhermt,narran))
IF(.not.ALLOCATED(Numax))ALLOCATE(numax(nvibrot), numin(nvibrot))
IF(.NOT.ALLOCATED(chinuj))ALLOCATE(chinuj(nvbrthrt))
IF(.not.ALLOCATED(chnow))THEN
  ALLOCATE(chlast(maxosc*maxosc,0:mxl,narran))
  ALLOCATE(chnow(maxosc*maxosc,0:mxl,narran))
ENDIF
!
!------------------------------------------------------------------
! CALL system routines to determine the date and time.
!------------------------------------------------------------------
CALL date_and_time(today, hour, curzone, dtvalues)
!------------------------------------------------------------------
! WRITE the titles for this particular run.
!------------------------------------------------------------------
WRITE(Out_Unit,1131)
WRITE(Out_Unit,*)'These results were made with ', 'APH3D version ',version
!-------------------------------------------------------------
!   WRITE the time and date to unit 6.
!-------------------------------------------------------------
CALL Print_Date_And_Time(Today, Hour, CurZone, DtValues, Out_Unit)
WRITE(Out_Unit,1132)Title1
WRITE(Out_Unit,1132)Title2
WRITE(Out_Unit,1132)Title3
!------------------------------------------------------------------
! The maximum energy will determine the wavelength
! of the most oscillatory state which will be used to set the
! step size in the log-derivative propagation.
!------------------------------------------------------------------
engmax=emin+(numengs-1)*ejump
!------------------------------------------------------------------
! If irdindep=true the energy independent information has already
!             been stored on the disk and hence iwrindep=false.
! If irdindep=false we usually WRITE the energy independent stuff
!             to disk, but if only doing one energy both irdindep
!             and iwrindep can be false.
!------------------------------------------------------------------
IF(irdindep) iwrindep=.false.

!                   *****************************
!                   *   loop for each energy.   *
!                   *****************************
kenergy = 1
!------------------------------------------------------------------
! When irdindep=true the energy independent information will be
! read from disk.  This implies that a previous run has been made with
! iwrindep=true, irdindep=false.
!------------------------------------------------------------------
WRITE(*,*)"iregion=", iregion,"  ioops=", ioops, "Loc1"
irdindep = .false.
iwrindep = .false.
IF(kenergy>1)THEN
   irdindep=.true.
   iwrindep=.false.
   OPEN(Unit=bas_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Del_Basis.bin',form='unformatted', status='old')
   REWIND bas_unit
ELSE
  IF(iwrindep)THEN
    OPEN(Unit=bas_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Del_Basis.bin',form='unformatted',status='unknown')
  ELSEIF(irdindep)THEN
    OPEN(Unit=bas_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Del_Basis.bin', form='unformatted',status='old')
  ENDIF
ENDIF
!------------------------------------------------------------------
! Determine the Total collision energy.
!------------------------------------------------------------------
Etot=emin+(kenergy-1)*ejump
eredev=Etot*autoev
WRITE(Out_Unit,184)Etot,eredev
  184 FORMAT(//,'energy (Ha) =',f10.5,5x,'energy (ev) =',f10.5)
!--------------------------------------------------------------------
! Obtain the asymptotic Jacobi basis functions.
!--------------------------------------------------------------------
WRITE(Out_Unit,*)' G E N E R A T E    J A C O B I    B A S I S'
CALL JacBasis (nchanl, vecold, tstore, xsq, xk, finish, chinuj, nnuj, njacobi)
!--------------------------------------------------------------------
! Propagate the coupled equations in the APH region from start
! to endaph if this region has nonzero lenght.
!--------------------------------------------------------------------
iregion='aph    '
iaph=.true.
!--------------------------------------------------------------------
! Get the quantum numbers associated with each channel.
!--------------------------------------------------------------------
inext=.false.
! PORT (2026-10-10): with CalcSFunc=.False. the ABM_APH stage is skipped, so
! Neigmin_Sum is never accumulated and the APH channel count must come from
! the $basis_n namelist (previously hard-wired to 100, which silently
! truncated or over-ran the propagation space whenever lam_n/=100).
IF (.NOT. ALLOCATED(lam_n)) THEN
   ALLOCATE(lam_n(0:maxl_nm))
   lam_n = 0
   OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
   READ(In_Unit,NML=basis_n,IOSTAT=ierr_lam)
   CLOSE(Unit=In_Unit)
   IF(ierr_lam/=0) STOP 'aph3d: missing Namelist Basis_n'
END IF
IF (Neigmin_Sum == 0) Neigmin_Sum = SUM(lam_n)
naph=Neigmin_Sum
CALL aphget(naph)
!---------------------------------------------------------------------
! Calculate the Delves basis functions so that one can transform
! to Delves coordinates.
!---------------------------------------------------------------------
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)' G E N E R A T E    D E L V E S    B A S I S'
CALL upsiln (vecold, tstore, xsq, xk, endaph, chinuj)
!----------------------------------------------------------------------
! Now determine the U matrix that transforms the R matrix from the
! aph to the delves basis.
!----------------------------------------------------------------------
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'A P H    T O    D E L V E S    ', 'P R O J E C T I O N'
CALL aphdelbf (endaph, naph, njacobi, tstore, vecnew, w, nidm)
CLOSE(Unit=bas_unit,status='keep')
!-------------------------------------------------------------------------
!DEALLOCATE memory
!-------------------------------------------------------------------------
!DEALLOCATE(xpth,wpth)
!DEALLOCATE(numax, numin)
!DEALLOCATE(chinuj)
!DEALLOCATE(chlast)
!DEALLOCATE(chnow)
!
WRITE(Out_unit,*)"APH3D Completed"
RETURN

 1131 FORMAT(//)
 1132 FORMAT(1x,a70)
!----------------------***END-aph3d----------------------------------
ENDSUBROUTINE aph3d
