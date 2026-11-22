SUBROUTINE aph3d1 (rmat, eye11, w, vecnew, tstore, vecold, xsq, xk, nidm)
!

!            P U R P O S E   O F    S U B R O U T I N E
! This is the Main program for the APH3D1 propagator.

!            I N P U T    A R G U M E N T S

!            O U T P U T    A R G U M E N T S

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
USE FileUnits_Module
USE InputFile_Module
USE TotJ_Module      ! PORT (2026-10-02): jtot for the URU record
USE nstate_Module
USE Narran_Module
USE parms_Module
USE Masses_Module
USE TotalEng_Module
USE Regins_Module
USE Title1_Module
USE region_MODULE
USE Storage_Module
USE convrsns_MODULE
USE aphchl_MODULE
USE VFunc_Module
USE NumNuj_Module
USE chltot_MODULE
USE engpro_MODULE
USE APHchl_Module, ONLY: Neigmin_Sum
IMPLICIT NONE
REAL(Kind=WP_Kind) version, eredev
REAL(Kind=WP_Kind) rmat(nstate,nstate), eye11(nstate,nstate), w(nstate,nstate)
REAL(Kind=WP_Kind) vecnew(nstate,nstate), tstore(nstate,nstate), vecold(nstate,nstate)
REAL(Kind=WP_Kind) xsq(nstate), xk(nstate)
INTEGER ithcall, ithsub, nidm, njacobi,i
  INTEGER :: irow, icol   ! PORT: loop indices for the Delves R-matrix dump
  INTEGER, ALLOCATABLE :: lam_n(:)
  INTEGER :: ierr_lam, lamfst, lamlst
  NAMELIST/basis_n/ lam_n, lamfst, lamlst
  ! PORT (2026-10-02): self-consistent Delves wiring -- the Delves channel
  ! space is defined by the U-matrix (BinOut/Umat.bin) written by the
  ! APH_to_Delves stage.  Read its header to obtain ndelve and the channel
  ! quantum numbers, propagate the classic Delves region in that space and
  ! write the R-matrix at rf as a 'URU_Matrix' record (the file Kmatrx_Main
  ! consumes; UMatph is skipped when the Delves region is active).
  INTEGER :: umat_unit2, naph_u, ndelve_u, ios_u
  INTEGER, ALLOCATABLE :: kchan3(:), elect3(:), mvib3(:), jrot3(:), lorb3(:)
  REAL(Kind=WP_Kind), ALLOCATABLE :: en3(:)
  REAL(Kind=WP_Kind) :: uru_r(nstate,nstate)
  CHARACTER(LEN=20) :: proj_u, approx_u, symm_u, type_u
  INTEGER, PARAMETER :: uru_unit = 85
  LOGICAL :: uru_written
CHARACTER(LEN=10) hour
CHARACTER(LEN=8) today
CHARACTER(LEN=5) curzone
INTEGER dtvalues(8)
LOGICAL iaph, idelve, ijacobi, little, medium, full
EXTERNAL popt, date, time, JacBasis, vsets, aphprop, delprop, jacprop
DATA iaph/.false./, idelve/.false./, ijacobi/.false./

!------------------------------------------------------------------
! Determine printing options.
!------------------------------------------------------------------
DATA little/.false./,medium/.false./,full/.false./, ithcall/0/,ithsub/0/,version/3.3/

CALL popt('aph3d1   ',little,medium,full,ithcall,ithsub)
!------------------------------------------------------------------
! Call system routines to determine the date and time.
!------------------------------------------------------------------
CALL dateh(today)
CALL timeh(hour)

!------------------------------------------------------------------
! WRITE the titles for this particular run.
!------------------------------------------------------------------
WRITE(Out_Unit,1131)
WRITE(Out_Unit,*)'These results were made with ', 'APH3D1 version ',version
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
IF(kenergy>1)THEN
    irdindep=.true.
    iwrindep=.false.
    OPEN(Unit=Basis_Bin_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Basis.bin',form='unformatted',status='old')
    REWIND frst_unit
ELSE
IF(iwrindep)THEN
    OPEN(Unit=Basis_Bin_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Basis.bin',form='unformatted', status='unknown')
WRITE(frst_unit)usys2
ELSEIF(irdindep)THEN
    OPEN(Unit=Basis_Bin_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Basis.bin',form='unformatted',status='old')
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
!      CALL JacBasis (nchanl, vecold, tstore, xsq, xk, finish, chinuj, nnuj, njacobi)
WRITE(Out_Unit,*)'finished jacobi basis'
!--------------------------------------------------------------------
! Set the initial R-matrix to essentially a zero matrix.
! The diagonal elements are set to a very small number to
! avoid overflow problems in inverting the R-matrix to
! obtain the Log-Derivative matrix.
!--------------------------------------------------------------------
! PORT (2026-10-10): see aph3d.f90 -- naph must follow $basis_n when the
! ABM surface-function stage is skipped (Neigmin_Sum left at zero).
IF (Neigmin_Sum == 0) THEN
   ALLOCATE(lam_n(0:2000)); lam_n = 0
   OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
   READ(In_Unit,NML=basis_n,IOSTAT=ierr_lam)
   CLOSE(Unit=In_Unit)
   IF(ierr_lam/=0) STOP 'aph3d1: missing Namelist Basis_n'
   Neigmin_Sum = SUM(lam_n)
   DEALLOCATE(lam_n)
END IF
naph=Neigmin_Sum
! PORT (2026-10-02): njacobi was never assigned in this routine (garbage ->
! automatic arrays with a garbage size -> SIGSEGV in overlap.f90:48).
! The Delves channel space is the one defined by the U-matrix.
njacobi = naph
uru_written = .FALSE.
IF(endaph < enddelve)THEN
   ALLOCATE(kchan3(nstate), elect3(nstate), mvib3(nstate), jrot3(nstate), lorb3(nstate), en3(nstate))
   OPEN(Unit=umat_unit2, File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Umat.bin', &
        Form='UNFORMATTED', Status='OLD', IOSTAT=ios_u)
   IF(ios_u == 0)THEN
      READ(umat_unit2, IOSTAT=ios_u) naph_u, ndelve_u
      IF(ios_u == 0 .AND. ndelve_u > 0 .AND. ndelve_u <= nstate)THEN
         READ(umat_unit2) (en3(i), i=1,ndelve_u)
         READ(umat_unit2) (kchan3(i), i=1,ndelve_u)
         READ(umat_unit2) (elect3(i), i=1,ndelve_u)
         READ(umat_unit2) (mvib3(i), i=1,ndelve_u)
         READ(umat_unit2) (jrot3(i), i=1,ndelve_u)
         READ(umat_unit2) (lorb3(i), i=1,ndelve_u)
         njacobi = ndelve_u
         WRITE(Out_Unit,'("aph3d1: njacobi set from Umat.bin = ",I6)') njacobi
      ELSE
         WRITE(Out_Unit,'("aph3d1: could not use Umat.bin (ios=",I4," ndelve=",I6,")")') ios_u, ndelve_u
      ENDIF
      CLOSE(umat_unit2)
   ELSE
      WRITE(Out_Unit,*)'aph3d1: Umat.bin not found; using njacobi=naph'
   ENDIF
ENDIF
WRITE(Out_Unit,'("aph3d1: naph=",I6," njacobi=",I6)') naph, njacobi
rmat=0.0d0
DO i=1,naph
   rmat(i,i)=1d-30
ENDDO

IF(naph>nstate)THEN
    Write(Out_Unit,*)'Error.  naph greater than nstate'
    STOP 'aph3d1'
ENDIF
!--------------------------------------------------------------------
! Propagate the coupled equations in the APH region from start
! to endaph if this region has nonzero lenght.
!--------------------------------------------------------------------
WRITE(Out_Unit,*)start,endaph
IF(start<=endaph)THEN
    iregion='aph    '
    iaph=.true.
    WRITE(Out_Unit,*)' A P H    P R O P A G A T I O N'
    CALL aphprop (naph,rmat,eye11,w,vecnew, tstore,vecold,xsq,start,endaph,njacobi,xk,nidm)
    WRITE(Out_Unit,*)'END APH PROPAGATION'
ENDIF
!--------------------------------------------------------------------
! Propagate the coupled equations in the Delves region from
! endaph to enddelve if this region has nonzero length.
!--------------------------------------------------------------------
WRITE(Out_Unit,*)endaph,enddelve
IF(endaph<enddelve)THEN
    iregion='delves '
    idelve=.true.
    WRITE(Out_Unit,*)' D E L V E S    P R O P A G A T I O N'
    CALL delprop (njacobi, rmat, eye11, w,  vecnew, tstore, vecold, xsq, endaph, enddelve)
    WRITE(Out_Unit,*)'END DELVE PROPAGATION'
    ! PORT (2026-10-02): dump the Delves-basis R-matrix at enddelve so the
    ! result of the Delves propagation is available as an artifact (the
    ! downstream asymptotic stage still consumes the APH-region R-matrix).
    OPEN(Unit=86, File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/R_Matrix_Delves.bin', &
         Form='UNFORMATTED', Status='UNKNOWN')
    WRITE(86) njacobi, enddelve, Etot
    WRITE(86) ((rmat(irow,icol), irow=1,njacobi), icol=1,njacobi)
    CLOSE(86)
    WRITE(Out_Unit,'("Wrote BinOut/R_Matrix_Delves.bin: njacobi=",I5," r=",ES12.4)') njacobi, enddelve
    ! PORT (2026-10-02): self-consistent wiring -- write the Delves-basis
    ! R-matrix at rf as a 'URU_Matrix' record (format of Write_Matrix2Old),
    ! which is what Kmatrx_Main reads.  The delprop propagation carries the
    ! space of the J*U-matrix, so no U-transformation is needed afterwards
    ! (BasisForScattering skips UMatph when endaph<enddelve).
    uru_block: BLOCK
       proj_u = 'delves              '
       approx_u = 'delves              '
       symm_u = 'none                '
       ! PORT (2026-10-03): TypeOfMat must be a padded 20-character variable.
       ! Passing the 10-character literal 'URU_Matrix' to the CHARACTER(LEN=20)
       ! dummy of Write_Matrix2Old (implicit interface) copies 20 bytes from the
       ! literal's storage, so the record header came out as
       ! 'URU_Matrix\0END JACOB' and Read_Matrix2 rejected the record.
       type_u = 'URU_Matrix          '
       uru_r = 0.d0
       DO icol = 1, njacobi
          DO irow = 1, njacobi
             uru_r(irow,icol) = rmat(irow,icol)
          ENDDO
       ENDDO
       IF(.NOT.uru_written)THEN
          OPEN(Unit=uru_unit, File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/URU_Matrix.bin', &
               Form='UNFORMATTED', Status='UNKNOWN')
          uru_written = .TRUE.
       ENDIF
       ! PORT (2026-10-03): the record index must be 1 (Read_Matrix2 requires
       ! jthmat == ithmat+1 with ithmat=0 on the first pass); 0 was written
       ! before and made Kmatrx_Main stop with 'STOP Read_Matrix2'.
       CALL Write_Matrix2Old(uru_r, en3, kchan3, elect3, mvib3, jrot3, lorb3, njacobi, &
                             Etot, 1, type_u, jtot, proj_u, approx_u, symm_u, uru_unit)
       WRITE(Out_Unit,'("Wrote BinOut/URU_Matrix.bin record: njacobi=",I5," Etot=",ES14.6)') njacobi, Etot
    END BLOCK uru_block
ENDIF
!--------------------------------------------------------------------
! Propagate the jacobi region from enddelve to finish if this
! region has nonzero length.
!--------------------------------------------------------------------
WRITE(Out_Unit,*)enddelve,finish
IF(enddelve<finish)THEN
    iregion='jacobi '
    ijacobi=.true.
    WRITE(Out_Unit,*)' J A C O B I    P R O P A G A T I O N'
    CALL jacprop (njacobi,rmat,eye11,w,vecnew, tstore,vecold,xsq,enddelve,finish)
    WRITE(*,*)'END JACOBI PROPAGATION'
    WRITE(Out_Unit,*)'END JACOBI PROPAGATION'
ENDIF
Write(*,*)'End APH3D1 (Matrix Elements)'
Write(Out_Unit,*)'End APH3D1 (Matrix Elements)'
Write(Out_Unit,*)'End APH3D1 (Matrix Elements)'
!--------------------------------------------------------------------
! End of the energy do loop.
!--------------------------------------------------------------------
CLOSE(Unit=Basis_Bin_Unit,status='keep')
! PORT (2026-10-03): flush the URU record written above.  Kmatrx_Main reads
! it back on a different unit; without this CLOSE a record smaller than the
! Fortran I/O buffer would still be unflushed and the reader would hit EOF.
IF(uru_written) CLOSE(Unit=uru_unit, status='keep')
RETURN

19 FORMAT(1x,'Time: ',a8)
18 FORMAT(1x,'Date: ',a8)
1131 FORMAT(///)
1132 FORMAT(1x,a70)
!----------------------***end-aph3d1----------------------------------
ENDSUBROUTINE aph3d1
