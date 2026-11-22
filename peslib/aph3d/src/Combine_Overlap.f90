SUBROUTINE Combine_Overlap(CallFEM,CallDVR,CallPDAF,CallABM, NCurves)
USE Numeric_Kinds_Module
USE FileUnits_Module
USE Parms_Module, ONLY: MaxDim
IMPLICIT NONE
LOGICAL there, CallFEM,CallDVR,CallPDAF,CallABM
CHARACTER(LEN=20) infile1, infile3
CHARACTER(LEN=27)  infile2
!INTEGER, PARAMETER :: maxdim=600
INTEGER ithrho, jbasis, ibasis, ithFEM, neigmin, NCurves
REAL(Kind=WP_Kind) s(maxdim,maxdim), rholast, rho


IF(CallFEM)THEN
   infile1 = 'BinOut/Ovrlp_FEM.bin'
   infile2 = 'BinOut/Ovrlp_FEM_to_ABM.bin'
   infile3 = 'BinOut/Ovrlp_ABM.bin'
ELSEIF(CallDVR)THEN
   infile1 = 'BinOut/Ovrlp_DVR.bin'
   infile2 = 'BinOut/Ovrlp_DVR_to_ABM.bin'
   infile3 = 'BinOut/Ovrlp_ABM.bin'
ELSEIF(CallPDAF)THEN
   infile1 = 'BinOut/Ovrlp_PDAF.bin'
   infile2 = 'BinOut/Ovrlp_PDAF_to_ABM.bin'
   infile3 = 'BinOut/Ovrlp_ABM.bin'
ELSEIF(CallABM)THEN
   WRITE(Out_Unit,*)'Only ABM is used'
   Print *,'Only ABM is used'
   infile3 = 'BinOut/Ovrlp_ABM.bin'
ELSE
   WRITE(Out_unit,*)'ERROR: Invalid type'
   STOP 'Invalid type'
ENDIF

OPEN(Unit=Ovrlp_Bin_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Ovrlp.bin',form='unformatted',status='unknown')
ithrho=0
ithFEM=0
!-----------------------------------------------------------------------
! READ in DVR or FEM matrix elements.
!-----------------------------------------------------------------------
IF(CallFEM.or.CallDVR.or.CallPDAF)THEN
   WRITE(Out_unit,*)'Reading in DVR or FEM matrix elements'
   OPEN(Unit=DVR_OR_FEM_Ovrlp_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//infile1,form='unformatted',status='old')
10 CONTINUE
   READ(DVR_OR_FEM_Ovrlp_Unit) ithrho, rholast, rho
   READ(DVR_OR_FEM_Ovrlp_Unit,END=100) neigmin, rholast, rho
   IF(neigmin>maxdim)THEN
      WRITE(Out_unit,*) 'maxdim too small, change it in Combine_Overlap.f Loc1'
      WRITE(Out_unit,*) 'neigmin,maxdim',neigmin,maxdim
      WRITE(*,*) 'maxdim too small, change it in Combine_Overlap.f: Loc1'
      WRITE(*,*) 'neigmin,maxdim',neigmin,maxdim
      STOP 'Combine_Overlap'
   ENDIF
   WRITE(Ovrlp_Bin_Unit) ithrho, rholast, rho
   WRITE(Ovrlp_Bin_Unit) neigmin, rholast, rho
   WRITE(Out_unit,*) ithrho, rholast, rho, neigmin
   DO jbasis=1,neigmin
      READ(DVR_OR_FEM_Ovrlp_Unit)(s(ibasis,jbasis),ibasis=1,neigmin)
   ENDDO
   DO jbasis=1,neigmin
      WRITE(Ovrlp_Bin_Unit)(s(ibasis,jbasis),ibasis=1,neigmin)
   ENDDO
   goto 10
100 CONTINUE
   CLOSE(DVR_OR_FEM_Ovrlp_Unit)
   ithrho=ithrho-1
   ithFEM = ithrho

   !-----------------------------------------------------------------------
   ! READ in DVR_to_ABM or FEM_to_ABM matrix elements
   !-----------------------------------------------------------------------
   WRITE(Out_unit,*)'Reading in DVR_to_ABM or FEM_to_ABM matrix elements'
   OPEN(Unit=ABM_OR_FEM_TO_ABM_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//infile2,form='unformatted',status='old')
   READ(ABM_OR_FEM_TO_ABM_Unit) ithrho, rholast, rho
   READ(ABM_OR_FEM_TO_ABM_Unit) neigmin, rholast, rho
   IF(neigmin>maxdim)THEN
      WRITE(Out_unit,*) 'maxdim too small, change it in Combine_Overlap.f: Loc2'
      WRITE(Out_unit,*) 'neigmin,maxdim',neigmin,maxdim
      WRITE(*,*) 'maxdim too small, change it in Combine_Overlap.f Loc2'
      WRITE(*,*) 'neigmin,maxdim',neigmin,maxdim
      STOP 'Combine_Overlap'
   ENDIF
   WRITE(Ovrlp_Bin_Unit) ithFEM+1, rholast, rho
   WRITE(Ovrlp_Bin_Unit) neigmin, rholast, rho
   WRITE(Out_unit,*) ithFEM+1, rholast, rho, neigmin
   DO jbasis=1,neigmin
      READ(ABM_OR_FEM_TO_ABM_Unit)(s(ibasis,jbasis),ibasis=1,neigmin)
   ENDDO
   DO jbasis=1,neigmin
      WRITE(Ovrlp_Bin_Unit)(s(ibasis,jbasis),ibasis=1,neigmin)
   ENDDO
   CLOSE(ABM_OR_FEM_TO_ABM_Unit)
   ithFEM = ithFEM + 1
ENDIF

!-----------------------------------------------------------------------
! READ in ABM matrix elements.
!-----------------------------------------------------------------------
WRITE(Out_unit,*)'Reading in ABM matrix elements'
OPEN(Unit=PotMatrx_ABM_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//infile3,form='unformatted',status='old')
20    CONTINUE
READ(PotMatrx_ABM_Unit) ithrho, rholast, rho
READ(PotMatrx_ABM_Unit,END=200) neigmin, rholast, rho

IF(neigmin>maxdim)THEN
   WRITE(Out_unit,*) 'maxdim too small, change it in Combine_Overlap.f Loc3'
   WRITE(Out_unit,*) 'neigmin,maxdim',neigmin,maxdim
      WRITE(*,*) 'maxdim too small, change it in Combine_Overlap.f: Loc3'
      WRITE(*,*) 'neigmin,maxdim',neigmin,maxdim
   STOP
ENDIF
WRITE(Ovrlp_Bin_Unit) ithrho+ithFEM, rholast, rho
WRITE(Ovrlp_Bin_Unit) neigmin, rholast, rho
WRITE(Out_unit,*) ithrho+ithFEM, rholast, rho, neigmin
DO jbasis=1,neigmin
   READ(PotMatrx_ABM_Unit)(s(ibasis,jbasis),ibasis=1,neigmin)
ENDDO
DO jbasis=1,neigmin
   WRITE(Ovrlp_Bin_Unit)(s(ibasis,jbasis),ibasis=1,neigmin)
ENDDO
GOTO 20
200   CONTINUE
WRITE(Ovrlp_Bin_Unit) ithrho+ithFEM, rholast, rho
CLOSE(PotMatrx_ABM_Unit)

CLOSE(Ovrlp_Bin_Unit)
NCurves=NeigMin

RETURN
ENDSUBROUTINE Combine_Overlap
