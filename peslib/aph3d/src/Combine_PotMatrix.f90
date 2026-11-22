SUBROUTINE Combine_PotMatrix(CallFEM,CallDVR,CallPDAF,CallABM)
USE Numeric_Kinds_Module
USE FileUnits_Module
USE Boundary_Module
IMPLICIT NONE
LOGICAL mid_zero, there, CallFEM,CallDVR,CallPDAF,CallABM
CHARACTER(LEN=24) infile
INTEGER krho, neigmin, nrho, ibasis, jbasis, ir, IthRho
REAL(Kind=WP_Kind) rhocent, eigen, s, rhos
INTEGER, PARAMETER :: maxdim=1000   ! PORT 2026-10-04: raised 600->1000
DIMENSION rhos(maxdim),eigen(maxdim),s(maxdim,maxdim)

OPEN(Unit=PotMatrix_Bin_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/PotMatrx.bin',form='unformatted',status='unknown')
IF(CallDVR)THEN
   infile = 'BinOut/PotMatrx_DVR.bin'
ELSEIF(CAllFEM)THEN
   infile = 'BinOut/PotMatrx_FEM.bin'
ELSEIF(CallPDAF)THEN
   infile = 'BinOut/PotMatrx_PDAF.bin'
ENDIF

!-----------------------------------------------------------------------
! READ potential matrix elements calculated using the FEM method.
!-----------------------------------------------------------------------
IF(CallDVR.or.CallFEM.or.CallPDAF)THEN
   OPEN(Unit=PotMatrx_FEM_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//TRIM(infile),form='unformatted',status='old')
   WRITE(Out_unit,*)'Starting to READ FEM or DVR matrix elements.'
   DO IthRho=IthStart_DVR,IthEnd_DVR
      WRITE(Msg_Unit,*)"IthRho=",IthRho," IthStart_DVR=",IthStart_DVR," IthEnd_DVR=",IthEnd_DVR
      WRITE(Out_Unit,*)"IthRho=",IthRho," IthStart_DVR=",IthStart_DVR," IthEnd_DVR=",IthEnd_DVR
      READ(PotMatrx_FEM_Unit,END=100) rhocent, neigmin, nrho, (rhos(krho),krho=1,nrho),mid_zero
      READ(PotMatrx_FEM_Unit) (eigen(ibasis), ibasis=1, neigmin)
      WRITE(Out_Unit,*)"rhocent=",rhocent," neigmin=",neigmin," nrho=",nrho
      WRITE(Out_Unit,'(5ES15.7)')(rhos(krho),krho=1,nrho)
      WRITE(Out_Unit,'(5ES15.7)')(eigen(ibasis), ibasis=1, neigmin)
      WRITE(Msg_Unit,*)"rhocent=",rhocent," neigmin=",neigmin," nrho=",nrho
      WRITE(Msg_Unit,'(5ES15.7)')(rhos(krho),krho=1,nrho)
      WRITE(Msg_Unit,'(5ES15.7)')(eigen(ibasis), ibasis=1,MIN(neigmin,20))
      WRITE(PotMatrix_Bin_Unit) rhocent, neigmin, nrho, (rhos(krho),krho=1,nrho), mid_zero
      WRITE(PotMatrix_Bin_Unit) (eigen(ibasis), ibasis=1, neigmin)
      WRITE(Out_unit,*) rhocent, neigmin, nrho, (rhos(krho),krho=1,nrho)
      DO ir=1,3,2
         DO jbasis=1,neigmin
            READ(PotMatrx_FEM_Unit)(s(ibasis,jbasis),ibasis=jbasis,neigmin)
         ENDDO
         DO jbasis=1,neigmin
            WRITE(PotMatrix_Bin_Unit)(s(ibasis,jbasis),ibasis=jbasis,neigmin)
         ENDDO
      ENDDO
   ENDDO
100 CONTINUE
   CLOSE(PotMatrx_FEM_Unit)
ENDIF

!-----------------------------------------------------------------------
! READ potential matrix elements calculated using the ABM method.
!-----------------------------------------------------------------------
OPEN(Unit=PotMatrx_ABM_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/PotMatrx_ABM.bin',form='unformatted',status='old')
WRITE(Out_unit,*)'Starting to READ ABM matrix elements.'
DO IthRho=IthStart_ABM,IthEnd_ABM
   WRITE(Msg_Unit,*)"IthRho=",IthRho," IthStart_ABM=",IthStart_ABM," IthEnd_ABM=",IthEnd_ABM
   WRITE(Out_Unit,*)"IthRho=",IthRho," IthStart_ABM=",IthStart_ABM," IthEnd_ABM=",IthEnd_ABM
   READ(PotMatrx_ABM_Unit,END=200) rhocent, neigmin, nrho, (rhos(krho),krho=1,nrho), mid_zero
   READ(PotMatrx_ABM_Unit) (eigen(ibasis), ibasis=1, neigmin)
   WRITE(Out_Unit,*)"rhocent=",rhocent," neigmin=",neigmin," nrho=",nrho
   WRITE(Out_Unit,'(5ES15.7)')(rhos(krho),krho=1,nrho)
   WRITE(Out_Unit,'(5ES15.7)')(eigen(ibasis), ibasis=1, neigmin)
   WRITE(Msg_Unit,*)"rhocent=",rhocent," neigmin=",neigmin," nrho=",nrho
   WRITE(Msg_Unit,'(5ES15.7)')(rhos(krho),krho=1,nrho)
   WRITE(Msg_Unit,'(5ES15.7)')(eigen(ibasis), ibasis=1,MIN(neigmin,20))
   WRITE(PotMatrix_Bin_Unit) rhocent, neigmin, nrho, (rhos(krho),krho=1,nrho), mid_zero
   WRITE(PotMatrix_Bin_Unit) (eigen(ibasis), ibasis=1, neigmin)
   WRITE(Out_unit,'(ES15.7,2I5,8ES15.7)') rhocent, neigmin, nrho, (rhos(krho),krho=1,nrho)
   DO ir=1,3,2
      DO jbasis=1,neigmin
         READ(PotMatrx_ABM_Unit)(s(ibasis,jbasis),ibasis=jbasis,neigmin)
         WRITE(PotMatrix_Bin_Unit)(s(ibasis,jbasis),ibasis=jbasis,neigmin)
      ENDDO
   ENDDO
ENDDO
200   CONTINUE
CLOSE(PotMatrx_ABM_Unit)

CLOSE(PotMatrix_Bin_Unit)
RETURN
ENDSUBROUTINE Combine_PotMatrix
