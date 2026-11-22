SUBROUTINE open_matph (pmatfile, ovrfile)
!-----------------------------------------------------------------------
! This routine opens MOST files used by the code.  A few must be opened
! and closed in the code itself.  Those routines should have a
! commented out open and some comments here.  A few more are in the
! program statement and are not described here.
!-----------------------------------------------------------------------
USE FileUnits_Module
USE Storage_Module
USE dip_Module
USE DoEnlvls_Module
IMPLICIT NONE
CHARACTER(LEN=*) pmatfile, ovrfile
!#include <Inc.unitno>
!#include <Inc.storage>
!#include <Inc.dip>
!
!-----------------------------------------------------------------------
! UNIT=PotMatrx_Unit Pmatrx is an unformatted file
!-----------------------------------------------------------------------
WRITE(Out_Unit,*)'PotMatrx_Unit pmatfile=',pmatfile
WRITE(Out_Unit,*)'unit7 ovrfile =',ovrfile
OPEN(Unit=PotMatrx_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//pmatfile,form='unformatted',status='old')
!-----------------------------------------------------------------------
! UNIT=frst_unit=firsteng is an unformatted filed created in the aph3d
!         programs which stores all of the energy independent
!         information so that subsequent runs at different energies
!         will execute quickly.
!-----------------------------------------------------------------------
IF(iwrindep)THEN
   OPEN(Unit=frst_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Firsteng.bin', form='unformatted',status='unknown')
ELSEIF(irdindep)THEN
   OPEN(Unit=frst_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Firsteng.bin',form='unformatted',status='old')
ENDIF
!-----------------------------------------------------------------------
! UNIT=Ovr_Unit=Ovrlp is the overlap matrix elements in the APH region.
!         This file is created in program ovr and used in the
!         aph3d program.
!-----------------------------------------------------------------------
OPEN(Unit=Ovr_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//ovrfile,form='unformatted',status='old')
!-----------------------------------------------------------------------
! UNIT=elvl_unit formatted file that contains the energy correlation
!            data in a form suitable for plotting
! created by enlvls
!            in the propagator.
!-----------------------------------------------------------------------
IF(lenlvls)THEN
   OPEN(Unit=elvl_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'elevels',form='formatted',status='unknown')
ENDIF
!--------------------------------------------------------------------
! Unit=tmutp_unit=tmutp is unformatted file created by DVR_dipole
!         which contains the dipole moment between surface functions
! Unit inho_unit=inho contains is an unformatted
!          file created by numbmc for the
!          lower state
! Unit fph_unit=firstengph is an unformatted file containing the energy
!         independent inhomogeneous term in photodissociation
!---------------------------------------------------------------------
IF(lphoto)THEN
    OPEN(Unit=tmutp_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Tmutp.bin', form='unformatted',status='old')
    OPEN(Unit=inho_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Inho.bin',form='unformatted',status='old')
    REWIND tmutp_unit
    REWIND inho_unit
    OPEN(Unit=fph_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Firstengph.bin',form='unformatted',status= 'unknown')
ENDIF
!

RETURN
ENDSUBROUTINE open_matph
