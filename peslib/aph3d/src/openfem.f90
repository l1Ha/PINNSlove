SUBROUTINE OpenFEM(jtot)
USE FileUnits_Module
USE Storage_Module
!
!   $RCSfile: openfem.f,v $ $Revision: 1.15 $
!   $Date: 89/11/16 11:03:16 $
!   $State: Stable $
!

!            P U R P O S E   O F    S U B R O U T I N E
!-----------------------------------------------------------------------
! This routine opens MOST files used by the code.  A few must be opened
! and closed in the code itself.  Those routines should have a
! commented out OPEN and some comments here.  A few more are in the
! program statement and are not described here.
!-----------------------------------------------------------------------

IMPLICIT NONE
INTEGER jtot

!-----------------------------------------------------------------------
! UNIT=Pmatrx_FEM_Unit Pmatrx_FEM is an unformatted file
!-----------------------------------------------------------------------
OPEN(Unit=Pmatrx_FEM_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/PotMatrx_FEM.bin',form='unformatted', status='unknown')

IF(jtot/=0)THEN
!-----------------------------------------------------------------------
! UNIT=Diag_Jtot_FEM_Unit diag_jtot_FEM is an unformatted file for storing the coupling
!        matrix elements.
!-----------------------------------------------------------------------
   OPEN(Unit=Diag_Jtot_FEM_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Diag_Jtot_FEM.bin',  form='unformatted', status='unknown')
!-----------------------------------------------------------------------
! UNIT=Diag_Lambda_FEM_Unit diag_lambda_FEM is an unformatted file for storing the coupling
!        matrix elements.
!-----------------------------------------------------------------------
   OPEN(Unit=Diag_Lambda_FEM_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Diag_Lambda_FEM.bin', form='unformatted', status='unknown')
!-----------------------------------------------------------------------
! UNIT=Asym_Top_FEM asym_top_FEM is an unformatted file for storing the coupling
!        matrix elements.
!-----------------------------------------------------------------------
   OPEN(Unit=Asym_Top_FEM_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Asym_Top_FEM.bin', form='unformatted',  status='unknown')
!-----------------------------------------------------------------------
! UNIT=Coriolis_FEM_Unit, coriolis_FEM is an unformatted file for storing the coupling
!        matrix elements.
!-----------------------------------------------------------------------
   OPEN(Unit=Coriolis_FEM_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Coriolis_FEM.bin', form='unformatted',  status='unknown')
ENDIF

!-----------------------------------------------------------------------
! UNIT=FEM_Ovrlp_Unit Ovrlp is the overlap matrix elements in the APH region. This
!         file is created in program ovr and used in the aph3d program.
!-----------------------------------------------------------------------
OPEN(Unit=FEM_Ovrlp_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Ovrlp_FEM.bin',form='unformatted', status='unknown')
!-----------------------------------------------------------------------
! UNIT=ElmDef_Unit elmdef is the unformatted file created in the sfunc program
!         that stores finite element locations. This file should
!         contain only the last hyperradius set of element definitions
!         when running APH3D.
!-----------------------------------------------------------------------
OPEN(Unit=ElmDef_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Elmdef_FEM.bin',form='unformatted', status='unknown')
!-----------------------------------------------------------------------
! UNIT=Nmodes_FEM_Bin_Unit nmodes is the unformatted file created in the sfunc program
!         that stores the surface functions at each hyperradius.
!          When running APH3d this file should contain only the
!          last hyperradius set of surface functions.
!-----------------------------------------------------------------------
OPEN(Unit=Nmodes_FEM_Bin_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Nmodes_FEM.bin',form='unformatted', status='unknown')
!-----------------------------------------------------------------------
! UNIT=eigfuns_Unit eigfuns is an unformatted file created in routine VibFun_Old or VibFun_New that
!         stores the eigenfunctions. This is a energy independent file
!         made on the first energy and needed on subsequent energies.
!-----------------------------------------------------------------------
!     OPEN(Unit=eigfuns_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'eigfuns',form='formatted',status='unknown')
!-------------------------------------------------------------------
! UNIT=Hamil_FEM_Unit hamil is an unformatted file created in sfunc that stores
!         the surface function hamiltonian matrix.
!-------------------------------------------------------------------
OPEN(Unit=Hamil_FEM_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Hamil_FEM.bin',form='unformatted', status='unknown')
!-------------------------------------------------------------------
! UNIT=SOvrlp_Unit sOvrlp is an unformatted file created in sfunc that stores
!         the surface function overlap matrix.
!-------------------------------------------------------------------
OPEN(Unit=SOvrlp_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/SOvrlp_FEM.bin',form='unformatted',status='unknown')
!-------------------------------------------------------------------
! UNIT=Temp1_FEM_Unit temp1 is an unformatted file created in sspace that stores
!         temporary DATA used in the subspace iteration method.
!-------------------------------------------------------------------
OPEN(Unit=Temp1_FEM_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Temp1_FEM.bin',form='unformatted', status='unknown')
!-------------------------------------------------------------------
! UNIT=Temp2_FEM_Unit temp2 is an unformatted file created in sspace that stores
!         temporary DATA used in the subspace iteration method.
!-------------------------------------------------------------------
OPEN(Unit=Temp2_FEM_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Temp2_FEM.bin',form='unformatted', status='unknown')
!-------------------------------------------------------------------
! UNIT=Relax_Unit relax is an unformatted file for storing overrelaxation
!         vectors in routine rapid.
!-------------------------------------------------------------------
OPEN(Unit=Relax_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Relax_FEM.bin',form='unformatted', status='unknown')
!-------------------------------------------------------------------
! UNIT=Tidi_Unit tidi is an unformatted file created in sfunc that stores
!         the vector idi temporarily.
!-------------------------------------------------------------------
OPEN(Unit=Tidi_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Tidi_FEM.bin',form='unformatted', status='unknown')
!-----------------------------------------------------------------------
! UNIT=sf33 Surface function DATA. Another rho.
!-----------------------------------------------------------------------
OPEN(Unit=Sf33_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Surf33_FEM.bin',form='unformatted', status='unknown')
!----------------------------------------------------------------------
! UNIT=sf34 Surface function DATA.
!----------------------------------------------------------------------
OPEN(Unit=sf34_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Surf34_FEM.bin',form='unformatted',  status='unknown')
!----------------------------------------------------------------------
! UNIT=Msher_Bin_Unit the nodedata is an unformatted file created by msher to store nodal DATA to be used by adini.
!----------------------------------------------------------------------
OPEN(Unit=Msher_Bin_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Nodedata_FEM.bin',form='unformatted',status='unknown')
!-------------------------------------------------------------------
! UNIT=Msglcz_Unit msglcz is a formatted file created by lanczo that contains
!         diagnostic and error messages from the lanczos method.
!-------------------------------------------------------------------
OPEN(Unit=Msglcz_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/Msglcz_FEM.txt',form='formatted',  status='unknown')
!-------------------------------------------------------------------
! UNIT=Intrp1_FEM_Unit intrp1 is a unformatted file created by interpol that contains
!         the map between the surface fns between two rhos.
!-------------------------------------------------------------------
OPEN(Unit=Intrp1_FEM_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Intrp1_FEM.bin',form='unformatted', status='unknown')
!-------------------------------------------------------------------
! UNIT=Intrp2_FEM_Unit intrp2 is a unformatted file created by interps that contains
!         the interpolated values of the surface function.
!-------------------------------------------------------------------
OPEN(Unit=Intrp2_FEM_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Intrp2_FEM.bin',form='unformatted', status='unknown')
!----------------------------------------------------------------------
! UNIT=TempNode_Bin_Unit tempnode is an unformatted file created by ncommon. This
!         is a temporary READ and should be passed directly.
!----------------------------------------------------------------------
OPEN(Unit=TempNode_Bin_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Tempnode_FEM.bin',form='unformatted',  status='unknown')
!-------------------------------------------------------------------
! UNIT=OldMsh_Unit oldmsh is an unformatted file created in mesher that stores
!         the old mesh.
!-------------------------------------------------------------------
OPEN(Unit=OldMsh_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Oldmsh_FEM.bin',form='unformatted', status='unknown')
!-------------------------------------------------------------------
! UNIT=Mapmsh_Bin_Unit mapmsh is an unformatted file created in mesher that stores
!         the map mesh.
!-------------------------------------------------------------------
OPEN(Unit=Mapmsh_Bin_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Mapmsh_FEM.bin',form='unformatted', status='unknown')
!-------------------------------------------------------------------
! UNIT=Zromsh_Bin_Unit zromsh is an unformatted file created in mesher that stores
!         the zero mesh.
!-------------------------------------------------------------------
OPEN(Unit=Zromsh_Bin_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Zromsh_FEM.bin',form='unformatted', status='unknown')
!-------------------------------------------------------------------
! UNIT=ZroMeg_Unit zromeg is an unformatted file created in wrmod that stores
!         the mega=megamin eigenvectors for use at next rho.
!        It is READ by previous.
!-------------------------------------------------------------------
OPEN(Unit=ZroMeg_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Zromeg_FEM.bin',form='unformatted',  status='unknown')
!-------------------------------------------------------------------
! UNIT=Nzromg_Unit nzromg is an unformatted file created in wrmod that stores
!     the most recent eigenvectors for use at next mega at same rho.
!  UNITS Lancos80_Unit and Lancos81_Unit are used instead of ZroMeg and Nzromg_Unit when Lanczos is used.
!-------------------------------------------------------------------
OPEN(Unit=Nzromg_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Nzromg_FEM.bin',form='unformatted', status='unknown')
!-----------------------------------------------------------------------
!  UNIT=PhiRms1_Unit phirms1 is an unformatted file used to store new phirms.
!  ---------------------------------------------------------------------
OPEN(Unit=PhiRms1_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Rmsphi1_FEM.bin',form='unformatted',status='unknown')
!  ---------------------------------------------------------------------
!  UNIT=PhiRms2_Unit phirms2 is an unformatted file used to store old phirms.
!  ---------------------------------------------------------------------
OPEN(Unit=PhiRms2_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Rmsphi2_FEM.bin',form='unformatted', status='unknown')
!  ---------------------------------------------------------------------
!  UNIT=Wrmod_Unit Sfelevl_FEM is an formatted file used to store the surface
!                  function energy levels in EV with vzero removed.
!  ---------------------------------------------------------------------
OPEN(Unit=Wrmod_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'GraphicsOut/Sfelevl_FEM.rbw',form='formatted',status='unknown')
RETURN
END
