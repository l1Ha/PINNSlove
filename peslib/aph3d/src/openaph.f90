SUBROUTINE openaph
!
!-----------------------------------------------------------------------
!This routine opens MOST files used by the code.  A few must be opened
! and closed in the code itself.  Those routines should have a
! commented out open and some comments here.  A few more are in the
! program statement and are not described here.
!-----------------------------------------------------------------------
! TRP  02/18/1998  Added code to read filename vec_fname from file .vecunit
!-----------------------------------------------------------------------
!
!this routine is called by:
!         aph_to_delves
!this routine calls
!
!-----------------------------------------------------------------------


USE FileUnits_Module
USE Storage_Module

IMPLICIT NONE
!CHARACTER(LEN=100) vec_fname
!NAMELIST/vecunit/vec_fname
!-----------------------------------------------------------------------
! UNIT=progre_unit=progre is a formatted file
!-----------------------------------------------------------------------
OPEN(Unit=progre_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/progre.txt',   form='formatted', status='unknown')
!-----------------------------------------------------------------------
! UNIT=Plmst_Unit is an unformatted file used to store the Legendre polynomials
!         that are calculated in routine AphDel_Old or AphDel_New.
!-----------------------------------------------------------------------
!     OPEN(UnitPlmst_Unit ,File=OutDIR(1:LEN(TRIM(OutDIR)))//'plmst',    form='unformatted',
!    >      status='unknown')
!-----------------------------------------------------------------------
! UNIT=umat_unit=umat is an unformatted file containing the
!         u-matrix which
!         is used to transform the R-matrix from APH to DELVES labels.
!-----------------------------------------------------------------------
!     OPEN(Unit=umat_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'umat',     form='unformatted',
!     >      status='unknown')
!-----------------------------------------------------------------------
! UNIT=bas_unit=basis is an unformatted file that contains basis functions.
! A energy independent file made during the first energy by the
! propagator and needed on subsequent energies.
!-----------------------------------------------------------------------
!     OPEN(Unit=bas_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/basis', form='unformatted', status='unknown')
!-----------------------------------------------------------------------
! UNIT=elm_unit=elmdef is the unformatted file
! created in the sfunc program
!         that stores finite element locations. This file should
!         contain only the last hyperradius set of element definitions
!        when running APH3D.
!-----------------------------------------------------------------------
OPEN(Unit=elm_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Elmdef_FEM.bin',   form='unformatted')
!-----------------------------------------------------------------------
! UNIT=nmode_unit=nmodes is the unformatted file
!          created in the sfunc program
!          that stores the surface functions at each hyperradius.
!          When running APH3d this file should contain only the
!          last hyperradius set of surface functions.
!-----------------------------------------------------------------------
OPEN(Unit=nmode_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Nmodes_FEM.bin',form='unformatted')
!----------------------------------------------------------------------
! UNIT=intrp_unit Made in interp, has the data for
! interpolating surface fns.
!----------------------------------------------------------------------
OPEN(Unit=intrp_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Interp_FEM.bin', form='unformatted', status='unknown')
!  --------------------------------------------------------------
! UNIT vec_unit=vector contains the coefficients
!         of the surface functions in
!         the primitive basis of the ABM.
!  --------------------------------------------------------------
OPEN(Unit=vec_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Vector_ABM.bin', form='unformatted', status='unknown')

RETURN
END
