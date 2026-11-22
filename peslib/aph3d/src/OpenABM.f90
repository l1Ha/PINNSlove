SUBROUTINE OpenABM(jtot)
!-----------------------------------------------------------------------
! This routine opens MOST files used by the code.
!-----------------------------------------------------------------------
USE FileUnits_Module
IMPLICIT NONE
INTEGER jtot
!-----------------------------------------------------------------------
! UNIT=pmat_unit=Pmatrx_ABM is an unformatted file for
!        storing the coupling matrix elements.
!-----------------------------------------------------------------------
OPEN(Unit=pmat_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/PotMatrx_ABM.bin',form='unformatted',status='unknown')

IF(jtot/=0)THEN
!-----------------------------------------------------------------------
! UNIT=jtot_unit=diag_jtot_ABM is an unformatted file
!        for storing the coupling matrix elements.
!-----------------------------------------------------------------------
   OPEN(Unit=jtot_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Diag_Jtot_ABM.bin', form='unformatted', status='unknown')
!-----------------------------------------------------------------------
! UNIT=lam_unit=diag_lambda_ABM is an unformatted file for
!        storing the coupling matrix elements.
!-----------------------------------------------------------------------
   OPEN(Unit=lam_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Diag_Lambda_ABM.bin', form='unformatted', status='unknown')
!-----------------------------------------------------------------------
! UNIT=asym_unit=asym_top_ABM is an unformatted file for
!        storing the coupling matrix elements.
!-----------------------------------------------------------------------
   OPEN(Unit=asym_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Asym_Top_ABM.bin', form='unformatted', status='unknown')
!-----------------------------------------------------------------------
! UNIT=cor_unit=coriolis_ABM is an unformatted file
!        for storing the coupling matrix elements.
!-----------------------------------------------------------------------
   OPEN(Unit=cor_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Coriolis_ABM.bin',form='unformatted', status='unknown')
!-----------------------------------------------------------------------
! UNIT=Ovr_Unit=Ovrlp_ABM is an unformatted file for storing the overlap
!        matrix elements.
!-----------------------------------------------------------------------
ENDIF

OPEN(Unit=Ovr_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Ovrlp_ABM.bin', form='unformatted', status='unknown')

!-----------------------------------------------------------------------
! UNIT=sf33_unit=surf33 is an unformatted file for storing the coeffs
!                of the surface fcns for use by ovrbas and matbas.
!-----------------------------------------------------------------------
OPEN(Unit=sf33_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Surf33_ABM.bin',   form='unformatted', status='unknown')

!-----------------------------------------------------------------------
! UNIT=sf34_unit surf34 is an unformatted file for storing the coeffs
!                of the surface fcns for use by ovrbas and matbas.
!                Used at alternating rho with unit 33.
!-----------------------------------------------------------------------
OPEN(Unit=sf34_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Surf34_ABM.bin',   form='unformatted',  status='unknown')

!-----------------------------------------------------------------------
! UNIT=phip_unit=phiprim is an unformatted file for storing the primitive basis for use by matbas.
!-----------------------------------------------------------------------
OPEN(Unit=phip_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/PhiPrim_ABM.bin',   form='unformatted', status='unknown')

!-----------------------------------------------------------------------
! UNIT=phid_unit=phideri is an unformatted file for storing the chi
!             derivative of the primitive basis for use by matbas.
!-----------------------------------------------------------------------
OPEN(Unit=phid_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/PhiDeri_ABM.bin',   form='unformatted', status='unknown')

!  ---------------------------------------------------------------------
!  UNIT=vec_unit=vector is an unformatted file used to store the surface
!           function coefficients at the most recent rho.  It is
!           opened and closed in sfunbas.
!  ---------------------------------------------------------------------


!  ---------------------------------------------------------------------
!  UNIT=sflvl_unit=Sfelevl_ABM.rbw is a formatted file used to
!           store the surface function energy levels in eV.
!  ---------------------------------------------------------------------
OPEN(Unit=sflvl_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'GraphicsOut/Sfelevl_ABM.rbw',form='formatted',   status='unknown')
OPEN(Unit=sflvl_unit2,File=OutDIR(1:LEN(TRIM(OutDIR)))//'GraphicsOut/Sfelevl_ABM.csv', form='formatted',  status='unknown')
RETURN
ENDSUBROUTINE OpenABM
