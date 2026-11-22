SUBROUTINE Bound3D
USE Numeric_Kinds_Module
USE FileUnits_Module
USE InputFile_Module
USE PotZero_Module
!-----------------------------------------------------------------------
! This routine was written by G. A. Parker
! If you find an error or have an improvement please send a messge to
! Parker@ou.edu
!   This routine will provide a daf fit to a set of data on a
!   uniform grid. See papers by David Hoffman and Donald J. Kouri
!   for details. One useful paper is by David K. Hoffman, Mark Arnold
!   Donald J. Kouri. J. Phys. Chem 96, 6539-6545 1992.
! This program requires the following input data:
! nrho=number of rho values.
! ntheta=number of theta values.
! nchi=number of chi values.
! ndaf=largest desired derivative.
! mval=order of the daf fit. This is M. To improve the
!      fit increase M.
! NBound3D=number of channels if f(x) is a vector or matrix.
!
! This program seems to work for several test cases.
!
!-----------------------------------------------------------------------
IMPLICIT NONE
CHARACTER(LEN=3) :: system, group
LOGICAL :: there, SetPotToZero, direct, morph
INTEGER :: diag_opt, jtot, nsize_save, imorph, start_morph, nmorph, IERR
REAL(Kind=WP_Kind) :: usys2, drho, dtheta, deltachi, morph_val
INTEGER :: ndaf, NBound3D, mval, nd_daf, nrho, ntheta, nchi
INTEGER :: nsize,ibug
INTEGER, ALLOCATABLE :: iblock (:)
REAL(Kind=WP_Kind), ALLOCATABLE::rho_val (:), tht_val (:), chi_val (:)
REAL(Kind=WP_Kind), ALLOCATABLE::dafx (:), hermit (:)
REAL(Kind=WP_Kind), ALLOCATABLE::wavefunction (:), toe (:), eig (:)
REAL(Kind=WP_Kind), ALLOCATABLE::v_pot (:)
REAL(Kind=WP_Kind), ALLOCATABLE::t_rho (:), t_theta (:)
REAL(Kind=WP_Kind), ALLOCATABLE::t_chi (:), t_chir (:)
REAL(Kind=WP_Kind), ALLOCATABLE::a (:), b (:), c (:), d (:)
REAL(Kind=WP_Kind), ALLOCATABLE::t_chi1 (:), t_chir1 (:)
REAL(Kind=WP_Kind), ALLOCATABLE::tblock (:)
REAL(Kind=WP_Kind), ALLOCATABLE::rhom2 (:), hamil (:)


NAMELIST/DAF_Bound3d/ ndaf, NBound3D, mval, diag_opt, SetPotToZero, direct, jtot, nd_daf

WRITE(Out_Unit,*)'Called    Bound3D'
CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/Bound3D.txt',Form='FORMATTED')

WRITE(Out_Unit, * )
WRITE(Out_Unit, * ) 'Begin(Bound3D) Main Program'
morph = .false.
start_morph = 1
nmorph = 1
IF(morph)THEN
   start_morph = 1
   nmorph = 3
ENDIF
nd_daf = 20
!-----------------------------------------------------------------------
! set up default data for for all input.
!
! ndaf=largest desired derivative.
! mval=order of the daf fit. This is M.
! NBound3D=number of desired eigenvalues and wavefunctions.
!-----------------------------------------------------------------------
direct = .true.
NBound3D = 1
mval = 20

ndaf = 2
DO imorph = start_morph, nmorph
   IF(imorph==1)THEN
      morph_val = 2.2
   ENDIF
   IF(imorph==2)THEN
      morph_val = 4.975
   ENDIF
   IF(imorph==3)THEN
      morph_val = 9.5
   ENDIF
   OPEN(Unit=pA_mot_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'GraphicsOut/pA.mot')
   OPEN(Unit=pB_mot_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'GraphicsOut/pB.mot')
   OPEN(Unit=pC_mot_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'GraphicsOut/pC.mot')
   CALL calc_pot (morph, morph_val)
   CLOSE(unit=pA_mot_Unit)
   CLOSE(unit=pB_mot_Unit)
   CLOSE(unit=pC_mot_Unit)

   CALL v_read (nrho, ntheta, nchi, v_pot, rho_val, tht_val, chi_val, usys2, .false., nsize, system, group)

   INQUIRE(file = TRIM(InputDIR)//TRIM(InputFile), exist = there)
   IF(there)THEN
      OPEN(Unit=In_Unit,File=TRIM(InputDIR)//TRIM(InputFile), status = 'old')
      READ(In_Unit,NML=DAF_Bound3d, IOSTAT=IERR)
      IF(IERR==-1)THEN
         CALL NameList_Default(Out_Unit, 'DAF_Bound3d')
         nd_daf=40
         ndaf=2
         mval=32
         NBound3D=30
         diag_opt=2
         SetPotToZero=.true.
         direct=.true.
         jtot=0
      ELSEIF(IERR/=0)THEN
         WRITE(Out_Unit,*)'ERROR with Namelist DAF_Bound3d'
         WRITE(Msg_Unit,*)'ERROR with Namelist DAF_Bound3d'
         STOP 'Bound3D: DAF_Bound3d'
      ENDIF
      WRITE(Out_Unit, * )
      WRITE(Out_Unit, * ) 'Begin(NAMELIST Input DATA=DAF_Bound3d)'
      WRITE(Out_Unit,NML=DAF_Bound3d)
      WRITE(Out_Unit, * ) 'END(NAMELIST Input DATA=DAF_Bound3d)'
      WRITE(Out_Unit, * )
      CLOSE(unit=In_Unit, status = 'keep')
   ELSE
      WRITE(Msg_Unit, * ) 'Error: .bound does not exist'
      WRITE(Out_Unit, * ) 'Error: .bound does not exist'
      STOP 'Bound3D'
   ENDIF

   WRITE(Out_Unit, '(A,ES22.15)' ) ' usys2=', usys2
   !-----------------------------------------------------------------------
! check input integers to make sure that ridiculous data is not
   ! present.
   !-----------------------------------------------------------------------

   CALL check_int (nsize, ndaf, NBound3D, mval)
   !-----------------------------------------------------------------------
! This routine dynamically allocates memory for each of the arrays.
   !-----------------------------------------------------------------------
   ALLOCATE (rho_val (nrho), tht_val (ntheta), chi_val (nchi + 1) )
   ALLOCATE (dafx ( (nd_daf + 1) * (ndaf + 1) ), hermit (mval + ndaf + 2) )
!   ALLOCATE (wavefunction ( (nrho * ntheta * nchi) **2) )
   ALLOCATE (wavefunction (2) )
   ALLOCATE (toe (nd_daf), eig (nsize) )
   ALLOCATE (v_pot (nrho * ntheta * nchi * 2) )
   ALLOCATE (t_rho (nrho * nrho), t_theta (ntheta * ntheta) )
   ALLOCATE (t_chi (nchi * nchi), t_chir (nchi * nchi) )
   ALLOCATE (a (ntheta), b (ntheta), c (ntheta), d (ntheta) )
   ALLOCATE (t_chi1 (nchi * nchi), t_chir1 (nchi * nchi) )
   ALLOCATE (iblock (nchi * nchi), tblock (nchi * nchi) )
   ALLOCATE (rhom2 (nrho) )
   nsize_save = nsize


   WRITE(Out_Unit, * ) 'nsize_save=', nsize_save
   !-----------------------------------------------------------------------
   ! READ the potential
   !-----------------------------------------------------------------------

   CALL v_read (nrho, ntheta, nchi, v_pot, rho_val, tht_val, chi_val, usys2, .true., nsize, system, group)
   IF(SetPotToZero)THEN
      WRITE(Out_Unit, * ) 'WARNING setting potential to zero'
      WRITE(Out_Unit, * ) 'nsize=', nrho * ntheta * (nchi + 1)
      v_pot = 0.d0
   ENDIF
   !-----------------------------------------------------------------------
   ! make sure that the DATA is uniformally spaced.
   !-----------------------------------------------------------------------
   IF(nrho>1)THEN
      CALL checkbound (rho_val, nrho, drho)
   ENDIF
   IF(ntheta>1)THEN
!Notice here that the theta points are constructed in DVR method
      !CALL checkbound (tht_val, ntheta, dtheta)
   ENDIF
   IF(nchi>1)THEN
      CALL checkbound (chi_val, nchi, deltachi)
   ENDIF
   !-----------------------------------------------------------------------
   ! This is where all of the work is done. A DAF of order M is used to
   ! construct the kinetic energy operator operator.
   !-----------------------------------------------------------------------
   CALL bound (nrho, rho_val, ntheta, tht_val, nchi, chi_val, drho, &
    dtheta, deltachi, wavefunction, mval, dafx, hermit, ndaf, NBound3D, &
    nd_daf, toe, eig, v_pot, diag_opt, nsize, t_rho, t_theta, t_chi, &
    t_chir, t_chi1, t_chir1, system, group, a, b, c, d, rhom2, jtot, &
    usys2, direct, tblock, iblock)
   !-----------------------------------------------------------------------
! release all of the allocated memory and check to make sure the
   ! dimensions were not exceeded.
   !-----------------------------------------------------------------------
   WRITE(Out_Unit, * ) 'nsize_save=', nsize_save
   DEALLOCATE (rho_val, tht_val, chi_val)
   DEALLOCATE (dafx, hermit)
   DEALLOCATE (wavefunction, toe, eig)
   DEALLOCATE (v_pot)
   DEALLOCATE (t_rho, t_theta)
   DEALLOCATE (t_chi, t_chir)
   DEALLOCATE (a, b, c, d)
   DEALLOCATE (t_chi1, t_chir1)
   DEALLOCATE (iblock, tblock)
   DEALLOCATE (rhom2)
   !      CLOSE(Unit=Out_Unit)
ENDDO
WRITE(Out_Unit, * ) 'END(Main_Bound) Calculations Competed'
CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/BasisForScatt.txt',Form='FORMATTED',ACCESS='APPEND')
WRITE(Out_Unit,*)'Completed Bound3D'
WRITE(Out_Unit,*)

RETURN
ENDSUBROUTINE Bound3D
