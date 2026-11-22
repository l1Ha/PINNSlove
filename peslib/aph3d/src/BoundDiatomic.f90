SUBROUTINE BoundDiatomic
USE Numeric_Kinds_Module
USE Numbers_Module
USE FileUnits_Module
USE Narran_Module
USE Masses_Module
USE InputFile_Module
USE Region_Module
USE Oops_Module
USE DiatomicPot_Module,  ONLY: DiatomicPot
USE Convrsns_Module,     ONLY: autoev
USE OneDim_Parms_Module, ONLY: nmax
!-----------------------------------------------------------------------
! This routine was written by G. A. Parker
! IF you find an error or have an improvement please send a messge to
! Parker@ou.edu
!   This routine will provide a daf fit to a set of DATA on a
!   uniform grid. See papers by David Hoffman and Donald J. Kouri
!   for details. One useful paper is by David K. Hoffman, Mark Arnold
!   Donald J. Kouri. J. Phys. Chem 96, 6539-6545 1992.
! This program requires the following input DATA:
! nradial=number of radial values.
! ndaf=largest desired derivative.
! mval=order of the daf fit. This is M. To improve the
!      fit increase M.
! NBoundDiatomic=number of channels IF f(x) is a vector or matrix.
!
! This program seems to work for several test cases.
!
!-----------------------------------------------------------------------
IMPLICIT NONE
CHARACTER(LEN=1)  :: ChLabel
CHARACTER(LEN=3)  :: system, group
CHARACTER(LEN=15) :: Identifier
CHARACTER(LEN=30) :: PotFilename
LOGICAL :: there, SetPotToZero, direct
INTEGER :: diag_opt, jtot, nsize_save, IArran, IRadial, jval, NBound, IMult, numrho, kmax, iblock(1), IERR
REAL(Kind=WP_Kind) :: dradial, dtheta, deltachi
REAL(Kind=WP_Kind) :: tht_val(1), chi_val(2), t_theta(1), t_chi(2), t_chir(2)
REAL(Kind=WP_Kind) :: a(2), b(2), c(2), d(2), t_chi1(2), t_chir1(2), r(3), tblock(1)
REAL(Kind=WP_Kind) :: RadialMax, VLast, VInfinity, VDiff, umass, umass2
INTEGER :: ndaf, NBoundDiatomic, mval, nd_daf, nradial, nsize
INTEGER, PARAMETER:: ntheta=1, nchi=1
REAL(Kind=WP_Kind), ALLOCATABLE::radial_val (:)
REAL(Kind=WP_Kind), ALLOCATABLE::dafx (:), hermit (:)
REAL(Kind=WP_Kind), ALLOCATABLE::wavefunction (:), toe (:), eig (:)
REAL(Kind=WP_Kind), ALLOCATABLE::v_pot (:), v_eff(:), centerm(:)
REAL(Kind=WP_Kind), ALLOCATABLE::t_radial (:)
REAL(Kind=WP_Kind), ALLOCATABLE::radialm2 (:), hamil (:)


NAMELIST/DAF_BoundDiatomic/ ndaf, NBoundDiatomic, mval, diag_opt, SetPotToZero, direct, jtot, nd_daf

NRadial=nmax
RadialMax=10.0d0
dradial=RadialMax/(NRadial-1)
NSize=NRadial

WRITE(Out_Unit, * )
WRITE(Out_Unit, * ) 'Begin(BoundDiatomic) Main Program'
nd_daf = 20

!-----------------------------------------------------------------------
! set up default DATA for for all input.
!
! ndaf=largest desired derivative.
! mval=order of the daf fit. This is M.
! NBoundDiatomic=number of desired eigenvalues and wavefunctions.
!-----------------------------------------------------------------------
direct = .true.
NBoundDiatomic = 1
mval = 20

ndaf = 2
Inquire (file = InputDIR(1:LEN(TRIM(InputDIR)))//InputFile, exist = there)
IF(there)THEN
   OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile, status = 'old')
   READ(In_Unit, NML=DAF_BoundDiatomic, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'DAF_BoundDiatomic')
      nd_daf=40
      ndaf=2
      mval=32
      NBoundDiatomic=13
      diag_opt=2
      SetPotToZero=.False.
      direct=.true.
      jtot=0
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist DAF_BoundDiatomic'
      WRITE(Msg_Unit,*)'ERROR with Namelist DAF_BoundDiatomic'
      STOP 'ReadAllData: DAF_BoundDiatomic'
   ENDIF
   WRITE(Out_Unit, * )
   WRITE(Out_Unit, * ) 'Begin(NAMELIST Input DATA=DAF_BoundDiatomic)'
   WRITE(Out_Unit,NML=DAF_BoundDiatomic)
   WRITE(Out_Unit, * ) 'END(NAMELIST Input DATA=DAF_BoundDiatomic)'
   WRITE(Out_Unit, * )
   CLOSE(unit=In_Unit, status = 'keep')
ELSE
   WRITE(Msg_Unit, * ) 'Error: .bound does not exist'
   WRITE(Out_Unit, * ) 'Error: .bound does not exist'
   STOP 'BoundDiatomic'
ENDIF

!-----------------------------------------------------------------------
! check input integers to make sure that ridiculous DATA is not present.
!-----------------------------------------------------------------------

CALL check_int (nsize, ndaf, NBoundDiatomic, mval)
!-----------------------------------------------------------------------
! This routine dynamically ALLOCATEs memory for each of the arrays.
!-----------------------------------------------------------------------
nsize=NRadial
ALLOCATE (radial_val (nradial) )
ALLOCATE (dafx ( (nd_daf + 1) * (ndaf + 1) ), hermit (mval + ndaf + 2) )
ALLOCATE (wavefunction ( nradial **2) )
ALLOCATE (toe (nd_daf), eig (nsize) )
ALLOCATE (v_pot (nradial * 2), v_eff(nradial), centerm(nradial))
ALLOCATE (t_radial (nradial * nradial))
ALLOCATE (radialm2 (nradial), hamil (nsize * nsize) )
nsize_save = nsize

WRITE(Out_Unit, * ) 'nsize_save=', nsize_save
!-----------------------------------------------------------------------
! Calculate the diatomic potential for each arrangement channel
!-----------------------------------------------------------------------
DO IArran=1,Narran
   umass=U_Diatom(IArran)
   umass2=2.d0*umass
   WRITE(Out_Unit,*)'masses=',mass(1),mass(2),mass(3)
   WRITE(Out_Unit,*)'Two Body reduced mass: umass=',umass
   WRITE(ChLabel,'(I1)')IArran
   PotFilename='DiatomicPlots/'//TRIM(DiatomicPot)//'_VChannel.csv'
   OPEN(Unit=VEE_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//PotFilename,Form='Formatted')
   ! Print Heading
   WRITE(VEE_unit,'(3(A16,","))')"r(au)","V_Pot(Ha)","v_Pot(eV)"
   DO IRadial=1,NRadial
      r(1)=2.154434690031884d12  ! PORT FIX: was Ten**(RANGE(One)/3)=1e102 in double precision -> PES overflow (NaN)
      r(2)=2.154434690031884d12  ! PORT FIX: was Ten**(RANGE(One)/3)=1e102 in double precision -> PES overflow (NaN)
      r(3)=2.154434690031884d12  ! PORT FIX: was Ten**(RANGE(One)/3)=1e102 in double precision -> PES overflow (NaN)
      Radial_Val(IRadial)=(IRadial-1)*dradial
      r(IArran)=Radial_Val(IRadial)
      IF(IRadial/=1)THEN
         Centerm(IRadial)=One/(umass2*r(IArran)**2)
      ELSE
         Centerm(IRadial)=Two/(umass2*dradial**2)
      ENDIF
      CALL Surface(v_pot(IRadial),r)
      IF(IRadial>=3)THEN
         IF(IRadial==3)THEN
            V_Pot(1)=V_Pot(2)+V_Pot(2)-V_Pot(3) !Avoid potentialsingularities at r(Irran)=0
            WRITE(VEE_unit,'(3(ES16.7,","))')Radial_Val(1), V_Pot(1), V_Pot(1)*autoev
            WRITE(VEE_unit,'(3(ES16.7,","))')Radial_Val(2), V_Pot(2), V_Pot(2)*autoev           
            WRITE(VEE_unit,'(3(ES16.7,","))')Radial_Val(3), V_Pot(3), V_Pot(3)*autoev
         ELSEIF(IRadial>3)THEN  
            WRITE(VEE_unit,'(3(ES16.7,","))')Radial_Val(IRadial), V_Pot(IRadial), V_Pot(IRadial)*autoev
         ENDIF
      ENDIF
   ENDDO
   CLOSE(Unit=VEE_unit)
   VLast=V_Pot(NRadial)
   VDIFF=10.0d0
   Imult=2
   DO WHILE (ABS(VDiff)>1.d-12)
      r(IArran)=Imult*RadialMax
      CALL Surface(VInfinity,r)
      VDIff=VInfinity-VLast
      VLast=VInfinity
      IMult=IMult+1
   ENDDO
   WRITE(Out_Unit,'(" r(Infinity) is less than or equal to ",ES16.7," and  VInfinity=",ES24.16)')r(IArran),VInfinity
   !-----------------------------------------------------------------------
   ! make sure that the DATA is uniformally spaced.
   !-----------------------------------------------------------------------
   IF(nradial>1)THEN
      CALL checkbound (radial_val, nradial, dradial)
   ENDIF
   !-----------------------------------------------------------------------
   ! This is where all of the work is done. A DAF of order M is used to
   ! construct the kinetic energy operator operator.
   !-----------------------------------------------------------------------

   PotFilename='DiatomicPlots/'//TRIM(DiatomicPot)//'_VChIntrp.csv'
   OPEN(Unit=VEE_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//PotFilename,Form='Formatted')
   t_radial=0.d0
   Identifier="_BoundDiatomic_"
   CALL get_daf (NRadial, radial_val, dradial, mval, dafx, hermit, ndaf,  nd_daf, toe, umass2, kmax, Identifier)
   CALL kin_rho (NRadial, toe, t_radial, Radial_Val, kmax)
   jval=0
   NBOUND=NRadial
   iregion='jacobi'
   ! Print Heading
   !WRITE(VEE_unit,'(4(A16,","))')"r(au)","v(Ha)","vp(Ha)","vpp(Ha)"
   DO WHILE (NBound>0)
      V_eff=V_Pot+jval*(jval+1)*Centerm
      group='C1'
      t_theta(1)=fortyfive
      t_chi(1)=30.d0
      CALL VProperties(NRadial, radial_val, V_eff, VInfinity, IArran, jval)
      CALL boundold (nradial, radial_val, ntheta, tht_val, nchi, chi_val, dradial, &
       dtheta, deltachi, wavefunction, mval, dafx, hermit, ndaf, NBoundDiatomic, &
       nd_daf, toe, eig, v_eff, diag_opt, nsize, t_radial, t_theta, t_chi, &
       t_chir, t_chi1, t_chir1, system, group, a, b, c, d, radialm2, jtot, &
       umass2, direct, tblock, iblock, hamil)
      jval=jval+1
      !WRITE(VEE_Unit,*)"&"
      IF(jval>0)NBound=0
   ENDDO
   CLOSE(VEE_Unit)
ENDDO
!-----------------------------------------------------------------------
! release all of the ALLOCATED memory and check to make sure the
! dimensions were not exceeded.
!-----------------------------------------------------------------------
WRITE(Out_Unit, * ) 'nsize_save=', nsize_save
DEALLOCATE(radial_val)
DEALLOCATE(dafx, hermit)
DEALLOCATE(wavefunction, toe, eig)
DEALLOCATE(v_pot, v_eff, centerm)
DEALLOCATE(t_radial)
DEALLOCATE(radialm2, hamil)

WRITE(Out_Unit, * ) 'END(BoundDiatomic) Calculations Competed'
WRITE(Out_Unit, * )

RETURN
ENDSUBROUTINE BoundDiatomic
