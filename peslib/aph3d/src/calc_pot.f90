SUBROUTINE calc_pot (morph, morph_val)
USE Numeric_Kinds_Module
USE FileUnits_Module
USE PES_MODULE
USE InputFile_Module
IMPLICIT NONE
INTEGER :: nrho, ntheta, nchi, IERR
REAL(Kind=WP_Kind), ALLOCATABLE::rho_val (:), tht_val (:), chi_val (:), v_pot (:)
LOGICAL :: morph
REAL(Kind=WP_Kind) :: morph_val

NAMELIST/pot/ nrho, ntheta, nchi
WRITE(Out_Unit, * )
WRITE(Out_Unit, * ) 'Begin(Calc_Pot)'
WRITE(Out_Unit, * )
WRITE(Out_Unit, * ) 'Begin(READ data)', InputDIR(1:LEN(TRIM(InputDIR)))//TRIM(InputFile)

OPEN(Unit=In_Unit,File=TRIM(InputDIR)//TRIM(InputFile))
READ(In_Unit,NML=Pot, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'Pot')
      NRho=15
      NTheta=21
      NChi=25
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Pot'
      WRITE(Msg_Unit,*)'ERROR with Namelist Pot'
      STOP 'Calc_Pot: Pot'
   ENDIF
CLOSE(In_Unit)
WRITE(Out_Unit,NML=pot)

WRITE(Out_Unit, * ) PES_Name
WRITE(Out_Unit, * ) 'END(READ data):', InputDIR(1:LEN(TRIM(InputDIR)))//InputFile


WRITE(Out_Unit, * )
ALLOCATE (rho_val (nrho), tht_val (ntheta), chi_val (nchi + 1), v_pot (nrho * ntheta * nchi) )
CALL calc_pot2 (nrho, ntheta, nchi, rho_val, tht_val, chi_val, v_pot, morph, morph_val)
DEALLOCATE (rho_val, tht_val, chi_val, v_pot)
!      CLOSE(Unit=Out_Unit)
WRITE(Out_Unit, * ) 'END(Calc_Pot)'

WRITE(Out_Unit, * )
RETURN
ENDSUBROUTINE calc_pot
