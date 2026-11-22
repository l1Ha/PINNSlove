PROGRAM Program_Test_Propagation
USE FileUnits_Module
USE TotJ_Module
IMPLICIT NONE
LOGICAL, PARAMETER:: H3_Calc=.True.
REAL(Kind=WP_Kind) rho
rho=0.d0

OutDir=TRIM(OutputRoot)//"Testing/Test_Propagation/"
InputDir=TRIM(InputRoot)//"Testing/Test_Propagation/"

CALL ReadAllData
OPEN(Unit=SBFrst_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/SetBasis_First.bin', form='unformatted')
WRITE(SBFrst_Unit)
WRITE(SBFrst_Unit)
WRITE(SBFrst_Unit)
WRITE(SBFrst_Unit)
WRITE(SBFrst_Unit)jtot, parity, mega, megamin, megamax, symmetry, jeven
WRITE(SBFrst_Unit)rho
CLOSE(Unit=SBFrst_Unit)
WRITE(0,*)"Test_Propagation not working"
CALL APHlogder

STOP 'Test_Propagation'
END PROGRAM Program_Test_Propagation