PROGRAM Program_Test_KMatrix
USE FileUnits_Asymptotic_Module
USE QAsy_Numbers_Module
USE Masses_Module
IMPLICIT NONE
LOGICAL, PARAMETER:: H3_Calc=.True.
INTEGER i
INTEGER, PARAMETER :: NOpen=90

!------Required Information--------------------------------
OutDir=TRIM(OutputRoot)//"Testing/Test_KMatrix/"       ! Output directory 
InputDir=TRIM(InputRoot)//"Testing/Test_KMatrix/"  ! Input directory
!NQuant=5
NConserved=4
NAtoms=3
ALLOCATE(Conserved(NConserved))
ALLOCATE(ConLabels(NConserved))
ALLOCATE(QLabels(NQuant))
CALL ReadAllData
Conserved=[0,1,3,1]
ConLabels=['Jtot','Parity','C2v','A1']
QLabels=['chanl', 'elec', 'nvib', 'jrot', 'lorb']
WRITE(0,*)"Test_KMatrix not working"
CALL KMatrixTest

STOP 'Test_KMatrix'
END PROGRAM Program_Test_KMatrix