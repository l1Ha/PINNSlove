PROGRAM Program_Test_Asymptotic
USE FileUnits_Asymptotic_Module
USE QAsy_Numbers_Module
USE Masses_Module
IMPLICIT NONE
LOGICAL, PARAMETER:: H3_Calc=.True.
INTEGER i

!------Required Information--------------------------------
InputDir=TRIM(InputRoot)//"Testing/Test_Asymptotic/"  ! Input directory
OutDir  =TRIM(OutputRoot)//"Testing/Test_Asymptotic/" ! Output directory 
NQuant=5
NConserved=4
NAtoms=3
ALLOCATE(Conserved(NConserved))
ALLOCATE(ConLabels(NConserved))
ALLOCATE(QLabels(NQuant))
AtomicSymbol=['F','H','H']
MassNumber=[19,1,1]
DO I=1,NAtoms
   CALL AtomicWeights(AtomicSymbol(i), MassNumber(i), AtomicNumber(i), AtomicMass(i), Abundance(i), AtomicWeight(i), Notes(i))
ENDDO
Conserved=[0,1,3,1]
ConLabels=['Jtot','Parity','C2v','A1']
QLabels=['chanl', 'elec', 'nvib', 'jrot', 'lorb']
CALL AsymptoticTest

STOP 'Test_Asymptotic'
END PROGRAM Program_Test_Asymptotic