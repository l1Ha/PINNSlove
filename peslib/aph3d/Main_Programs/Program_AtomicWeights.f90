PROGRAM Program_AtomicWeights
USE Numeric_Kinds_Module       ! Set numeric kinds
USE FileUnits_Module           ! Fileunit numbers
   !
   ! Test Program for AtomicWeights.
   !

IMPLICIT NONE
INTEGER, PARAMETER:: NAtoms=3
CHARACTER(LEN=3) AtomicSymbol(NAtoms)
CHARACTER(LEN=100) Line(NAtoms), Notes(NAtoms)
INTEGER I
INTEGER AtomicNumber(NAtoms), MassNumber(NAtoms)
REAL(KIND=WP_Kind) AtomicMass(NAtoms), Abundance(NAtoms), AtomicWeight(NAtoms)

!CALL NumericInformation
AtomicSymbol=['Li','F','H']
MassNumber  =[   7, 19,  1]
DO I=1,NAtoms
   CALL AtomicWeights(AtomicSymbol(I), MassNumber(I), AtomicNumber(I), AtomicMass(I), Abundance(I), AtomicWeight(I), Notes(I))
   WRITE(Msg_Unit, '(A,A3)')'Atomic Symbol = ', TRIM(AtomicSymbol(I))
   WRITE(Msg_Unit, '(A,I3)')'Atomic Number = ', AtomicNumber(I)
   WRITE(Msg_Unit, '(A,I3)')'Mass Number = '  , MassNumber(I)
   WRITE(Msg_Unit, '(A,1PES23.15)')'Relative Atomic Mass = '  , AtomicMass(I)
   WRITE(Msg_Unit, '(A,1PES23.15)')'Isotopic Composition = '  , Abundance(I)
   WRITE(Msg_Unit, '(A,1PES23.15)')'Standard Atomic Weight = ', AtomicWeight(I)
   WRITE(Msg_Unit, '(A,A)')'Notes = ', TRIM(Notes(I))
   WRITE(Msg_Unit,*)
ENDDO


WRITE(Msg_Unit,*)"Successful Calculation of AtomicWeights"
STOP "End of Program_AtomicWeights"
END PROGRAM Program_AtomicWeights