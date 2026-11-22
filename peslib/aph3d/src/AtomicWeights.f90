SUBROUTINE AtomicWeights(AtomicSymbol, MassNumber, AtomicNumber, AtomicMass, Abundance, AtomicWeight, Notes)
USE Numeric_Kinds_Module       ! Set numeric kinds
USE FileUnits_Common_Module    ! Fileunit numbers
   !
   ! Input: Atomic Symbol and Mass Number
   ! Output: Mass Number, Relative Atomic Mass, Isotopic Composition, Standard Atomic Weight, and Notes
   ! DATA was compiled by NIST
   !

IMPLICIT NONE
SAVE
INTEGER, PARAMETER :: NumMasses=3180

! Input: 
CHARACTER(LEN=3),    INTENT(IN) :: AtomicSymbol
INTEGER,             INTENT(IN) :: MassNumber

! Output:
CHARACTER(LEN=100), INTENT(OUT) :: Notes
INTEGER,            INTENT(OUT) :: AtomicNumber
REAL(KIND=WP_Kind), INTENT(OUT) :: AtomicMass
REAL(KIND=WP_Kind), INTENT(OUT) :: Abundance
REAL(KIND=WP_Kind), INTENT(OUT) :: AtomicWeight

LOGICAL ReadData, IDBug
CHARACTER(LEN=3) ASymbolTable(NumMasses)
CHARACTER(LEN=100) Line
CHARACTER(LEN=100) NoteTable(NumMasses)
INTEGER I
INTEGER AtomicNumberTable(NumMasses)
INTEGER MassNumberTable(NumMasses)
REAL(KIND=WP_Kind) AtomicMassTable(NumMasses)
REAL(KIND=WP_Kind) AtomicWeightTable(NumMasses)
REAL(KIND=WP_Kind) AbundanceTable(NumMasses)

DATA ReadData/.True./, IDBug/.False./

IF(ReadData)THEN
    OPEN(In_Unit, File=TRIM(InputRoot)//'AtomicWeights/AtomicWeights.txt',Form='FORMATTED')
READ(In_Unit,'(A)') Line
IF(IDBug)WRITE(Msg_Unit,*) Line
DO I=1, NumMasses
   READ(In_Unit,'(A16,I2)') Line, AtomicNumberTable(I)
   IF(IDBug)WRITE(Msg_Unit,*) Line, AtomicNumberTable(I)
   READ(In_Unit,'(A16,A)') Line, ASymbolTable(I)
   IF(IDBug)WRITE(Msg_Unit,*) Line, ASymbolTable(I)
   READ(In_Unit,'(A14,I10)') Line, MassNumberTable(I)
   IF(IDBug)WRITE(Msg_Unit,*) Line, MassNumberTable(I)
   READ(In_Unit,'(A23,D23.15)') Line, AtomicMassTable(I)
   IF(IDBug)WRITE(Msg_Unit,*) Line, AtomicMassTable(I)
   READ(In_Unit,'(A23,D23.15)') Line, AbundanceTable(I)
   IF(IDBug)WRITE(Msg_Unit,*) Line, AbundanceTable(I)
   READ(In_Unit,'(A25,D23.15)') Line, AtomicWeightTable(I)
   IF(IDBug)WRITE(Msg_Unit,*) Line, AtomicWeightTable(I)
   READ(In_Unit,'(A8,A)') Line, NoteTable(I)
   IF(IDBug)WRITE(Msg_Unit,*) Line, NoteTable(I)
   READ(In_Unit,'(A)') Line
   IF(IDBug)WRITE(Msg_Unit,*) Line
   IF(AtomicMassTable(I)==0.d0)AtomicMassTable(I)=AtomicWeightTable(I)
ENDDO
READ(In_Unit,'(A)') Line
IF(IDBug)WRITE(Msg_Unit,*) Line
ReadData=.False.
ENDIF

DO I=1, NumMasses
   IF(MassNumber==MassNumberTable(I).and.TRIM(ADJUSTL(AtomicSymbol))==TRIM(ADJUSTL(ASymbolTable(I))))THEN
      AtomicNumber=AtomicNumberTable(I)
      AtomicMass=AtomicMassTable(I)
      Abundance=AbundanceTable(I)
      AtomicWeight=AtomicWeightTable(I)
      Notes=NoteTable(I)
      EXIT
   ENDIF
   IF(I==NumMasses)THEN
      WRITE(Msg_Unit,*)"Incorrect Mass Number or Atomic Symbol", MassNumber, AtomicSymbol
      STOP "AtomicWeights"
   ENDIF
ENDDO

RETURN
END SUBROUTINE AtomicWeights
