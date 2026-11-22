SUBROUTINE Read_K2(n, amat, InUnit, OutUnit)
USE Numeric_Kinds_Module
USE ReadK_Module
USE QAsy_Numbers_Module, ONLY: AsymptEnergies
USE Masses_Module
IMPLICIT NONE
LOGICAL PrintOpt/.False./
INTEGER n, i, k, NQuant_TMP, NConserved_Read
INTEGER InUnit, OutUnit
REAL(Kind=WP_Kind) amat(n, n)
CHARACTER(LEN=21) end_label
CHARACTER(LEN=21) start_label
CHARACTER(LEN=7), PARAMETER:: ProcName='Read_K2'
CHARACTER(LEN=6) Print_Flag


CALL PoptAsy(ProcName, Print_Flag)
IF(Print_Flag=='Full')PrintOpt=.True.
!
!Read in the number of conserved quantum numbers
start_label = 'Start NConserved     '
end_label   = 'End NConserved       '
READ(InUnit)start_label, NConserved_Read, end_label
IF(PrintOpt)WRITE(OutUnit,'(A,I5,A)')start_label, NConserved_Read, end_label
!
! Read Conserved
start_label = 'Start Conserved      '
end_label =   'End Conserved        '
READ(InUnit)start_label, Conserved, end_label
IF(PrintOpt)WRITE(OutUnit,*)start_label, Conserved, end_label
!
! Read Conserved Labels
start_label = 'Start ConLabels      '
end_label =   'End ConLabels        '
READ(InUnit)start_label, ConLabels, end_label
IF(PrintOpt)WRITE(OutUnit,*)start_label, ConLabels, end_label
!
! Read Energies
start_label = 'Start Energies       '
end_label =   'End Energies         '
READ(InUnit)    start_label, (AsymptEnergies(i),i=1,n), end_label
end_label =   'End Energies         '   !tmpmodgregparker should fix write statement
IF(PrintOpt)WRITE(OutUnit,*)start_label, (AsymptEnergies(i),i=1,n), end_label
!
! Write in the energies.
start_label = 'Start NAtoms         '
end_label =   'End NAtoms           '
READ(InUnit)start_label, NAtoms, end_label
IF(PrintOpt)WRITE(OutUnit,*)start_label, NAtoms, end_label
!
! Read information about the masses
start_label = 'Start MassInformation'
end_label =   'End MassInformation  '
READ(InUnit)start_label, AtomicNumber, AtomicSymbol, MassNumber, AtomicMass, abundance, AtomicWeight,Notes, end_label
IF(PrintOpt)WRITE(OutUnit,*)start_label, AtomicNumber, end_label
IF(PrintOpt)WRITE(OutUnit,*)start_label, AtomicSymbol, end_label
IF(PrintOpt)WRITE(OutUnit,*)start_label, MassNumber, end_label
IF(PrintOpt)WRITE(OutUnit,*)start_label, AtomicMass, end_label
IF(PrintOpt)WRITE(OutUnit,*)start_label, Abundance, end_label
IF(PrintOpt)WRITE(OutUnit,*)start_label, AtomicWeight, end_label
IF(PrintOpt)WRITE(OutUnit,*)start_label, Notes, end_label
!
! Read two times the system reduced mass
start_label = 'Start usys2          '
end_label =   'End usys2            '
READ(InUnit)start_label, usys2, end_label
IF(PrintOpt)WRITE(OutUnit,*)start_label, usys2, end_label
!
! Read wavevector
start_label = 'Start wavevector     '
end_label =   'End wavevector       '
READ(InUnit)start_label, (wvec(i),i=1,n), end_label
IF(PrintOpt)WRITE(OutUnit,*)start_label, (wvec(i),i=1,n), end_label
!
! Read xksq
start_label = 'Start xksq           '
end_label =   'End xksq             '
READ(InUnit)start_label, (xksq(i),i=1,n), end_label
end_label =   'End xksq             '  !tmpmodgregparker should fix write statement
IF(PrintOpt)WRITE(OutUnit,*)"pos1: ",start_label, (xksq(i),i=1,n), end_label
!
! Read the Channel quantum numbers.
start_label = 'Start Channel        '
end_label =   'End Channel          '
READ(InUnit)    start_label, NQuant_TMP, end_label
IF(NQuant/=NQuant_TMP)THEN
    WRITE(*,*)'Stopping in Read_K2: NQuant/=NQuant_TMP', NQuant, NQuant_TMP
    STOP 'NQuant/=NQuant_TMP'
ENDIF
IF(PrintOpt)WRITE(OutUnit,*)start_label, NQuant, end_label
!
! Read the Channel quantum numbers.
start_label = 'StartQLabels         '
end_label =   'End Channel          '
READ(InUnit)    start_label, QLabels, end_label
IF(PrintOpt)WRITE(OutUnit,*)start_label, QLabels, end_label

QNumbr=0
!
! Read the Electronic quantum numbers.
start_label = 'Start Electronic     '
end_label =   'End Electronic       '
READ(InUnit)    start_label, ((QNumbr(k,i),k=1,NQuant),i=1,n), end_label
IF(PrintOpt)THEN
   WRITE(OutUnit,*)start_label
   DO k=1,NQuant
      WRITE(OutUnit,'(A,I5)')'QNumbr: k=',k
      WRITE(OutUnit,'(20I5)') (QNumbr(k,i),i=1,n)
   ENDDO
   WRITE(OutUnit,*)End_label
ENDIF
!
! Read the Matrix.
start_label = 'Start Matrix         '
end_label =   'End Matrix           '
READ(InUnit)start_label, amat, end_label
IF(PrintOpt)WRITE(OutUnit,*)start_label, amat, end_label

RETURN
ENDSUBROUTINE Read_K2
