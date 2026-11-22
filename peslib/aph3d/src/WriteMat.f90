SUBROUTINE WriteMat(amat, n, outunit)
USE Numeric_Kinds_Module
USE FileUnits_Asymptotic_Module
USE QAsy_Numbers_Module
USE Masses_Module
IMPLICIT NONE
INTEGER n, k, i, outunit
REAL(Kind=WP_Kind) amat(n, n)
CHARACTER(LEN=21) end_label
CHARACTER(LEN=21) start_label
CHARACTER(LEN=21), PARAMETER:: ProcName='WriteMat'
CHARACTER(LEN=6) Print_Flag

CALL PoptAsy(ProcName, Print_Flag)
!
! Write number of conserved quantum numbers
start_label = 'Start NConserved     '
end_label   = 'End NConserved       '
WRITE(outunit)Start_Label, NConserved, end_label
!
! Write Conserved
start_label = 'Start Conserved      '
end_label =   'End Conserved        '
WRITE(outunit)Start_Label, Conserved, end_label
!
! Write Conserved
start_label = 'Start ConLabels      '
end_label =   'End ConLabels        '
WRITE(outunit)Start_Label, ConLabels, end_label
!
! Write the energies.
start_label = 'Start Energies       '
end_label =   'End Energies         '
WRITE(outunit)Start_Label, AsymptEnergies, end_label
!
! Write NAtoms.
start_label = 'Start NAtoms         '
end_label =   'End NAtoms           '
WRITE(outunit)Start_Label, NAtoms, end_label
!
! Write information about the masses
start_label = 'Start MassInformation'
end_label =   'End MassInformation  '
WRITE(outunit)Start_Label, AtomicNumber, AtomicSymbol, MassNumber, AtomicMass, abundance, AtomicWeight, Notes, end_label
!
! Write two times the system reduced mass
start_label = 'Start usys2          '
end_label =   'End usys2            '
WRITE(outunit)Start_Label, usys2, end_label
!
! Write wavevector
start_label = 'Start wavevector     '
end_label =   'End wavevector       '
WRITE(outunit)Start_Label, (wvec(i),i=1,n), end_label
!
! Write xksq
start_label = 'Start xksq           '
end_label =   'End xksq             '
WRITE(outunit)"pos2 ",Start_Label, (xksq,i=1,n), end_label
!
! Write in the Channel quantum numbers.
start_label = 'Start NQuant         '
end_label =   'End Nquant           '
WRITE(outunit)Start_Label, NQuant, end_label
!
! Write Quantum number labels
start_label = 'Start QLabels        '
end_label =   'End QLabels          '
WRITE(outunit)Start_Label, QLabels, end_label
!
! Write in the Electronic quantum numbers.
start_label = 'Start QNumbr         '
end_label =   'End QNumbr           '
WRITE(outunit)Start_Label, ((QNumbr(k,i),k=1,NQuant),i=1,n), end_label
!
! Write in the Matrix.
start_label = 'Start Matrix         '
end_label =   'End Matrix           '
WRITE(outunit)Start_Label, amat, end_label

RETURN
ENDSUBROUTINE WriteMat
