SUBROUTINE WriteMatOld(amat, energies, chanl, elect, nvib, jrot, lorb, n, outunit)
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
!#include <Inc.unitno>
INTEGER n, outunit, chanl(n), elect(n), nvib(n), jrot(n), lorb(n)
REAL(Kind=WP_Kind) amat(n, n),  energies(n)
CHARACTER(LEN=21) end_label
CHARACTER(LEN=21) start_label
!
! Write in the energies.
!
start_label = 'Start Energies       '
end_label =   'End Energies         '
WRITE(outunit)start_label, energies, end_label
!
! Write in the Channel quantum numbers.
!
start_label = 'Start Channel        '
end_label =   'End Channel          '
WRITE(outunit)start_label, chanl,    end_label
!
! Write in the Electronic quantum numbers.
!
start_label = 'Start Electronic     '
end_label = 'End Electronic         '
WRITE(outunit)start_label, elect,    end_label
!
! Write in the Vibrational quantum numbers.
!
start_label = 'Start Vibrational    '
end_label = 'End Vibrational        '
WRITE(outunit)start_label, nvib,     end_label
!
! Write in the Rotational quantum numbers.
!
start_label = 'Start Rotational     '
end_label = 'End Rotational         '
WRITE(outunit)start_label, jrot,     end_label
!
! Write in the Orbital quantum numbers.
!
start_label = 'Start Orbital        '
end_label = 'End Orbital            '
WRITE(outunit)start_label, lorb,     end_label
!
! Write in the Matrix.
!
start_label = 'Start Matrix         '
end_label = 'End Matrix             '
WRITE(outunit)start_label, amat,     end_label

RETURN
ENDSUBROUTINE WriteMatOld
