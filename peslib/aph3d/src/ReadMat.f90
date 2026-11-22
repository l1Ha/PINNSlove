SUBROUTINE ReadMat(amat, energies, chanl, elect, nvib, jrot, lorb, n, inunit)
USE FileUnits_Module
IMPLICIT NONE
INTEGER n, inunit, chanl(n), elect(n), nvib(n), jrot(n), lorb(n)
REAL(Kind=WP_Kind) amat(n, n),  energies(n)
CHARACTER(LEN=11) label
CHARACTER(LEN=21) end_label
CHARACTER(LEN=21) start_label

!
! Read in the energies.
!
!write(*,*)"n pos1=",n
READ(inunit)start_label, energies, end_label
!write(*,*)"n pos2=",n
!write(*,*)"energies=",energies
label='Energies   '
CALL check_label(start_label, end_label, label)

!
! Read in the Channel quantum numbers.
!
READ(inunit)start_label, chanl, end_label
label='Channel    '
CALL check_label(start_label, end_label, label)

!
! Read in the Electronic quantum numbers.
!
READ(inunit)start_label, elect, end_label
label='Electronic '
CALL check_label(start_label, end_label, label)

!
! Read in the Vibrational quantum numbers.
!
READ(inunit)start_label, nvib, end_label
label='Vibrational'
CALL check_label(start_label, end_label, label)

!
! Read in the Rotational quantum numbers.
!
READ(inunit)start_label, jrot, end_label
label='Rotational '
CALL check_label(start_label, end_label, label)

!
! Read in the Projection of rotational quantum numbers.
!
READ(inunit)start_label, lorb, end_label
label='Orbital    '
CALL check_label(start_label, end_label, label)

!
! Read in the Matrix.
!
READ(inunit)start_label, amat, end_label
label='Matrix     '
CALL check_label(start_label, end_label, label)

RETURN
ENDSUBROUTINE ReadMat
