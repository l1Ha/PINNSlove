SUBROUTINE sym_check (hamil, nsize)
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
INTEGER :: i, j, nsize
REAL(Kind=WP_Kind), PARAMETER::errtol = 1.d-12

REAL(Kind=WP_Kind) :: hamil (nsize, nsize), err, maxerr
WRITE(Out_Unit, * )
WRITE(Out_Unit, * ) 'Begin(Sym_Check) Symmetry Check'
maxerr = 0.d0
DO i = 1, nsize-1
   DO j = i + 1, nsize
      err = abs (hamil (i, j) - hamil (j, i) )
      IF(err>maxerr) maxerr = err
      IF(err>errtol)THEN
   WRITE(Out_Unit, * ) 'Error: i, j, err=', i, j, err
   WRITE(Out_Unit, * ) 'hamil(i,j)=', hamil (i, j)
   WRITE(Out_Unit, * ) 'hamil(j,i)=', hamil (j, i)
      ENDIF
   ENDDO

ENDDO
IF(maxerr<errtol)THEN
   WRITE(Out_Unit, * ) 'Hamiltonian is symmetric'
ELSE
   WRITE(Out_Unit, * ) 'Hamiltonian is NOT symmetric'
   WRITE(Out_Unit, * ) 'Maximum asymmetry=', maxerr
!cc         STOP 'sym_check'
ENDIF
WRITE(Out_Unit, * ) 'END(Sym_Check)'

WRITE(Out_Unit, * )
RETURN
ENDSUBROUTINE sym_check
