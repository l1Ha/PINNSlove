SUBROUTINE mxoutf (a, nrow, ncol, nrowd, ncold)
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
INTEGER nrow, ncol, irow, jcol, nrowd, ncold
REAL(Kind=WP_Kind) a(nrowd,ncold)

WRITE(Out_Unit,*)'Number of Rows=',NRow,' Number of Columns=',ncol
WRITE(Out_Unit,*)'Dimension of Rows=',NRowd,' Dimension of Columns=',ncold
irow = 1
WRITE(Out_Unit,230)(jcol,jcol=1,min(ncol,10))
WRITE(Out_Unit,240)irow, (a(1,jcol),jcol=1,min(ncol,10))
DO irow = 2, min(nrow,10)
   WRITE(Out_Unit,240)irow, (a(irow,jcol),jcol=1,min(ncol,10))
ENDDO

RETURN
 230  FORMAT(/, 3x, '  row   ', 10(3x, i3, 6x))
 240  FORMAT(1x , 1x, 'col ', i3, 2x, 10(1pe12.4))
 250  FORMAT(/, 1x, 'col ', i3, 2x, 10(1pe12.4))
END
