SUBROUTINE title
USE FileUnits_Module
IMPLICIT NONE
CHARACTER(LEN=16), PARAMETER:: spc='                '
WRITE(Out_Unit,*)'----------Start Title----------'
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)spc, '   aa     pppp   h    h  33333   dddd'
WRITE(Out_Unit,*)spc, '  a  a   p    p  h    h 3     3  d   d'
WRITE(Out_Unit,*)spc, ' a    a  p    p  h    h       3  d    d'
WRITE(Out_Unit,*)spc, ' aaaaaa  pppppp  hhhhhh  33333   d    d'
WRITE(Out_Unit,*)spc, ' a    a  p       h    h       3  d    d'
WRITE(Out_Unit,*)spc, ' a    a  p       h    h 3     3  d   d'
WRITE(Out_Unit,*)spc, ' a    a  p       h    h  33333   dddd'
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)spc, spc, 'APH3D $Revision: 1.14 $'
WRITE(Out_Unit,*)spc, spc, '$Date: 89/11/30 11:28:52 $'
WRITE(Out_Unit,*)spc, spc, '$State: Stable $'
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)spc, spc, spc, ' by G.A. Parker'
WRITE(Out_Unit,*)spc, spc, spc, '    R.T  Pack'
WRITE(Out_Unit,*)spc, spc, spc, '    B.J. Archer'
WRITE(Out_Unit,*)spc, spc, spc, '    J.D. Kress'
WRITE(Out_Unit,*)spc, spc, spc, '    R.B. Walker'
WRITE(Out_Unit,*)spc, spc, spc, '    S.B. Padjker'
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'----------End Title------------'
RETURN
END
