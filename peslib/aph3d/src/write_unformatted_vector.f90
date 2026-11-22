subroutine WRITE_unformatted_vector(file,rho,a,n)
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
logical there
integer n,i,j,unit
REAL(Kind=WP_Kind) rho
REAL(Kind=WP_Kind) a(n)
CHARACTER(LEN=*) file

OPEN(Unit=UFMat_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//file,form='unformatted')
rewind(UFMat_Unit)
WRITE(UFMat_Unit) rho,a
call flush(UFMat_Unit)
close(UFMat_Unit)

return
end
