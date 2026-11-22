subroutine read_unformatted_matrix(file,rho,a,n,np)
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
logical there
integer n,i,j,unit,np
REAL(Kind=WP_Kind) rho
REAL(Kind=WP_Kind) a(n,np)
CHARACTER(LEN=*) file

OPEN(Unit=UFMat_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//file,form='unformatted')
rewind(UFMat_Unit)
read(UFMat_Unit) rho,a
close(UFMat_Unit)

return
end
