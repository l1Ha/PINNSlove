MODULE Dipole_Module
   USE Numeric_Kinds_Module
   IMPLICIT NONE
   SAVE
   INTEGER(KIND=IW_Kind):: NSurf=1   !  nsurf is the number of electronic potential energy surfaces used.
   NAMELIST/dipole/nsurf
END MODULE Dipole_Module
