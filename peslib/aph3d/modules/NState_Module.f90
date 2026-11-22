MODULE NState_Module
   USE Numeric_Kinds_Module
   IMPLICIT NONE
   SAVE
   INTEGER(KIND=IW_Kind), PARAMETER:: NState=3000  ! PORT: raised 200->1200->3000 for large CC spaces
   INTEGER(KIND=IW_Kind)::            Nopen_Max
END MODULE NState_Module