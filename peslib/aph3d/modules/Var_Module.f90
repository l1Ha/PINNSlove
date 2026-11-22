MODULE Var_Module
   USE Numeric_Kinds_Module
   SAVE
   INTEGER(KIND=IW_Kind) nnn(7)   !
   INTEGER(KIND=IW_Kind) kpri     !
   NAMELIST/var/ nnn, kpri
END MODULE Var_Module