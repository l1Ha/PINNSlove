MODULE TotalEng_Module
   USE Numeric_Kinds_Module
   SAVE
   INTEGER(KIND=IW_Kind) IthEnergy
   REAL(Kind=WP_Kind)Etot   ! Total scattering energy.
   NAMELIST/totaleng/ Etot
END MODULE TotalEng_Module
