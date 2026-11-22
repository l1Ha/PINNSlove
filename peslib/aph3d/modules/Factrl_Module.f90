MODULE Factrl_Module
   USE Numeric_Kinds_Module
   SAVE
   REAL(Kind=WP_Kind) ftl(101) ! Contains the first 100 factorials
   NAMELIST/factrl/ ftl
END MODULE Factrl_Module
