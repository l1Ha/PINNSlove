MODULE Fuzzy_Module
   USE Numeric_Kinds_Module
   IMPLICIT NONE
   SAVE
   REAL(Kind=WP_Kind) fuzz ! Set to a small number so that if two real numbers differ by
                           ! less than this amount they will be considered equal.
END MODULE Fuzzy_Module
