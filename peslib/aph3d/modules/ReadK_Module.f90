MODULE ReadK_Module
   USE Numeric_Kinds_Module
   USE TotalEng_Module
   USE QAsy_Numbers_Module
   IMPLICIT NONE
   SAVE
   INTEGER(Kind=IW_Kind) ithmat
   INTEGER(Kind=IW_Kind), ALLOCATABLE :: EminLoc(:), LastLoc(:), ElasLoc(:), InelasLoc(:), ReactLoc(:)
END MODULE ReadK_Module