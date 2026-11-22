MODULE Arrch_Module
   USE Numeric_Kinds_Module
   SAVE
   INTEGER(KIND=IW_Kind):: marran    ! Mth arrangement channel
   INTEGER(KIND=IW_Kind):: karran    ! Kth arrangement channel
   NAMELIST/Arrch/ marran, karran
END MODULE Arrch_Module
