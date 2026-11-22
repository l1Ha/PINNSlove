MODULE Boundary_Module
   USE Numeric_Kinds_Module
   IMPLICIT NONE
   SAVE
   INTEGER(KIND=IW_Kind), PARAMETER:: Max_Sector=3000     ! Maximum number of sectors
   CHARACTER(LEN=10):: Basis_Type(Max_Sector)             ! Type of basis method used
   INTEGER(KIND=IW_Kind):: IthStart_ABM                   ! Staring sector for ABM Method
   INTEGER(KIND=IW_Kind):: IthEnd_ABM                     ! Ending sector for ABM Method
   INTEGER(KIND=IW_Kind):: IthStart_DVR                   ! Starting sector for DVR Method
   INTEGER(KIND=IW_Kind):: IthEnd_DVR                     ! Ending sector for DVR method
   INTEGER(KIND=IW_Kind):: IthStart_FEM                   ! Staring sector for FEM method
   INTEGER(KIND=IW_Kind):: IthEnd_FEM                     ! Ending sector for FEM method
   INTEGER(KIND=IW_Kind):: IthStart_PDAF                  ! Staring sector for PDAF method
   INTEGER(KIND=IW_Kind):: IthEnd_PDAF                    ! Ending sector for PDAF method
   INTEGER(KIND=IW_Kind):: LastRho=0                      ! Last sector sector number
   INTEGER(KIND=IW_Kind):: CurrentSector=1                ! The current sector
   INTEGER(KIND=IW_Kind):: Kase_CalcParm=1                ! Used in CalcParm to determine sector boundaries
   REAL(Kind=WP_Kind) Basis_Dist(Max_Sector)              ! Distance for basis set calculations
   REAL(Kind=WP_Kind) Sector_Boundaries(Max_Sector+1)     ! Sector Boundaries
END MODULE Boundary_Module
