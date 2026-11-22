MODULE BasisRegions_Module
   USE Numeric_Kinds_Module
   USE Boundary_Module
   LOGICAL :: JacobiProp = .False.  ! Propagate in Jacobi coordinates
   LOGICAL :: FEMBasis   = .False.  ! Use FEM method to calculate surface functions
   LOGICAL :: DVRBasis   = .False.  ! Use DVR method to calculate surface functions
   LOGICAL :: PDAFBasis  = .False.  ! Use PDAF method to calculate surface functions
   INTEGER(KIND=IW_Kind) :: NStartSectors   =0      !
   INTEGER(KIND=IW_Kind) :: NAPH_ABMSectors =0      ! Number of ABM-APH sectors
   INTEGER(KIND=IW_Kind) :: NDelvesSectors  =0      ! Number of Delves sectors
   INTEGER(KIND=IW_Kind) :: NSectors     =31        ! Number to sectors
   REAL(Kind=WP_Kind) :: RhoStart        =2.1000D0  ! Minimum sector distance
   REAL(Kind=WP_Kind) :: RhoStartABM_APH =2.1000D0  ! Starting distance for ABM Method
   REAL(Kind=WP_Kind) :: RhoStartDelves  =8.1000D0  ! Starting distance for Delves propagation
   REAL(Kind=WP_Kind) :: RStartJacobi    =8.1000D0  ! Starting distance for Jacobi propagation
END MODULE BasisRegions_Module
