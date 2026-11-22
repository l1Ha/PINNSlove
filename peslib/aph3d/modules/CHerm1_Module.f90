MODULE CHerm1_Module
   USE Narran_Module
   USE Numeric_Kinds_Module
   SAVE
   REAL(Kind=WP_Kind), ALLOCATABLE :: chlast(:,:,:) ! Previous expansion coefficients
!      NAMELIST/cherm1/ chlast
END MODULE CHerm1_Module
