MODULE VFunc_Module
   USE Numeric_Kinds_Module
   USE Narran_Module
   SAVE
   REAL(Kind=WP_Kind), ALLOCATABLE ::  chinuj(:)  !Vibrational wavefunctions for each rotational state evaluated at the Gauss_Hermite
!      NAMELIST/vfunc/ chinuj
END MODULE VFunc_Module
