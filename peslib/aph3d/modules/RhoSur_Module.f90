MODULE RhoSur_Module
   USE Numeric_Kinds_Module
   SAVE
   REAL(Kind=WP_Kind),  ALLOCATABLE:: curve(:,:)    !curve(NCurves,NSectors)
   REAL(Kind=WP_Kind) rhosurf    ! Hyperradius at surface functions are being calculated.
   NAMELIST/rhosur/ rhosurf
END MODULE RhoSur_Module
