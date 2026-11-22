MODULE FBR_Module
   USE Numeric_Kinds_Module
   USE System_Module,Only: NSfunc
   SAVE
   INTEGER(KIND=IW_Kind) :: lmax=50    ! maximum order of Legendre polynomial in theta for dvr
   INTEGER(KIND=IW_Kind) :: mmax=100   ! maximum order of trignometric function in chi for dvr
   INTEGER(KIND=IW_Kind) :: nsym=1     ! symmetry of the trginometric functions in chi for dvr
   NAMELIST/fbr/ lmax, mmax, nsfunc, nsym
END MODULE FBR_Module
