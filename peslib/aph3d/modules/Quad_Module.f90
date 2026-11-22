MODULE Quad_Module
   USE Numeric_Kinds_Module
   SAVE
   INTEGER(KIND=IW_Kind) :: nquad = 5 ! Number of quadrature points
   REAL(Kind=WP_Kind) xg(96)          ! Gauss_Legendre quadrature points.
   REAL(Kind=WP_Kind) wgt(96)         ! Gauss_Legendre weights.
   NAMELIST/quad/ xg, wgt, nquad
END MODULE Quad_Module
