MODULE VIVS_Module
   USE Numeric_Kinds_Module
   IMPLICIT NONE
   SAVE
   INTEGER(KIND=4)    :: IAlpha = 4     ! Number of Steps per Interval
   REAL(Kind=WP_Kind) :: drnow  = 0.1d0 ! Current Step Size
   REAL(Kind=WP_Kind) :: toff   = 1.d-4 ! Off-diagonal tolerance
   NAMELIST/vivs/ drnow, toff, ialpha
END MODULE VIVS_Module
