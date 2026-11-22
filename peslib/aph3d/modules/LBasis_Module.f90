MODULE LBasis_Module
   USE Parms_Module
   IMPLICIT NONE
   SAVE
   INTEGER(KIND=IW_Kind) lmin(narran)   ! minimum value of l for Legendre polynomial.
   INTEGER(KIND=IW_Kind) lskip(narran)  ! should be 1 for channels with no symmetry and 2 for channels with symmetry.
END MODULE LBasis_Module