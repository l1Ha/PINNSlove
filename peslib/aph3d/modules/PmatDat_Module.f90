MODULE PmatDat_Module
USE Numeric_Kinds_Module
        IMPLICIT NONE
        SAVE
        REAL(Kind=WP_Kind) rhos(3)                  !
        REAL(Kind=WP_Kind), ALLOCATABLE::  pm1(:,:) !
        REAL(Kind=WP_Kind), ALLOCATABLE::  pm2(:,:) !
        REAL(Kind=WP_Kind), ALLOCATABLE::  e2d(:)   !
END MODULE PmatDat_Module