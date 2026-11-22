MODULE Point_Picture_Module
   USE Numeric_Kinds_Module
   IMPLICIT NONE
   SAVE
   INTEGER(KIND=IW_Kind), ALLOCATABLE:: igloss(:)  ! igloss(nelemx)
   INTEGER(KIND=IW_Kind), ALLOCATABLE:: kolr(:)    ! kolr(nelemx)
   REAL(Kind=WP_Kind), ALLOCATABLE:: center(:)  ! center(3*nelemx)
   REAL(Kind=WP_Kind), ALLOCATABLE:: rad(:)     ! rad(nelemx)
END MODULE Point_Picture_Module
