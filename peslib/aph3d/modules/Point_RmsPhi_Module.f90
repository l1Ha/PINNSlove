MODULE Point_RmsPhi_Module
   USE Numeric_Kinds_Module
   IMPLICIT NONE
   SAVE
   INTEGER(KIND=IW_Kind), ALLOCATABLE:: IntTheta(:) ! th(numnp)
   INTEGER(KIND=IW_Kind), ALLOCATABLE:: IntChi(:)   ! ch(numnp)
   REAL(Kind=WP_Kind),    ALLOCATABLE:: phi(:,:)    ! phi(numnp, nq)
   REAL(Kind=WP_Kind),    ALLOCATABLE:: phirms(:,:) ! phirms(ntheta, nchi)
END MODULE Point_RmsPhi_Module
