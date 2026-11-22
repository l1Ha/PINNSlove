MODULE Point_Quad_Module
   USE Numeric_Kinds_Module
   IMPLICIT NONE
   SAVE
   INTEGER(KIND=IW_Kind), ALLOCATABLE:: lambda(:)      ! lambda(nq)
   REAL(Kind=WP_Kind), ALLOCATABLE:: xpt(:)         ! xpt(nquad)
   REAL(Kind=WP_Kind), ALLOCATABLE:: wht(:)         ! wht(nquad)
   REAL(Kind=WP_Kind), ALLOCATABLE:: value(:, :)    ! value(nquad, nquad)
   REAL(Kind=WP_Kind), ALLOCATABLE:: h(:, :, :)     ! h(9, nquad, nquad)
   REAL(Kind=WP_Kind), ALLOCATABLE:: p(:, :, :, :)  ! p(2, 9, nquad, nquad)
   REAL(Kind=WP_Kind), ALLOCATABLE:: hcoef(:, :, :) ! hcoef(9, nquad, nquad)
   REAL(Kind=WP_Kind), ALLOCATABLE:: thetaval(:)    ! thetaval(ntheta)
   REAL(Kind=WP_Kind), ALLOCATABLE:: chivals(:)     ! chivals(nchi)
   REAL(Kind=WP_Kind), ALLOCATABLE:: ovrlapmx(:)    ! ovrlapmx(nq)
   REAL(Kind=WP_Kind), ALLOCATABLE:: energy(:)      ! energy(nq)
   REAL(Kind=WP_Kind), ALLOCATABLE:: energy1(:)     ! energy1(nq)
END MODULE Point_Quad_Module
