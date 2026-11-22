MODULE Point_Ovr_Module
   USE Numeric_Kinds_Module
   IMPLICIT NONE
   SAVE
   INTEGER(KIND=IW_Kind), ALLOCATABLE:: idlast(:,:)   ! idlast(9, nel1)
   INTEGER(KIND=IW_Kind), ALLOCATABLE:: igrid(:)      ! igrid(2*ntheta*nchi)
   REAL(Kind=WP_Kind), ALLOCATABLE:: fb(:, :, :)   ! fb(nq, nquad, nquad)
   REAL(Kind=WP_Kind), ALLOCATABLE:: fbp(:, :, :)  ! fbp(nq, nquad, nquad)
   REAL(Kind=WP_Kind), ALLOCATABLE:: sovrlap(:, :) ! sovrlap(nq, nq)
   INTEGER(KIND=IW_Kind), ALLOCATABLE:: thlast(:)     ! thlast(numnp1)
   INTEGER(KIND=IW_Kind), ALLOCATABLE:: chlast(:)     ! chlast(numnp1)
   REAL(Kind=WP_Kind), ALLOCATABLE:: philast(:)    ! philast(max(numnp,numnp1)*nq)
   REAL(Kind=WP_Kind), ALLOCATABLE:: f3(:)         ! f3(numnp)
   REAL(Kind=WP_Kind), ALLOCATABLE:: chi1(:)       ! chi1(numnp)
   REAL(Kind=WP_Kind), ALLOCATABLE:: ThetAPH(:)    ! ThetAPH(numnp)
END MODULE Point_Ovr_Module
