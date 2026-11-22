MODULE Point_Sfunc_Module
   USE Numeric_Kinds_Module
   USE Nidm_Module
   IMPLICIT NONE
   SAVE
   INTEGER(KIND=IW_Kind), ALLOCATABLE:: lm(:)      ! lm(18*nel)
   !INTEGER(KIND=IW_Kind), ALLOCATABLE:: idi(:,:)   ! idi(ndisce,nidm)
   INTEGER(KIND=IW_Kind), ALLOCATABLE:: maxa(:)    ! maxa(neq+1)
   INTEGER(KIND=IW_Kind), ALLOCATABLE:: mht(:)     ! mht(neq) column heights.
   INTEGER(KIND=IW_Kind), ALLOCATABLE:: nloc(:)    ! nloc(neq)
   INTEGER(KIND=IW_Kind), ALLOCATABLE:: nsit(:)    ! nsit(nq)
   REAL(Kind=WP_Kind), ALLOCATABLE:: yz(:)      ! yz(18*nel)
   REAL(Kind=WP_Kind), ALLOCATABLE:: bghaml(:)  ! bghaml(nwk)
   INTEGER(KIND=IW_Kind), ALLOCATABLE:: ar(:)      ! ar(nnc)
   REAL(Kind=WP_Kind), ALLOCATABLE:: br(:)      ! br(nnc)
   REAL(Kind=WP_Kind), ALLOCATABLE:: vec(:,:)   ! vec(nq,nq)
   REAL(Kind=WP_Kind), ALLOCATABLE:: eigv(:)    ! eigv(nq)
   REAL(Kind=WP_Kind), ALLOCATABLE:: d(:)       ! d(nq)
   REAL(Kind=WP_Kind), ALLOCATABLE:: tt(:)      ! tt(neq)
   REAL(Kind=WP_Kind), ALLOCATABLE:: w(:)       ! w(neq)
   REAL(Kind=WP_Kind), ALLOCATABLE:: rtolv(:)   ! rtolv(nq)
   REAL(Kind=WP_Kind), ALLOCATABLE:: evc1(:)    ! evc1(nq)
   REAL(Kind=WP_Kind), ALLOCATABLE:: evc2(:)    ! evc2(nq)
   REAL(Kind=WP_Kind), ALLOCATABLE:: ww(:)      ! ww(neq)
END MODULE Point_Sfunc_Module
