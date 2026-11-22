MODULE Point_Msher_Module
   USE Numeric_Kinds_Module
   USE Nidm_Module
   IMPLICIT NONE
   SAVE
   INTEGER(KIND=IW_Kind), ALLOCATABLE :: lnod(:)      ! lnod(2*nel)
   INTEGER(KIND=IW_Kind), ALLOCATABLE :: iscr1(:)     ! iscr1(ntheta*(nchi/2+1))
   INTEGER(KIND=IW_Kind), ALLOCATABLE :: iscr2(:)     ! iscr2(nel)
   INTEGER(KIND=IW_Kind), ALLOCATABLE :: nid(:)       ! nid(2*numnp/3)
   INTEGER(KIND=IW_Kind), ALLOCATABLE :: nodj(:)      ! nodj(nidm+1, 2*numnp/3)
   INTEGER(KIND=IW_Kind), ALLOCATABLE :: idteq(:)     ! idteq(numnp)
   INTEGER(KIND=IW_Kind), ALLOCATABLE :: idiv(:)      ! idiv(nel)
   INTEGER(KIND=IW_Kind), ALLOCATABLE :: lcor(:)      ! lcor(2*nel)
   INTEGER(KIND=IW_Kind), ALLOCATABLE :: mapr(:)      ! mapr(nel)
   INTEGER(KIND=IW_Kind), ALLOCATABLE :: idivo(:)     ! idivo(nel)
   INTEGER(KIND=IW_Kind), ALLOCATABLE :: lcoro(:)     ! lcoro(2*nel)
   INTEGER(KIND=IW_Kind), ALLOCATABLE :: mapnod(:)    ! mapnod(numnp)
   INTEGER(KIND=IW_Kind), ALLOCATABLE :: iwork33(:)   ! iwork33(2*nel)
   REAL(Kind=WP_Kind),    ALLOCATABLE :: scr1(:)      ! scr1(numnp)
   REAL(Kind=WP_Kind),    ALLOCATABLE :: scr2(:)      ! scr2(numnp)
   REAL(Kind=WP_Kind),    ALLOCATABLE :: beta(:,:)    ! beta(nidm,2*numnp/3)
   REAL(Kind=WP_Kind),    ALLOCATABLE :: sinlist(:)   ! sinlist(nel)
   REAL(Kind=WP_Kind),    ALLOCATABLE :: philist(:)   ! philist(nel)
   INTEGER, ALLOCATABLE:: Virtualx(:,:) ! Virtualx(ntheta, nchi)
END MODULE Point_Msher_Module
