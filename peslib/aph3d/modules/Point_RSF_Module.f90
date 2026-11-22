MODULE Point_RSF_Module
   USE Numeric_Kinds_Module
   USE Nidm_Module
   IMPLICIT NONE
   SAVE
   INTEGER(KIND=IW_Kind),   PARAMETER :: Nine=9        ! 9
   INTEGER(KIND=IW_Kind), ALLOCATABLE :: idi(:,:)      ! idi(9, nquad, nquad)
   !INTEGER(KIND=IW_Kind), ALLOCATABLE :: nid(:)        ! nid(ntheta)
   INTEGER(KIND=IW_Kind), ALLOCATABLE :: id(:,:)       ! id(9, nel)
   INTEGER(KIND=IW_Kind), ALLOCATABLE :: idni(:)       ! idni(numnp)
   REAL(Kind=WP_Kind),    ALLOCATABLE :: th(:,:)       ! th(9, nel)
   REAL(Kind=WP_Kind),    ALLOCATABLE :: ch(:,:)       ! ch(9, nel)
   REAL(Kind=WP_Kind),    ALLOCATABLE :: f(:,:)        ! f(numnp, nmode)
   !REAL(Kind=WP_Kind),    ALLOCATABLE :: beta(:,:)     ! beta(nidm, ndisce)
   REAL(Kind=WP_Kind),    ALLOCATABLE :: freq(:)       ! freq(nmode)
   REAL(Kind=WP_Kind),    ALLOCATABLE :: xn(:)         ! xn(numnp)
END MODULE Point_RSF_Module