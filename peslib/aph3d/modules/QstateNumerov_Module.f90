MODULE QstateNumerov_Module
   USE Numeric_Kinds_Module
   USE Qstate_Module
   IMPLICIT NONE
   SAVE
   REAL(KIND=WP_KIND), ALLOCATABLE :: EigNumerov(:), LeftTP(:), RightTP(:)
   INTEGER, ALLOCATABLE :: NthElec(:), Spin(:), Lambda(:), Nuvib(:), jrot(:), YIJ(:,:)
   LOGICAL, ALLOCATABLE :: Bound(:), QuasiB(:), Scatt(:)
END MODULE QstateNumerov_Module