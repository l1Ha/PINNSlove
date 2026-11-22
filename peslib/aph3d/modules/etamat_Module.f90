MODULE etamat_Module
!-----------------------------------------------------------------------
! Inc.etamat
!
!  etamx   Body-fixed Coriolis coupling matrix.
!-----------------------------------------------------------------------
USE Numeric_Kinds_Module
USE nstate_Module
IMPLICIT NONE
   SAVE
   REAL(Kind=WP_Kind)etamx
      COMMON / etamat / etamx(nstate*(nstate+1)/2)
END MODULE etamat_Module
