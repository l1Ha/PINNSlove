MODULE Mshdat_Module
   USE Numeric_Kinds_Module
   USE Extra_Module
   USE MData_Module
   USE Main_Module
   USE Converge_Module
   IMPLICIT NONE
   SAVE
   INTEGER(KIND=IW_Kind) :: nave1_FEM       !
   INTEGER(KIND=IW_Kind) :: nave2_FEM       !
   INTEGER(KIND=IW_Kind) :: Nconsidr  = 246        !
   REAL(Kind=WP_Kind)    :: tolexp                 !
   REAL(Kind=WP_Kind)    :: tolfn1     = 0.8D0     !
   REAL(Kind=WP_Kind)    :: tolfn2     = 5.0D-02   !
   REAL(Kind=WP_Kind)    :: aterm                  !
   REAL(Kind=WP_Kind)    :: bterm                  !
   REAL(Kind=WP_Kind)    :: xkmesh       = 12.22D0 !
   NAMELIST/mshdat/ nlt, nlc, ndiv, nfreq, nxtra, nodmax,tolexp, tolfn2, xkmesh, nblock, maxj, nlancz,            &
                         tolfn1, nconsidr, aterm, bterm, nave1_FEM, nave2_FEM
END MODULE Mshdat_Module
