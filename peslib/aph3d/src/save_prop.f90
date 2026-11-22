SUBROUTINE save_prop(rho, RN, VecOld, Tn, n, nprop)
USE Numeric_Kinds_Module
USE FileNames_Module
IMPLICIT NONE
INTEGER n, nprop
REAL(Kind=WP_Kind) rho
REAL(Kind=WP_Kind) RN(nprop,nprop), VecOld(n,nprop)
REAL(Kind=WP_Kind) Tn(n)

CALL WRITE_unformatted_matrix(work_dir//rnrec_name,rho,RN,nprop,nprop)
CALL WRITE_unformatted_matrix(work_dir//vecoldrec_name,rho,VecOld,n,nprop)
CALL WRITE_unformatted_vector(work_dir//tnrec_name,rho,Tn,n)

RETURN
END
