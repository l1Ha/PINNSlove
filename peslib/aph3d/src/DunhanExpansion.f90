SUBROUTINE DunhamExpansion(nvib_qn, jrot_qn, jth, jth_max, kth, kth_max, term, E_nu_j, Yjk)
USE Numeric_Kinds_Module
USE Numbers_Module
IMPLICIT NONE
CHARACTER(LEN=21), PARAMETER:: ProcName='DunhamExpansion'
INTEGER nvib_qn, jrot_qn, jth, kth, j, jth_max, k, kth_max
REAL(Kind=WP_Kind) :: term, E_nu_j, Yjk(0:jth_max, 0:kth_max)

!Single term in the Dunham Expansion
term=(nvib_qn+Half)**jth*(jrot_qn*(jrot_qn+One))**kth

!Sum of Coefficients (Yjk) times the Dunham term
IF(jth_max>0.and.kth_max>0)THEN
  DO j=0,jth_max
    DO k=0, kth_max
      E_nu_j=E_nu_j+Yjk(j, k)*(nvib_qn+Half)**j*(jrot_qn*(jrot_qn+One))**k
    ENDDO
  ENDDO
ENDIF

RETURN
ENDSUBROUTINE DunhamExpansion