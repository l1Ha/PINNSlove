SUBROUTINE Kmatrx_Calc(ndelves, R_mat, RE_mat, K_mat, KE_mat, xk,   &
                   chinuj, wvec, nbasis, usys2, etot, jtot,         &
                   chanl, elect, nvib, jrot, lorb, xksq,            &
                   energy, nchanl, narran, temp, rmatch,            &
                   Out_Unit, A_mat, B_mat, E_mat, F_mat,            &
                   AE_mat, BE_mat, EE_mat, FE_mat,                  &
                   nopen, calc_kmat, EDeriv, njacobi, gamma,        &
                   phor, phoi, temp1, temp2)
USE melm_Module
USE region_Module
USE dip_Module
!
IMPLICIT NONE
LOGICAL EDeriv, calc_kmat
LOGICAL, PARAMETER :: Full=.False.
INTEGER nopen, ndelves, njacobi
INTEGER narran, nchanl(narran), jtot, Out_Unit, i
INTEGER elect(ndelves),  nvib(ndelves),   jrot(ndelves)
INTEGER  lorb(ndelves), chanl(ndelves), nbasis(ndelves)
REAL(Kind=WP_Kind) usys2, rmatch, etot
REAL(Kind=WP_Kind)  energy(ndelves), wvec(njacobi)
REAL(Kind=WP_Kind)    xksq(njacobi),   xk(njacobi)
REAL(Kind=WP_Kind)  R_mat(ndelves,ndelves), RE_mat(ndelves,ndelves)
REAL(Kind=WP_Kind)  K_mat(njacobi,njacobi), KE_mat(njacobi,njacobi)
REAL(Kind=WP_Kind)  A_mat(ndelves,njacobi),  B_mat(ndelves,njacobi)
REAL(Kind=WP_Kind)  E_mat(ndelves,njacobi),  F_mat(ndelves,njacobi)
REAL(Kind=WP_Kind) AE_mat(ndelves,njacobi), BE_mat(ndelves,njacobi)
REAL(Kind=WP_Kind) EE_mat(ndelves,njacobi), FE_mat(ndelves,njacobi)
REAL(Kind=WP_Kind) temp(njacobi,njacobi), chinuj(*)
REAL(Kind=WP_Kind) temp1(njacobi,njacobi)
REAL(Kind=WP_Kind) temp2(njacobi,njacobi)
REAL(Kind=WP_Kind) gamma(njacobi), phor(njacobi), phoi(njacobi)
!-----------------------------------------------------------------
!     Calculate the asymptotic jacobi basis functions.
!-------------------------------------------------------------------
IF(njacobi/=ndelves)THEN
   WRITE(Out_Unit,*)'currently njacobi must equal ndelves'
   WRITE(Out_Unit,*)'You must have exactly the same basis'
   WRITE(Out_Unit,*)'in the delves and jacobi regions'
   WRITE(Out_Unit,*)' Stopped in boundary'
   WRITE(Out_Unit,*)' ndelves= ', ndelves,'  njacobi= ',njacobi
   WRITE(*,*)'currently njacobi must equal ndelves'
   WRITE(*,*)'You must have exactly the same basis'
   WRITE(*,*)'in the delves and jacobi regions'
   WRITE(*,*)' Stopped in boundary'
   WRITE(*,*)' ndelves= ', ndelves,'  njacobi= ',njacobi
   STOP 'Stopping Inside Kmatrx_Calc'
ENDIF
   
nopen=njacobi
iregion='jacobi '
IF(Full)THEN
  WRITE(Out_Unit,*)'Quantum Numbers just before calling upsiln'
  CALL Quant_Out(ndelves, chanl, elect, nvib, jrot, lorb, energy, xksq,'KMatrx_Calc_1')
ENDIF

WRITE(Out_Unit,*)'Calculate jacobi basis functions in Kmatrx_Calc'
CALL upsiln(temp, A_mat, xksq, xk, rmatch, chinuj)
iregion='delves '
DO i=1,ndelves
   xksq(i)=usys2*(Etot-energy(i))
ENDDO
!-------------------------------------------------------------------
!     Calculate the K-matrix.
!-------------------------------------------------------------------
IF(Full)THEN
  WRITE(Out_Unit,*)'Inside Kmatrx_Calc Calling Quant_Out: Kmatrx_Full will be called'
  CALL Quant_Out(ndelves, chanl, elect, nvib, jrot, lorb, energy, xksq,'KMatrx_Calc_2')
ENDIF

WRITE(Out_Unit,*)'Inside Kmatrx_Calc Calling Kmatrx_Full: NJacobi=',njacobi
CALL Kmatrx_Full (ndelves, R_mat, RE_mat, K_mat, KE_mat, wvec,   &
                  chanl, elect, nvib, jrot, lorb, jtot, etot,    &
                  energy, usys2, Out_Unit, temp,                 &
                  A_mat,  B_mat,  E_mat,  F_mat,                 &
                  AE_mat, BE_mat, EE_mat, FE_mat,                &
                  calc_kmat, EDeriv,njacobi,gamma,phor,          &
                  phoi,temp1,temp2)
WRITE(Out_Unit,*)'Inside Kmatrx_Calc returned from calling Kmatrx_Full'
!-------------------------------------------------------------------
!     Determine the number of open channels and then reduce the
!     K-matrix to only include the open channels.
!-------------------------------------------------------------------

DO i=1,ndelves  !do I need this tmpmodgregparker
   xksq(i)=usys2*(Etot-energy(i))
ENDDO
IF(Full)THEN
  WRITE(Out_Unit,*)'Inside Kmatrx_Calc Calling Quant_Out: Kmatrx_Open will be called'
  CALL Quant_Out(ndelves, chanl, elect, nvib, jrot, lorb, energy, xksq,'KMatrx_Calc_3')
ENDIF

WRITE(Out_Unit,*)'Inside Kmatrx_Calc Calling Kmatrx_Open: NJacobi=',njacobi
CALL Kmatrx_Open (njacobi, nopen, K_mat, KE_mat, wvec, nbasis,    &
                  chanl, elect, nvib, jrot, lorb, xksq,           &
                  energy, nchanl, narran, Out_Unit, calc_kmat,    &
                  EDeriv, phor, phoi)
!-------------------------------------------------------------------
!     on return
!     nopen  = no. of open channels.
!     K_mat  = K_matrix
!     KE_mat = derivative of K_matrix with respect to energy.
!-------------------------------------------------------------------
RETURN
ENDSUBROUTINE Kmatrx_Calc
