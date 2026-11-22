SUBROUTINE Coupling_Matrix(nchanl, wmat, sfpotmat, eigen_new,ksquared, usys2, rho, rho_basis)
USE Numeric_Kinds_Module
USE numbers_Module
IMPLICIT NONE
CHARACTER(LEN=21), PARAMETER:: ProcName='Coupling_Matrix'
INTEGER nchanl, i, j
REAL(Kind=WP_Kind) wmat(nchanl,nchanl), sfpotmat(nchanl,nchanl), fac
REAL(Kind=WP_Kind) ksquared, eigen_new(nchanl), usys2, rho, rho_basis
!  add zeroth order energies to coupling matrix and subtract the
!  1/8murho**2 term of 98/4/30-3 notes.
!  all these are 2*mu times the usual energies
!      Write(Msg_Unit,*)'rho_basis,rho,usys2=',rho_basis,rho,usys2
fac=rho_basis**2/rho**2
DO i=1,nchanl
   DO j=1,i
      wmat(i,j)=sfpotmat(i,j)
      IF(i==j)THEN
      wmat(i,j)=wmat(i,j)+eigen_new(i)*usys2*fac-one/(four*rho**2) -ksquared
      ENDIF
      wmat(j,i)=wmat(i,j)
   ENDDO
ENDDO

wmat=-wmat

RETURN
ENDSUBROUTINE Coupling_Matrix
