SUBROUTINE sorter (eig, eig_im, ntheta, psi)
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
INTEGER :: ntheta, i, k, kmin, k1, kp

REAL(Kind=WP_Kind) :: psi (ntheta, ntheta), eig (ntheta), eig_im (ntheta), esave, psave, emin
WRITE(Out_Unit, * )
WRITE(Out_Unit, * ) 'Begin(Sorter)'
!-----------------------------------------------------------------------
! order quantum numbers according to energy.
!-----------------------------------------------------------------------
DO k = 1, ntheta - 1
   emin = eig (k)
   kmin = k
   k1 = k + 1
   DO 80 kp = k1, ntheta
      IF(emin<=eig (kp) ) goto 80
      emin = eig (kp)
      kmin = kp
80 ENDDO
   esave = eig (k)
   eig (k) = eig (kmin)
   eig (kmin) = esave
   esave = eig_im (k)
   eig_im (k) = eig_im (kmin)
   eig_im (kmin) = esave
   DO i = 1, ntheta
      psave = psi (i, k)
      psi (i, k) = psi (i, kmin)
      psi (i, kmin) = psave
   ENDDO
ENDDO
WRITE(Out_Unit, * ) 'END(Sorter)'

WRITE(Out_Unit, * )
RETURN


ENDSUBROUTINE sorter
