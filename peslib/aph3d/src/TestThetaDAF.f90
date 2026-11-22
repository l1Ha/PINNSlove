SUBROUTINE TestThetaDAF (ntheta, kinetic, theta, DafTheta, nd_daf)
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
INTEGER, PARAMETER::nleg = 8
INTEGER :: ntheta, itheta, ktheta, jleg, nd_daf
REAL(Kind=WP_Kind) :: x
REAL(Kind=WP_Kind) :: kinetic (ntheta, ntheta), theta (ntheta)
REAL(Kind=WP_Kind) :: DafTheta ( - nd_daf:nd_daf, 0:2)
REAL(Kind=WP_Kind), ALLOCATABLE::eig (:), Lmat (:, :)
REAL(Kind=WP_Kind), ALLOCATABLE::pn (:, :), pnp (:, :), pnpp (:, :)

REAL(Kind=WP_Kind), ALLOCATABLE::pmat (:, :), daf0 (:, :), daf1 (:, :), daf2 (:, :)
ALLOCATE (eig (0:nleg), pn (ntheta, 0:nleg), pnp (ntheta, 0:nleg), &
 pnpp (ntheta, 0:nleg), pmat (ntheta, 0:nleg), Lmat (0:nleg, 0: &
 nleg) )

ALLOCATE (daf0 (ntheta, ntheta), daf1 (ntheta, ntheta), daf2 ( &
 ntheta, ntheta) )
WRITE(Out_Unit, * )

WRITE(Out_Unit, * ) 'Begin(TestThetaDAF)'
WRITE(Out_Unit, * )
WRITE(Out_Unit,  * ) '      Theta                    Legendre Polynomials'
DO itheta = 1, ntheta
x = cos (2.d0 * theta (itheta) )
pn (itheta, 0) = 1.d0
pn (itheta, 1) = x
DO jleg = 1, nleg - 1
pn (itheta, jleg + 1) = 2.d0 * x * pn (itheta, jleg) - pn (itheta,jleg - 1) &
   - (x * pn (itheta, jleg) - pn (itheta, jleg - 1) )/ (jleg +1)
ENDDO
WRITE(Out_Unit, 50) theta (itheta), (pn (itheta, jleg), jleg = 0, NLeg)

ENDDO
WRITE(Out_Unit, * )
WRITE(Out_Unit,  * ) '      Theta                 D(Legendre Polynomials)'
DO itheta = 1, ntheta
x = cos (2.d0 * theta (itheta) )
pnp (itheta, 0) = 0.d0
DO jleg = 1, nleg
pnp (itheta, jleg) = - jleg * (x * pn (itheta, jleg) - &
                     pn (itheta,jleg - 1) ) / (1 - x**2)
ENDDO
WRITE(Out_Unit, 50) theta (itheta), (pnp (itheta, jleg), jleg = 0, NLeg)

ENDDO
Lmat = 0.d0
DO jleg = 0, nleg
Lmat (jleg, jleg) = jleg * (jleg + 1)

ENDDO
DO itheta = 1, Min (ntheta, nd_daf)
DO ktheta = 1, Min (ntheta, nd_daf)
daf0 (itheta, ktheta) = DafTheta (itheta - ktheta, 0)
daf1 (itheta, ktheta) = DafTheta (itheta - ktheta, 1)
daf2 (itheta, ktheta) = DafTheta (itheta - ktheta, 2)
ENDDO


ENDDO

   50 FORMAT(1x,11es12.4)
WRITE(Out_Unit, * )
WRITE(Out_Unit, * ) 'Daf(0)'

CALL MxOut(Daf0, NTheta, NTheta)
WRITE(Out_Unit, * )
WRITE(Out_Unit, * ) 'Daf(1)'

CALL MxOut(Daf1, NTheta, NTheta)
WRITE(Out_Unit, * )
WRITE(Out_Unit, * ) 'Daf(2)'

CALL MxOut(Daf2, NTheta, NTheta)
WRITE(Out_Unit, * )
WRITE(Out_Unit, * ) 'Kinetic'

CALL MxOut(Kinetic, NTheta, NTheta)
      !Test Pn-Daf(0)*Pn
WRITE(Out_Unit, * )
WRITE(Out_Unit, * ) 'Test Pn-Daf(0)*Pn'
Pmat = Pn - Matmul (Daf0, Pn)

CALL MxOut(Pmat, NTheta, NLeg + 1)
      !Test Pnp-Daf(1)*Pn
WRITE(Out_Unit, * )
WRITE(Out_Unit, * ) 'Test Pnp-Daf(1)*Pn'
Pmat = Pnp - Matmul (Daf1, Pn)

CALL MxOut(Pmat, NTheta, NLeg + 1)
      !Test L*(L+1)*Pn-Daf(1)*Sin[2*Theta]*Daf(1)*Pn
WRITE(Out_Unit, * )
WRITE(Out_Unit, * ) 'Test L*(L+1)*Pn-Daf(1)*Sin[2*Theta]*Daf(1)*Pn'
Pmat = Matmul (Pn, Lmat) - Matmul (Daf1, Pn)

CALL MxOut(Pmat, NTheta, NLeg + 1)
      !Test Number L*(L+1)*Pn-Kinetic*Pn
WRITE(Out_Unit, * )
WRITE(Out_Unit, * ) 'Test L*(L+1)*Pn-Kinetic*Pn'
Pmat = Matmul (Pn, Lmat) - Matmul (Kinetic, Pn)

CALL MxOut(Pmat, NTheta, NLeg + 1)
WRITE(Out_Unit, * ) 'END(TestThetaDAF)'

WRITE(Out_Unit, * )
RETURN

ENDSUBROUTINE TestThetaDAF
