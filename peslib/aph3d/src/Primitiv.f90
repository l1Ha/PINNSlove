SUBROUTINE Primitiv (nang, nbasis, mega, oscil,  pn, parityfn, weight, phi, wptg, wpth, bfaci)
!

!            P U R P O S E   O F    S U B R O U T I N E
!   This routine sets up the primitive basis set.

!            I N P U T    A R G U M E N T S
!   nang     number of angles.
!   nbasis   number of basis functions.
!   narran   number of arrangement channels.
!   oscil    an array containing the oscillators evaluated at each
!             of the angles.
!   pn       an array containing the Legendre polynomials evaluated
!            at each of the angles.
!   parityfn an array containing the cosine of the APH chi angle
!            for odd parity functions and 1 for even parity
!            functions.
!   weight   an array containing the weights for the Gaussian
!            quadrature.
!   wptg     Gauss_Legendre weights.
!   wpth     Gauss_Hermite weights.
!   bfaci    a factor that goes with the Jacobian and sho functions.

!            O U T P U T    A R G U M E N T S
!   phi      an array of the primitive basis evaluated at each of
!            the quadrature points.  Each column is a basis function
!            evaluated at each of the quadrature angles.

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

USE Parms_Module
USE FileUnits_Module
USE Quantb_Module
USE QCase_Module
IMPLICIT NONE

!            I N T E G E R S
INTEGER  iang, nang, ibasis, nbasis, marran, n, l, mega

!            R E A L S
REAL(Kind=WP_Kind) sqrt
REAL(Kind=WP_Kind) phi(nang, nbasis), pn(0:mxl, mxang, narran)
REAL(Kind=WP_Kind) oscil(0:mxn, mxang, narran)
REAL(Kind=WP_Kind) parityfn(mxang,narran), weight(mxang)
REAL(Kind=WP_Kind) wptg(mxglegn, narran), wpth(mxhermt, narran)
REAL(Kind=WP_Kind) bfaci(mxang,narran)

!            I N T R I N S I C    F U N C T I O N S
INTRINSIC sqrt

!-----------------------------------------------------------------------
! Determine printing options.
!-----------------------------------------------------------------------
LOGICAL little, medium, full
INTEGER ithcall, ithsub
DATA ithcall /0/, ithsub /0/
DATA little /.false./, medium /.false./, full /.false./
CALL popt ('primitiv', little, medium, full, ithcall, ithsub)

!  ----------------------------------------------------------------
!  build total primitive as a product of a legendre poly, a sho,
!  a parity function, a jacobian factor, and a sqrt of a weight.
!  loop 4 runs over the quadrature points
!  marran runs over the function arrangement channels.
!  loops 1 and 2 run over the functions.
!  ----------------------------------------------------------------
DO iang = 1, nang
   ibasis = 0
   DO marran = 1, narran
      DO n = minvib(marran), maxvib(marran)
         DO l = jmin(n,marran,mega), jmax(n,marran),  jskip(marran)
            ibasis = ibasis + 1
            phi(iang,ibasis) = pn(l, iang, marran)*oscil(n, iang, marran)*parityfn(iang,marran)*weight(iang)*bfaci(iang,marran)
         ENDDO
      ENDDO
   ENDDO
ENDDO

IF(ibasis/=nbasis)THEN
   WRITE(Out_Unit,*)'ibasis, nbasis=', ibasis, nbasis
   STOP 'primitiv'
ENDIF


!-----------------------------------------------------------------------
! IF desired print out the basis.
!-----------------------------------------------------------------------
IF(little)THEN
   WRITE(Out_Unit,*)'routine primitiv'
   WRITE(Out_Unit,*)'nang, nbasis = ',nang, nbasis
ENDIF
IF(medium.AND..NOT.full)THEN
   DO iang = 1, min(100,nang)
      WRITE(Out_Unit,*)(phi(iang,ibasis),ibasis=1,min(10,nbasis))
   ENDDO
ENDIF
IF(full)THEN
   WRITE(Out_Unit,*)'primitiv basis functions'
   CALL MxOut(phi, nang, nbasis)
ENDIF

RETURN
ENDSUBROUTINE Primitiv
