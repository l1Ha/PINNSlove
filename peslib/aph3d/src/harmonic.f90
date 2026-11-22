SUBROUTINE harmonic (norder, hn, x)
USE Numeric_Kinds_Module

!  P U R P O S E   O F    S U B R O U T I N E
!  This routine calculates normalized harmonic oscillator basis
!  functions at x.


!  I N P U T    A R G U M E N T S
!  norder   maximum order of the harmonic oscillator desired.
!  x        argument of the harmonic oscillator.

!  O U T P U T    A R G U M E N T S
!  hn       array containing the harmonic oscillators evaluated at x.

USE FileUnits_Module
USE Narran_Module
USE Parms_Module
IMPLICIT NONE
SAVE

!  I N T E G E R S
INTEGER  n, norder, nsave

!  R E A L S
REAL(Kind=WP_Kind) x, hn(0:mxn), prefactr, twon, pi2, nfac, zero, one, half, two, four
REAL(Kind=WP_Kind), ALLOCATABLE:: fnorm(:)
EXTERNAL popt, vsets, hermite

!  P A R A M E T E R S
PARAMETER (zero=0.0d0,half=0.5d0,one=1.0d0,two=2.0d0,four=4.0d0)

!-----------------------------------------------------------------------
!  Determine printing options.
!-----------------------------------------------------------------------
LOGICAL little, medium, full
INTEGER ithcall, ithsub
DATA ithcall /0/, ithsub /0/
DATA little /.false./, medium /.false./, full /.false./
DATA nsave /-1/
CALL popt ('harmonic', little, medium, full, ithcall, ithsub)
!-----------------------------------------------------------------------
!  Check for input argument consistency.
!-----------------------------------------------------------------------
IF(norder < 0)THEN
   WRITE(Out_Unit,*)'Incorrect arguments to harmonic'
   WRITE(Out_Unit,*)'norder = ', norder, '  x =', x
   STOP 'harmonic'
ENDIF
IF(.NOT.ALLOCATED(fnorm))ALLOCATE(fnorm(0:mxn))

!----------------------------------------------------------------------
!  calculate normalization factors.
!----------------------------------------------------------------------
IF(norder <= nsave) GOTO 2
nsave=norder
pi2=sqrt(four*atan(one))
nfac=one
twon=one
fnorm(0)=one/sqrt(pi2)
IF(norder > 0)THEN
   DO 1 n=1, norder
      nfac=n*nfac
      twon=two*twon
      fnorm(n)=one/sqrt(pi2*twon*nfac)
 1       CONTINUE
ENDIF
IF(medium)WRITE(Out_Unit,*)'fnorm=', (fnorm(n), n=0, norder)
 2    CONTINUE

!-----------------------------------------------------------------------
!  Evaluate the exponential prefactor. IF the prefactor is zero THEN
!  all of the harmonic oscillators must be zero. Hence, zero the array
!  and RETURN.
!-----------------------------------------------------------------------
prefactr = exp(-half*x*x)
IF(prefactr == 0)THEN
   CALL vsets (norder+1, hn(0), 1, zero)
   RETURN
ENDIF

!-----------------------------------------------------------------------
!  Evaluate the hermite polynomials.
!-----------------------------------------------------------------------
CALL hermite (norder, hn(0), x)

!-----------------------------------------------------------------------
!  multiply prefactor, norm factor and polynomial to make sho
!  functions normalized on x.
!-----------------------------------------------------------------------
DO 3 n=0, norder
   hn(n)=fnorm(n)*hn(n)*prefactr
 3    CONTINUE

!-----------------------------------------------------------------------
!  IF desired WRITE out the Harmonic oscillator.
!-----------------------------------------------------------------------
IF(medium)THEN
   WRITE(Out_Unit,*)'routine Harmonic'
   WRITE(Out_Unit,*)'maximum order = ', norder
ENDIF
IF(full)THEN
   WRITE(Out_Unit,*)'Harmonic oscillators evaluated at x = ', x
   DO 5 n = 0, norder
      WRITE(Out_Unit,*)'n= ', n, '  Hn(x) = ', hn(n)
 5       CONTINUE
ENDIF

RETURN
END
