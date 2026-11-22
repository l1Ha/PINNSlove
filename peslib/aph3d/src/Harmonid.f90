SUBROUTINE Harmonid (norder, hn, dhndx, x, fnorm)
!

!  P U R P O S E   O F    S U B R O U T I N E
!  This routine calculates normalized harmonic oscillator basis
!  functions and their derivative at x.


!  I N P U T    A R G U M E N T S
!  norder   maximum order of the harmonic oscillator desired.
!  x        argument of the harmonic oscillator.

!  O U T P U T    A R G U M E N T S
!  hn       array containing the harmonic oscillators evaluated at x.
!  dhndx    array of first derivatives of the hn.

!  <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

USE Numeric_Kinds_Module
USE Parms_Module
USE FileUnits_Module
USE Numbers_Module
IMPLICIT NONE
INTEGER  n, norder, nsave
REAL(Kind=WP_Kind) x, prefactr, twon, pi2, nfac, exp
REAL(Kind=WP_Kind) hn(0:norder), dhndx(0:norder), fnorm(0:norder)
SAVE prefactr, twon, pi2, nfac
EXTERNAL popt, hermited
!-----------------------------------------------------------------------
!  Determine printing options.
!-----------------------------------------------------------------------
LOGICAL little, medium, full
INTEGER ithcall, ithsub
DATA ithcall /0/, ithsub /0/
DATA little /.false./, medium /.false./, full /.false./
DATA nsave /-1/
CALL popt ('harmonid', little, medium, full, ithcall, ithsub)
!-----------------------------------------------------------------------
!  Check for input argument consistency.
!-----------------------------------------------------------------------
IF(norder < 0)THEN
   WRITE(Out_Unit,*)'Incorrect arguments to harmonid'
   WRITE(Out_Unit,*)'norder = ', norder, '  x =', x
   STOP 'harmonid'
ENDIF

!----------------------------------------------------------------------
!  calculate normalization factors.
!----------------------------------------------------------------------
!IF(norder <= nsave) GOTO 2
nsave=norder
pi2=sqrt(four*atan(one))
nfac=one
twon=one
fnorm(0)=one/sqrt(pi2)
IF(norder > 0)THEN
   DO n=1, norder
      nfac=n*nfac
      twon=two*twon
      fnorm(n)=one/sqrt(pi2*twon*nfac)
   ENDDO
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
   hn=zero
   dhndx=zero
   RETURN
ENDIF

!-----------------------------------------------------------------------
!  Evaluate the hermite polynomials.
!  note:  the derivative is that of sho fcn, not just of
!  hermite polynomial.
!-----------------------------------------------------------------------
CALL hermited (norder, hn(0), dhndx(0), x)

!-----------------------------------------------------------------------
!  multiply prefactor, norm factor and polynomial to make sho
!  functions normalized on x.
!-----------------------------------------------------------------------
DO n=0, norder
   dhndx(n)=fnorm(n)*dhndx(n)*prefactr
   hn(n)=fnorm(n)*hn(n)*prefactr
ENDDO

!-----------------------------------------------------------------------
!  IF desired WRITE out the Harmonic oscillator.
!-----------------------------------------------------------------------
IF(medium)THEN
   WRITE(Out_Unit,*)'routine Harmonic'
   WRITE(Out_Unit,*)'maximum order = ', norder
ENDIF
IF(full)THEN
   WRITE(Out_Unit,*)'Harmonic oscillators evaluated at x = ', x
   DO n = 0, norder
      WRITE(Out_Unit,*)'n= ', n, '  Hn(x) = ', hn(n)
   ENDDO
ENDIF

RETURN
ENDSUBROUTINE Harmonid
