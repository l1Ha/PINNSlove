SUBROUTINE matgenph (u, z, w, n, rstart, rend, nstpl, xksq)
USE FileUnits_Module
USE region_MODULE   ! PORT (2026-10-02): for iregion gating of the storage writes
USE tlogd_Module
USE Storage_Module
USE engpro_Module
USE Masses_Module
USE Numbers_Module
USE dip_Module
IMPLICIT NONE
LOGICAL nsave, wrpotmat
INTEGER n, nstpl, nstep, istep, i
REAL(Kind=WP_Kind) rstart, rend, r, hstep, d1, d2, d4, drnow, wavel, h, d3
REAL(Kind=WP_Kind) u(*), z(n,n), w(n,n), xksq(n)   ! PORT: u assumed-size
REAL(Kind=WP_Kind) kmax2   ! PORT 2026-10-03: fallback wavenumber scale
EXTERNAL potmatph, smxinv
!            I N T R I N S I C    F U N C T I O N S
intrinsic sqrt
DATA hstep/1.d-4/

nsave=numder
numder=.true.
wrpotmat=.false.

!-----------------------------------------------------------------------
! Read in energy independent information if possible.
!-----------------------------------------------------------------------
wavel=2.d0*pi/sqrt(usys2*engmax)
! PORT (2026-10-03): engmax is the deck's TOTAL energy, which is negative in
! the absolute-reference convention (PORT_NOTES 5aa).  sqrt of a negative
! number gives NaN, and the step count below then becomes garbage - the
! Delves propagation ran without ever advancing the R-matrix.  When engmax
! is not a positive kinetic energy, take the fastest oscillation from the
! channel wavenumbers (xksq = 2 mu (Etot - E_channel)) instead.  For decks
! in the old positive-energy convention this branch is never taken, so the
! behaviour there is unchanged.
IF(.NOT.(wavel>0.d0))THEN
   kmax2=0.d0
   DO i=1,n
      IF(xksq(i)>kmax2) kmax2=xksq(i)
   ENDDO
   IF(kmax2>0.d0)THEN
      wavel=2.d0*pi/sqrt(kmax2)
   ELSE
      wavel=2.d0
   ENDIF
   WRITE(Out_Unit,'(A,ES12.4,A,ES12.4,A,ES12.4)') &
      ' matgenph: engmax=',engmax,' kmax2=',kmax2,' wavel=',wavel
ENDIF
IF(wavel>2.d0) wavel=2.d0
nstep=nstpl*(rend-rstart)/wavel+.9d0
IF(nstep==0) nstep=1
drnow=(rend-rstart)/nstep
h=drnow*0.5d0
d1=h*h/3.d0
d2=2.d0*d1
d4=-d1/16.d0
d3=2.d0*d2
IF(iwrindep .AND. iregion /= 'delves ')THEN   ! PORT: Delves records must not enter the APH storage
   WRITE(frst_unit)wavel,nstep,drnow,h,d1,d2,d4
   IF(lphoto)THEN
      WRITE(fph_unit)wavel,nstep,drnow,h,d1,d2,d3
   ENDIF
ENDIF

!-----------------------------------------------------------------------
! Calculate the interaction matrix.
!-----------------------------------------------------------------------
r=rstart
IF(wrpotmat)THEN
    CALL potmatph (n, r, u, hstep, -1)
ELSE
    IF(iregion /= 'delves ') WRITE(frst_unit)n,r,hstep
ENDIF
!-----------------------------------------------------------------------
! start integrating the differential equations.
!-----------------------------------------------------------------------
DO istep = 1, nstep
   r=rstart+h*(2*istep-1)
   IF(wrpotmat)THEN
      CALL potmatph (n,r,u,hstep,-1)
   ELSE
      IF(iregion /= 'delves ') WRITE(frst_unit)n,r,hstep
   ENDIF
   r=rstart+h*(2*istep)
   IF(wrpotmat)THEN
      CALL potmatph (n, r, u, hstep, -1)
   ELSE
      IF(iregion /= 'delves ') WRITE(frst_unit)n,r,hstep
   ENDIF
   IF(istep==nstep) d2=d1
ENDDO

numder=nsave

RETURN
END
