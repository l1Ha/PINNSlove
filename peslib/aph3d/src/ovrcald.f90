SUBROUTINE ovrcald (nunit2, nunit1, nsfunc, nthnch, nrho, rholast, rho,ii, medium)
USE Numeric_Kinds_Module
USE FileUnits_Module
USE DVR_Module
USE Parms_Module
USE TwoD_Module
!
IMPLICIT NONE
INTEGER NSfunc, NRhb, Nsfb, nthb, nchb, i, k, nthnch, nunit1, nrha, ntha, ncha
INTEGER j, m, n, l, nrho, ii, nunit2, nsfa
REAL(Kind=WP_Kind) rhoa, rhob, rholast, rho
LOGICAL medium
!  this calculates the overlap matrix elements.  All the arrays needed are supplied by twod.
!
!COMMON/twod/ csurf1(ndvptmax,nsfnmax), csurf2(ndvptmax,nsfnmax), ovlp(nsfnmax,nsfnmax)
!
IF(nsfunc > nsfnmax)THEN
  WRITE(Out_Unit,*) 'nsfunc greater than nsfnmax in ovrcalc'
  WRITE(Out_Unit,*) 'nsfunc = ',nsfunc,' nsfnmax = ',nsfnmax
  STOP
ENDIF
!
!  --------------------------------------------------------
! overlap matrix between surface functions at the current and
! the previous rho is calculated now
!
! the overlap matrix will be stored in the first nsfunc by nsfunc
! elements of the hall array.
!  -------------------------------------------------------
!
!     nstart = nsfunc*nsfunc
REWIND nunit2
READ(nunit2) nrhb,rhob,nsfb,nthb,nchb
DO i=1,nsfunc
   READ(nunit2) (csurf1(k,i), k=1,nthnch)
ENDDO

IF(Medium)THEN
   WRITE(Out_Unit,'(/2x,"surface function csurf1"/)')
   CALL MxOutf(csurf1, nthnch, nsfunc, ndvptmax, nsfnmax)
ENDIF
!
REWIND nunit1
READ(nunit1) nrha,rhoa,nsfa,ntha,ncha
DO j=1,nsfunc
   READ(nunit1) (csurf2(k,j), k=1,nthnch)
ENDDO

IF(Medium)THEN
   WRITE(Out_Unit,'(/2x,"surface function csurf2"/)')
   CALL MxOutf(csurf2, nthnch, nsfunc, ndvptmax, nsfnmax)
ENDIF
!
m=nsfunc
n=nthnch    ! n=ntheta*nchi
l=nsfunc
CALL dgemm('T','N',m ,l ,n ,1.d0 ,csurf1 ,ndvptmax ,csurf2 ,ndvptmax ,0.d0 ,ovlp  ,nsfnmax)
!
!  This set is for unformatted WRITEs.
WRITE(dvr18) nrho-1,rholast,rho                    !Be Careful
WRITE(dvr18) nsfunc,rholast,rho                    !Be Careful
WRITE(Out_Unit,*)'In ovrcald: ',nrho-1,rholast,rho,nsfunc
DO j=1,nsfunc
   WRITE(dvr18) (ovlp(i,j), i=1,nsfunc)
ENDDO

IF(Medium)THEN
   WRITE(Out_Unit,'(/2x,"Ovrcald: surface function overlap matrix elements"/)')
   CALL MxOutf(ovlp, nsfunc, nsfunc, nsfnmax, nsfnmax)
ENDIF

RETURN
END
