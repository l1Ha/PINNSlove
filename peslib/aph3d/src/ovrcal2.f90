SUBROUTINE ovrcal2 (nunitb, nunita, nsfunc, ntheta, nthetap,nchi, nchip, nrho, rholast, rho)
USE FileUnits_Module
USE DVR_Module
USE Prvios_Module
USE DVR2_Module
USE Parms_Module
USE TwoD_Module
IMPLICIT NONE
LOGICAL, PARAMETER:: Medium=.True.
INTEGER NSfunc, ntheta, nthetap, nchi, nchip, ialphap, l, ip, i, m, nunita, nrha, nsfa
INTEGER ntha, j, k, it, nunitb, ialpha, ncha, nrhb, nsfb, nthb, ncombp, itp
INTEGER nrho, nchb, icombp
REAL(Kind=WP_Kind) suml, summ, rhoa, sum, rhob, rholast, rho
!
!  this calculates the overlap matrix elements

!
!     csurf1,csurf2--surface function coefficients in the ray eigenvector
!                    basis
!     ovlp--overlap matrix
!
!COMMON/twod/ csurf1(ndvptmax,nsfnmax), csurf2(ndvptmax,nsfnmax), ovlp(nsfnmax,nsfnmax)
!
!     tchic,tthetac--dvr transformation matrices which are the product of
!                    the previous and current dvr transformation matrices
!
REAL(Kind=WP_Kind) tchic(nchimax,nchimax), tthetac(nthetmax,nthetmax)
!
! ----------------------------------------------------------------------------
!
IF(Ntheta>NthetMax.or.Nchi>NchiMax)THEN
   WRITE(*,*)"Ntheta or Nchi exceeds maximum allowed values"
   WRITE(*,*)"Ntheta=",Ntheta,"  NthetMax=",NthetMax
   WRITE(*,*)"Nchi=",Nchi,"  NchiMax=",NchiMax
   STOP "Ovrcal2"
ENDIF
IF(nsfunc > nsfnmax)THEN
  WRITE(Out_Unit,*) 'nsfunc greater than nsfnmax in ovrcalc'
  WRITE(Out_Unit,*) 'nsfunc = ',nsfunc,' nsfnmax = ',nsfnmax
  STOP
ENDIF
!
!  --------------------------------------------------------
! overlap matrix between surface functions at the current and
! the previous rho is calculated now
!  -------------------------------------------------------
!
! Calculate overlaps matrices for the case where the DVR orders
!  are NOT indentical for the two values of rho
!
! Calculate the composite transformation matrix for theta and chi.
!  (Notice that the tchi and ttheta matrices are being used in
!   their TRANSPOSED form: i.e., first index labels the angle,
!   second index labels the basis function).
!
WRITE(Out_Unit,*) 'ntheta,nthetap = ',ntheta,nthetap
WRITE(Out_Unit,*) 'nchi,nchip = ',nchi,nchip
!
DO ialphap=1,nthetap
   DO ialpha=1,ntheta
      suml = 0.0d0
      DO l=1,min(nthetap, ntheta)
         suml = suml + tthetap(ialphap,l)*ttheta(ialpha,l)
      ENDDO
      tthetac(ialphap,ialpha) = suml
   ENDDO
ENDDO
!     WRITE(Out_Unit,*) 'tthetac matrix'
!     CALL mxoutd2 (Out_Unit,nthetmax,tthetac,nthetap,ntheta,5)
!
DO ip=1,nchip
   DO i=1,nchi
      summ = 0.0d0
      DO m=1,min(nchip,nchi)
         summ = summ + tchip(ip,m)*tchi(i,m)
      ENDDO
      tchic(ip,i) = summ
    ENDDO
 ENDDO
!     WRITE(Out_Unit,*) 'tchic matrix'
!     CALL mxoutd2 (Out_Unit,nchimax,tchic,nchip,nchi,5)
!
! Calculate the composite indices for the first index of the
!  csurfx (x =1,2) vectors.
!
!     ndum = 0
!     DO 3040 ialphap=1,nthetap
!     DO 3040 ip=1,nchip
!       ndum = ndum + 1
!       ncmpp(ialphap,ip) = ndum
!3040 CONTINUE
!
!     ndum = 0
!     DO 3050 ialpha=1,ntheta
!     DO 3050 i=1,nchi
!       ndum = ndum + 1
!       ncmp(ialpha,i) = ndum
!3050 CONTINUE
!
! Calculate the overlap matrix.
!
! READ the second surface function into csurf2
REWIND nunita
READ(nunita) nrha,rhoa,nsfa,ntha,ncha
DO j=1,nsfunc
   READ(nunita) (csurf2(k,j), k=1,ntheta*nchi)
ENDDO

IF(Medium)THEN
   WRITE(Out_Unit,'(/2x,"surface function csurf1 at old dvr point"/)')
   CALL MxOutf(csurf2, ntheta*nchi, nsfunc, ntheta*nchi, nsfunc)       !csurf2 is actually csurf1 here
ENDIF

!
! Use csurf1 as temporary storage
DO it=1,nsfunc
   DO ialphap=1,nthetap
      DO ip=1,nchip
         !     DO 3110 icombp=1,ncombp
         !     ialphap=(icombp-1)/nchip + 1
         !     ip=icombp - (ialphap - 1)*nchip
         !
         sum = 0.0d0
         DO ialpha=1,ntheta
            DO i=1,nchi
               sum = sum + tchic(ip,i)*tthetac(ialphap,ialpha)* csurf2( ((ialpha - 1)*nchi + i ) , it ) !* csurf2(ncmp(ialpha,i), it)
            ENDDO
         ENDDO
         csurf1( ((ialphap - 1)*nchip + ip ) , it ) = sum           ! csurf1 at new dvr points
      ENDDO
   ENDDO
ENDDO

IF(Medium)THEN
   WRITE(Out_Unit,'(/2x,"surface function csurf1 at new dvr points"/)')
   CALL MxOutf(csurf1, nthetap*nchip, nsfunc, nthetap*nchip, nsfunc)
ENDIF

! READ the first surface function into csurf2
REWIND nunitb
READ(nunitb) nrhb,rhob,nsfb,nthb,nchb
DO i=1,nsfunc
   READ(nunitb) (csurf2(k,i), k=1,nthetap*nchip)
ENDDO

IF(Medium)THEN
   WRITE(Out_Unit,'(/2x,"surface function csurf2 at calculated dvr points"/)')
   CALL MxOutf(csurf2, nthetap*nchip, nsfunc, nthetap*nchip, nsfunc)
ENDIF

ncombp=nthetap*nchip
DO itp=1,nsfunc
   DO it=1,nsfunc
      sum = 0.0d0
      DO icombp=1,ncombp
         sum = sum + csurf1( icombp , it )* csurf2( icombp , itp )                               !* csurf1(ncmpp(ialphap,ip),itp)
      ENDDO
      ovlp(itp,it) = sum
    ENDDO
 ENDDO
!
! ---------------------------------------------------------------------
!
WRITE(Out_Unit,'(/2x,"Ovrcal2: surface function overlap matrix elements"/)')

WRITE(dvr18) nrho-1,rholast,rho
WRITE(dvr18) nsfunc,rholast,rho
DO j=1,nsfunc
   WRITE(dvr18) (ovlp(i,j), i=1,nsfunc)
ENDDO

IF(Medium)THEN
   WRITE(Out_Unit,'(/2x,"Ovrcal2: surface function overlap matrix elements"/)')
   CALL MxOutf(ovlp, nsfunc, nsfunc, nsfnmax, nsfnmax)
ENDIF

RETURN
END
