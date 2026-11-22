SUBROUTINE Sfpmatrix(jtot,jrot,lorb,nvib, nchan,bfpotmat,sfpotmat,jmax)
USE Numeric_Kinds_Module
USE numbers_Module
!      USE popt_Module

! Transforms from body-fixed matrix elements to space-fixed
! matrix elements and the results are stored in sfpotmat.

IMPLICIT NONE
INTEGER nchan, f, n, vf, jf, lf, vn, jn, ln, jmax, abslam
INTEGER jtot, jrot(nchan), lorb(nchan), nvib(nchan), lambda
CHARACTER(LEN=21), PARAMETER:: ProcName='sfpmatrix'
REAL(Kind=WP_Kind) bfpotmat(nchan,nchan,0:jmax)
REAL(Kind=WP_Kind) sfpotmat(nchan,nchan), cleb, cleb1, cleb2, prefac
LOGICAL little, medium, full
INTEGER ithcll, ithsub
data ithcll/0/,ithsub/0/
data little/.true./, medium/.false./, full/.false./


CALL popt(procname, little, medium, full, ithcll, ithsub)

! This corresponds to rtp98/4/30-5 Eq. 26.
! The body-fixed potential energy matrix elements stored in
! array bfpotmat are transformed into the space-fixed
! representation using the clebsch-gordan transformation.

DO f=1,nchan
   vf=nvib(f)
   jf=jrot(f)
   lf=lorb(f)
   DO n=1,nchan
      vn=nvib(n)
      jn=jrot(n)
      ln=lorb(n)
      prefac=sqrt(REAL(2*lf+1,WP_Kind)*REAL(2*ln+1,WP_Kind))/REAL(2*jtot+1,WP_Kind)
!  note:  bfpotmat is only calculated for lambda>=0.  Its values
!  for -lambda is the same as for lambda.  The following runs over
!  all lambda but uses this fact.
!  The Clebsch-Gordan coeffs do depend on the sign of lambda.
      sfpotmat(f,n)=0.d0
      DO lambda=-min(jf,jn,jtot),min(jf,jn,jtot)
         abslam=abs(lambda)
         cleb1=cleb(jf,lf,jtot,lambda,0,lambda)
         cleb2=cleb(jn,ln,jtot,lambda,0,lambda)
         sfpotmat(f,n)=sfpotmat(f,n)+prefac*cleb1*cleb2*bfpotmat(f,n,abslam)
      ENDDO
   ENDDO
ENDDO

! On return sfpotmat contains the space-fixed potential energy
! matrix elements.

RETURN
ENDSUBROUTINE Sfpmatrix
