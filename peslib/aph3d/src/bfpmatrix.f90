SUBROUTINE bfpmatrix(jrot,lorb,nvib,dt1,dt2,upsln,vleg,nchan,work,bfpotmat,jmax,jtot)
USE Numeric_Kinds_Module
USE Numbers_Module
USE int1d_Module

! Calculates body-fixed potential energy matrix elements and
! stores results in bfpotmat

IMPLICIT NONE
CHARACTER(LEN=21), PARAMETER:: ProcName='bfpmatrix'
INTEGER nchan, jmax
REAL(Kind=WP_Kind) upsln(nltheta+1,nchan)
REAL(Kind=WP_Kind) work(nltheta+1), dt1, dt2, bfpotmat(nchan,nchan,0:jmax)
REAL(Kind=WP_Kind) voops_int
REAL(Kind=WP_Kind) vleg(nltheta+1,0:jmax,0:jmax,0:jmax)
INTEGER f, n, vf, jf, lf, vn, jn, ln, jtot
INTEGER jrot(nchan), lorb(nchan), nvib(nchan), lambda
LOGICAL little, medium, full
INTEGER ithcll, ithsub
data ithcll/0/,ithsub/0/
data little/.true./, medium/.false./, full/.false./


CALL popt(procname, little, medium, full, ithcll, ithsub)

! This corresponds to rtp98/4/30-5 Eq. (26).  bfpotmat has both
! space-fixed labels f=(v,j,l) and n and the BF label Lambda.
! The Big-theta integrals have already been done and are stored
! in the array vleg. This array is a function of the little-theta
! coordinate. This routine uses the upsilons that have already
! been calculated and integrates the results using the
! voops_int routine.

DO f=1,nchan
   vf=nvib(f)
   jf=jrot(f)
   lf=lorb(f)
   DO n=1,f
      vn=nvib(n)
      jn=jrot(n)
      ln=lorb(n)
      DO lambda=0,min(jf,jn,jtot)
         bfpotmat(f,n,lambda)=voops_int(upsln(1,f),upsln(1,n),  &
         vleg(1,jf,jn,lambda),dt1,dt2,work)
         bfpotmat(n,f,lambda)=bfpotmat(f,n,lambda)
      ENDDO
   ENDDO
ENDDO

! On return from this routine the array bfpotmat contains
! the "body-fixed" potential energy matrix elements as
! written in rtp98/4/30-5 Eq. (26).

RETURN
END
