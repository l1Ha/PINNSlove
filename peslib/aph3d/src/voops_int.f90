FUNCTION voops_int(upsln1,upsln2,vleg,dt1,dt2,f)
USE Numeric_Kinds_Module
USE numbers_Module
USE fileunits_Module
USE int1d_Module
IMPLICIT NONE
CHARACTER(LEN=21), PARAMETER:: ProcName='voops_int'

! integration or quadrature routine for evaluating PES matrix
! elements. This is the little-theta integration.

INTEGER n
REAL(Kind=WP_Kind) voops_int, vleg(nltheta+1),dt1,dt2
REAL(Kind=WP_Kind) voops_int1, voops_int2
REAL(Kind=WP_Kind) upsln1(nltheta+1), upsln2(nltheta+1), f(nltheta+1)

! store integrand at quadrature points

DO n=1,nltheta+1
   f(n)=upsln1(n)*vleg(n)*upsln2(n)
!         WRITE(Out_Unit,*)' voops_int integrand=', f(n)
ENDDO
!
! trapezoidal rule implemented 98/5 by GAP who suggested improving it.
! DCM 6/18/98: Upgraded the integration method to Simpson's rule.  Note
! that n1,n2 must be even, otherwise the last panel gets counted twice!
! (Same comment goes for oops_int.f)
! revised by RTP 99/2/12 to speed it and correct index ranges.   
! The integration is divided into two segments.
! The first segment has n1 intervals. The second segment has nltheta-n1
! intervals.  Since index n starts at 1, the first segment uses points
! 1 thru n1+1, and the second uses points n1+1 thru nltheta+1.
! The boundary conditions on the upsilons require that they vanish at
! points 1 and nltheta+1, so those are omitted.  RTP 99/2/12
! Simpson's rule.  Let's try to minimize the number of multiplies
! and divides.
! first half of integration.
!
voops_int1=zero
!  add in all the interior points
DO n=2,n1
   voops_int1 = voops_int1 + f(n)
ENDDO
!  add in the more heavily weighted points again
DO n=2,n1,2
   voops_int1 = voops_int1 + f(n)
ENDDO
!  get weights up to four and two and finish the first integral
voops_int1 = dt1*(two*voops_int1 + f(n1+1))/three
!

! second half of integration.  Same pattern as first half.    
!
voops_int2 = zero
DO n=n1+2,nltheta
   voops_int2 = voops_int2 + f(n)
ENDDO
DO n=n1+2,nltheta,2
   voops_int2 = voops_int2 + f(n)
ENDDO
voops_int2 = dt2*(two*voops_int2 + f(n1+1))/three
!  add the two halves
voops_int = voops_int1 + voops_int2
!
RETURN
END
