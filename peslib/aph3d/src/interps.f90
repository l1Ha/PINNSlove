SUBROUTINE interps (npoints, IntChi, IntTheta, phi, th, ch, id, naph, f, numnp, nel)
USE Numeric_Kinds_Module
USE FileUnits_Module
!
!   $RCSfile: interps.f,v $ $Revision: 1.15 $
!   $Date: 89/11/16 11:03:01 $
!   $State: Stable $
!

!               P U R P O S E    O F   S U B R O U T I N E

!  Interpolates the surface functions onto the new mesh.

!               I N P U T    A R G U M E N T S

!  npoints   Number of points the surface functions must be interpolated
!            at.
!  IntChi       IntChi index for interpolation points.
!  IntTheta     IntTheta index for interpolation points.
!  th        List of IntTheta indices.
!  ch        List of IntChi indices.
!  id        List of node numbers.
!  naph      Number of surface functions.
!  f         Surface functions at the nodal points.
!  numnp     Number of nodal points.
!  nel       Number of elements.

!               O U T P U T    A R G U M E N T S

!  phi       Interpolated surface functions.

IMPLICIT NONE
INTEGER npoints, naph, numnp, nel, ipoint
INTEGER ntot, c, c1, c9, t, t1, t9, iaph, k, ielem, node1, node9
INTEGER id(9,nel), imap(9), IntChi(npoints), IntTheta(npoints), th(numnp), ch(numnp) 
REAL(Kind=WP_Kind) fi, r, s, det
REAL(Kind=WP_Kind) phi(naph)
REAL(Kind=WP_Kind) f(numnp,naph), h1(9), p(2,9), xx(2,9), xj(2,2)

EXTERNAL funct2
DATA imap / 2, 6, 3, 5, 9, 7, 1, 8, 4 /

REWIND Intrp2_FEM_Unit
ntot = 0
!-----------------------------------------------------------------------
! Loop over the points where the surface functions need to be
! determined.
!-----------------------------------------------------------------------
DO ipoint = 1,npoints
!         WRITE(Out_unit,*)'ipoint = ', ipoint, IntTheta(ipoint), IntChi(ipoint)
!-----------------------------------------------------------------------
! Loop over the elements to determine which element the point is
! within.
!-----------------------------------------------------------------------
   DO 50 ielem = 1,nel
      t = IntTheta(ipoint)
      c = IntChi(ipoint)
      node1 = id(1,ielem)
      node9 = id(9,ielem)
      t1 = th(node1)
      t9 = th(node9)
      c1 = ch(node1)
      c9 = ch(node9)
!-----------------------------------------------------------------------
! Test to see IF the point is within the element.
!-----------------------------------------------------------------------
      IF(t==1.AND.t1==1) GOTO 109
      IF(t<=t1.or.t>t9) GOTO 50
  109 IF(c==1.AND.c1==1) GOTO 114
      IF(c<=c1.or.c>c9) GOTO 50
!-----------------------------------------------------------------------
! This is reached IF the point is within the element hence interpolate
! the functions.
!-----------------------------------------------------------------------
  114 s = -REAL(2*t-t1-t9,WP_Kind)/REAL(t9-t1,WP_Kind)
      r =  REAL(2*c-c1-c9,WP_Kind)/REAL(c9-c1,WP_Kind)
      ntot = ntot+1
!            WRITE(Out_unit,*)'ntot = ',ntot
      CALL funct2 (r, s, h1, p, xj, det, xx, 1)
!-----------------------------------------------------------------------
! Loop over all of the surface functions.
!-----------------------------------------------------------------------
      DO iaph = 1,naph
         fi = 0.0
         DO k = 1,9  ! Interpolate surface function using the 9 nodes in the element.
            fi = fi+h1(imap(k))*f(id(k,ielem),iaph)
         ENDDO
         phi(iaph) = fi
      ENDDO
!            WRITE(Out_unit,*)ielem, ipoint
      WRITE(Intrp2_FEM_Unit) ipoint, (phi(iaph),iaph = 1,naph)
   50 CONTINUE
ENDDO
!-----------------------------------------------------------------------
! Make sure the number of interpolation points is correct.
!-----------------------------------------------------------------------
IF(ntot/=npoints)THEN
   WRITE(Out_unit,81)ntot, npoints
   81    FORMAT(1x,'***error*** ntot/=npoints',2i5)
   STOP 'interps'
ENDIF
!-----------------------------------------------------------------------
! Everything looks OK RETURN.
!-----------------------------------------------------------------------
RETURN
END
