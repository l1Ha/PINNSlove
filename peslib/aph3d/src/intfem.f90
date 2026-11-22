SUBROUTINE intfem (npoints, chi, theta, phi, th, ch, id, naph, f, numnp, nel)
USE FileUnits_Module
USE Numbers_Module
!
!
!               P U R P O S E    O F   S U B R O U T I N E
!  Interpolates the surface functions onto the new mesh.
!               I N P U T    A R G U M E N T S
!  npoints   Number of points the surface functions must be interpolated
!            at.
!  chi       Chi index for interpolation points.
!  theta     Theta index for interpolation points.
!  th        List of theta indices.
!  ch        List of chi indices.
!  id        List of node numbers.
!  naph      Number of surface functions.
!  f         Surface functions at the nodal points.
!  numnp     Number of nodal points.
!  nel       Number of elements.
!               O U T P U T    A R G U M E N T S
!  phi       Interpolated surface functions.
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
IMPLICIT NONE
INTEGER :: npoints, id, naph, numnp, nel, ipoint, ntot, iaph, k, ielem, imap
REAL(Kind=WP_Kind) :: chi, theta, th, ch, t, c, t1, t9, c1, c9, gtot, gtheta, gchi, abs
REAL(Kind=WP_Kind) :: phi, f, fi, h1, p, xj, xx, r, s, det
DIMENSION chi (npoints), theta (npoints), phi (npoints, naph), th (9, nel), ch (9, nel)
DIMENSION id (9, nel), f (numnp, naph), h1 (9), p (2, 9), imap (9), xx (2, 9), xj (2, 2)
EXTERNAL funct2
DATA imap / 2, 6, 3, 5, 9, 7, 1, 8, 4 /
WRITE(Msg_Unit, * ) 'Interpolating FEM surface functions, npoints=', npoints
ntot = 0
!-----------------------------------------------------------------------
! Loop over the points where the surface functions need to be
! determined.
!-----------------------------------------------------------------------
DO 55 ipoint = 1, npoints
!         WRITE(Out_Unit,*)'ipoint = ', ipoint, theta(ipoint), chi(ipoint)
!         WRITE(Out_Unit,*)'nel=',nel,'   naph=',naph
!-----------------------------------------------------------------------
! Loop over the elements to determine which element the point is
! within.
!-----------------------------------------------------------------------
   DO 50 ielem = 1, nel
      t = theta (ipoint)
      c = chi (ipoint)
!ccc            IF(c<0.d0)c=c+twopi
      IF(c<0.d0) c = c + pi
      t1 = th (1, ielem)
      t9 = th (9, ielem)
      c1 = ch (1, ielem)
      c9 = ch (9, ielem)
!cc            WRITE(Out_Unit,*) ielem,t1,t9,c1,c9
!-----------------------------------------------------------------------
! Test to see IF the point is within the element.
!-----------------------------------------------------------------------
      IF( (t<=t1.or.t>t9) .AND.abs ( (t - t1) / (t9 - t1) ) >1.d-4.AND.abs ( (t - t9) / (t9 - t1) ) >1.d-4) GOTO 50
      IF( (c<=c1.or.c>c9) .AND.abs ( (c - c1) / (c9 - c1) ) >1.d-4.AND.abs ( (c - c9) / (c9 - c1) ) >1.d-4) GOTO 50
!-----------------------------------------------------------------------
! This is reached IF the point is within the element hence interpolate
! the functions.
!-----------------------------------------------------------------------
  114       s = - (2 * t - t1 - t9) / (t9 - t1)
      r = (2 * c - c1 - c9) / (c9 - c1)
      CALL gtotal (t, c, gtot, gtheta, gchi)
      ntot = ntot + 1
!            WRITE(Out_Unit,*)'Found enclosing element: ntot = ',ntot
!            WRITE(Out_Unit,*) ielem,t1,t9,c1,c9
!            WRITE(Out_Unit,*)
      CALL funct2 (r, s, h1, p, xj, det, xx, 1)
!-----------------------------------------------------------------------
! Loop over all of the surface functions.
!-----------------------------------------------------------------------
      DO 30 iaph = 1, naph
         fi = 0.0
!-----------------------------------------------------------------------
! Interpolate surface function using the 9 nodes in the element.
!-----------------------------------------------------------------------
         DO 20 k = 1, 9
            fi = fi + gtot * h1 (imap (k) ) * f (id (k, ielem), iaph)
   20          ENDDO
         phi (ipoint, iaph) = fi
   30       ENDDO
!           WRITE(Out_Unit,*)ielem, ipoint
!           WRITE(Intrp2_FEM_Unit) ipoint, (phi(ipoint,iaph),iaph = 1,naph)
      GOTO 55
   50    ENDDO
   WRITE(Out_Unit, * ) 'Error Failed to find point within element.'
   WRITE(Out_Unit, * ) 'ipoint = ', ipoint, theta (ipoint) , chi ( ipoint) , c
   STOP 'intfem'
   55 ENDDO
!-----------------------------------------------------------------------
! Make sure the number of interpolation points is correct.
!-----------------------------------------------------------------------
IF(ntot/=npoints)THEN
   WRITE(Out_Unit, 81) ntot, npoints
   81 FORMAT   (1x,'***error*** ntot/=npoints',2i5)
   STOP 'interps'
ENDIF
!-----------------------------------------------------------------------
! Everything looks OK RETURN.
!-----------------------------------------------------------------------
WRITE(Msg_Unit, * ) 'Finished Interpolation'
RETURN
ENDSUBROUTINE intfem
