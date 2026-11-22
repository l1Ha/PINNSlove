SUBROUTINE interpol (f1, th1, ch1, id1, nel1, numnp1, th2, ch2, id2, nel2, &
   numnp2, naph, ntheta, nchi, igrid, f3, f4, IntChi, IntTheta, phi)
USE Numeric_Kinds_Module
USE FileUnits_Module
!
!   $RCSfile: interpol.f,v $ $Revision: 1.15 $
!   $Date: 89/11/16 11:03:00 $
!   $State: Stable $
!

!               P U R P O S E    O F   S U B R O U T I N E
!               I N P U T    A R G U M E N T S

!  f1        Surface functions at rho1.
!  th1       List of theta indices at rho1.
!  ch1       List of chi indeices at rho1.
!  id1       List of node numbers at rho1.
!  nel1      Number of elements at rho1.
!  numnp1    Number of nodal points at rho1.
!  th2       List of theta indices at rho2.
!  ch2       List of chi indices at rho2.
!  id2       List of node numbers at rho2.
!  nel2      Number of elements at rho2.
!  numnp2    Number of nodal points at rho2.
!  naph      Number of surface functions.
!  ntheta    Number of theta angles.
!  nchi      Number of chi angles.
!  igrid     Temporary array to determine which points are not common
!            to both meshes.
!  f3        Temporary array.
!  f4        Should be equivalenced to f1.
!  IntChi      Temporary array.
!  IntTheta   Temporary array.
!  phi       Temporary array.
!
!               O U T P U T    A R G U M E N T S
!
!  f4        Contains the interpolated values of the surface functions
!            at rho1 on the rho2 mesh.

IMPLICIT NONE

!      I N T E G E R S
INTEGER id1, id2, igrid, naph, nel1, nel2, ntheta, nchi, numnp1, numnp2, itheta, jchi, iel, k, node, inode, iaph
INTEGER k2, iel2, inode1, inode2, icount, ipoint, npoints, ch1, ch2, IntChi, th1, th2, IntTheta

!      R E A L S
REAL(Kind=WP_Kind) phi, f1, f3, f4

!      D I M E N S I O N S
DIMENSION IntChi(numnp2), IntTheta(numnp2), phi(naph), th1(numnp1),  ch1(numnp1), id1(9,nel1), th2(numnp2), ch2(numnp2)
DIMENSION id2(9,nel2), f1(numnp1,naph),  f3(numnp2), f4(numnp2,naph), igrid(2,nchi,ntheta)

!      E X T E R N A L S
EXTERNAL interps
!----------------------------------------------------------------------
! NOTE: In the program that calls this routine f4 should be
!       equivalenced to f1.
! Initialize igrid to a negative one.
!----------------------------------------------------------------------
DO 15 itheta = 1,ntheta
   DO 10 jchi = 1,nchi
      igrid(1,jchi,itheta) = -1
      igrid(2,jchi,itheta) = -1
   10    CONTINUE
   15 CONTINUE
!-----------------------------------------------------------------------
! Store element indices in igrid for the second rho value.
!-----------------------------------------------------------------------
DO 25 iel = 1,nel2
   DO 20 k = 1,9
      node = id2(k,iel)
      itheta = th2(node)
      jchi = ch2(node)
      igrid(1,jchi,itheta) = k
      igrid(2,jchi,itheta) = iel
   20    CONTINUE
   25 CONTINUE
!----------------------------------------------------------------------
! Initialize f3 to a very large negative number
!----------------------------------------------------------------------
DO 35 iel = 1,nel2
   DO 30 k = 1,9
      inode = id2(k,iel)
      f3(inode) = -1.0d+30
   30    CONTINUE
   35 CONTINUE
!----------------------------------------------------------------------
! For the nodes that are common in both grids set f3 to the
! value of the surface function at the first rho value, but indexed by
! the node number at the second rho value.
!----------------------------------------------------------------------
REWIND Intrp1_FEM_Unit
DO iaph = 1,naph
   DO iel = 1,nel1
      DO k = 1,9
         node = id1(k,iel)
         itheta = th1(node)
         jchi = ch1(node)
         IF(igrid(1,jchi,itheta)>0)THEN
            k2 = igrid(1,jchi,itheta)
            iel2 = igrid(2,jchi,itheta)
            inode2 = id2(k2,iel2)
            inode1 = id1(k,iel)
            f3(inode2) = f1(inode1,iaph)
         ENDIF
      ENDDO
   ENDDO
   WRITE(Intrp1_FEM_Unit) (f3(inode),inode = 1,numnp2)
ENDDO
!----------------------------------------------------------------------
! Find out the number of points that are in the new grid that were not
! in the old grid.
!----------------------------------------------------------------------
icount = 0
DO 60 inode = 1,numnp2
   IF(f3(inode)==-1.0d+30) icount = icount+1
   60 CONTINUE
WRITE(Out_unit,130) icount
!----------------------------------------------------------------------
! IF there are no new points RETURN after loading f4 with the values of
! the first rho surface functions.
!----------------------------------------------------------------------
IF(icount==0)THEN
REWIND Intrp1_FEM_Unit
DO 75 iaph = 1,naph
   READ(Intrp1_FEM_Unit) (f3(inode),inode = 1,numnp2)
   DO 70 inode = 1,numnp2
      f4(inode,iaph) = f3(inode)
   70    CONTINUE
   75 CONTINUE
RETURN
ENDIF
!----------------------------------------------------------------------
! Interpolate to find the functions at the new nodal points.
!----------------------------------------------------------------------
npoints = 0
DO 85 iel = 1,nel2
   DO 80 k = 1,9
      IF(f3(id2(k,iel))==-1.0d+30)THEN
         f3(id2(k,iel)) = -1.0d+29
         npoints = npoints+1
         inode2 = id2(k,iel)
         IntTheta(npoints) = th2(inode2)
         IntChi(npoints) = ch2(inode2)
      ENDIF
   80    CONTINUE
   85 CONTINUE
CALL interps (npoints, IntChi, IntTheta, phi, th1, ch1, id1, naph, f1, numnp1, nel1)
REWIND Intrp1_FEM_Unit
DO 95 iaph = 1,naph
   READ(Intrp1_FEM_Unit) (f3(inode),inode = 1,numnp2)
   DO 90 inode = 1,numnp2
      f4(inode,iaph) = f3(inode)
   90    CONTINUE
   95 CONTINUE
REWIND Intrp2_FEM_Unit
npoints = 0
DO 115 iel = 1,nel2
   DO 110 k = 1,9
   IF(f3(id2(k,iel))==-1.0d+30)THEN
      npoints = npoints+1
      READ(Intrp2_FEM_Unit) ipoint,(phi(iaph),iaph = 1,naph)
      IF(ipoint/=npoints)THEN
         WRITE(Out_unit,140) ipoint,npoints
         STOP 'interpol'
      ENDIF
      DO 100 iaph = 1,naph
         f4(id2(k,iel),iaph) = phi(iaph)
  100       CONTINUE
      f3(id2(k,iel)) = phi(1)
   ENDIF
  110    CONTINUE
  115 CONTINUE
icount = 0
DO 120 inode = 1,numnp2
   IF(f3(inode)==-1.0d+30) icount = icount+1
  120 CONTINUE
IF(icount==0) RETURN
WRITE(Out_unit,150)
STOP 'interpol'
!
  130 FORMAT(' there were ',i5,' nodes in the new grid that were not in  ',' the old grid')
  140 FORMAT(/,'***error*** ipoint/=npoints',2i5)
  150 FORMAT(/,'***error***: icount is not equal to zero.')
END
