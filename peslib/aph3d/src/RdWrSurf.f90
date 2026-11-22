SUBROUTINE rdwrsurf (id, th, thetaval, ch, chivals, phi, energy,  nelem, &
  numnp, naph, ntheta, nchi, iunit, kase, mxnode, mxelem, lambda, lammin)
USE rhosur_Module
USE rhochk_Module
USE FileUnits_Module
!
!   $RCSfile: rdwrsurf.f,v $ $Revision: 1.15 $
!   $Date: 89/11/16 11:03:23 $
!   $State: Stable $
!

!               P U R P O S E    O F    S U B R O U T I N E

!  Reads or WRITEs surface function DATA.

!               I N P U T    A R G U M E N T S

!  lambda   Project of total angular momentum on the Body-fixed z-axis.
!  lammin   Minimum lambda allowed.
!  kase     IF kase is equal to 1 READ surface function DATA.
!           IF kase is not equal to 1 WRITE surface function DATA.
!  mxnode   Maximum number of nodes.
!  mxelem   Maximum number of elements.
!  naph     Number of surface functions.
!  ntheta   Number of theta angles
!  nchi     Number of chi angles
!  iunit    unit number for surface functions.

!               O U T P U T    A R G U M E N T S

!  id       List of node numbers
!  th       List of theta indices.
!  thetaval List of theta values.
!  ch       List of chi indices.
!  chivals  List of chi values.
!  phi      Surface functions.
!  energy   Surface function eigenenergies.
!  nelem    Number of elements.
!  numnp    Number of nodal points.

IMPLICIT NONE
INTEGER id, th, ch, nelem, numnp, naph, lambda, k, iel, node, kase, ntheta, nchi, lammin, index, lam, mxnode, mxelem, iunit
REAL(Kind=WP_Kind) thetaval, chivals, phi, energy
DIMENSION id(9,mxelem), th(mxnode), ch(mxnode), thetaval(ntheta), chivals(nchi), phi(mxnode*naph), energy(naph)

IF(kase==1)THEN
!-----------------------------------------------------------------------
! Position the Disk head to the correct position.
!-----------------------------------------------------------------------
WRITE(Out_unit,*) lammin, lambda
   REWIND iunit
   IF(lambda>lammin)THEN
      DO 1 lam = lammin,lambda-1
         READ(iunit)
         READ(iunit)
         READ(iunit)
         READ(iunit)
         READ(iunit)
         READ(iunit)
         READ(iunit)
         READ(iunit)
    1       CONTINUE
   ENDIF
!-----------------------------------------------------------------------
! READ in surface function DATA
!-----------------------------------------------------------------------
   READ(iunit) lam, nelem, numnp
       WRITE(Out_unit,*)'mxelem,nelem=',mxelem,nelem
   READ(iunit) ((id(k,iel),k=1,9),iel=1,nelem)
    WRITE(Out_unit,*)'mxnode,numnp=',mxnode,numnp
   READ(iunit) (th(node),node=1,numnp)
   READ(iunit) thetaval
   READ(iunit) (ch(node),node=1,numnp)
   READ(iunit) chivals
   READ(iunit) (phi(index),index=1,numnp*naph)
   READ(iunit) energy,rhor
!-----------------------------------------------------------------------
! WRITE surface function DATA.
!-----------------------------------------------------------------------
ELSE
   IF(lambda==lammin) REWIND iunit
   WRITE(iunit) lambda, nelem, numnp
WRITE(Out_unit,*) 'naph = ',naph
       WRITE(Out_unit,*)'mxelem,nelem=',mxelem,nelem
   WRITE(iunit) ((id(k,iel),k=1,9),iel=1,nelem)
   WRITE(iunit) (th(node),node=1,numnp)
   WRITE(iunit) thetaval
   WRITE(iunit) (ch(node),node=1,numnp)
   WRITE(iunit) chivals
   WRITE(iunit) (phi(index),index=1,numnp*naph)
   WRITE(iunit) energy,rhosurf
 ENDIF
 RETURN
 END
