SUBROUTINE znodes(ntheta,nchi,idt,parity,jeven,symmetry, mega)
USE FileUnits_Module
USE FEM_Module
!
!     $RCSfile: znodes.f,v $ $Revision: 1.14 $
!     $Date: 89/10/18 14:18:18 $
!     $State: Stable $
!

!     I N P U T    A R G U M E N T S
!     ntheta
!     nchi
!     idt      Identifies each node as free, fixed, constrained, unique, nonexistant
!     parity
!     jeven
!     symmetry
!     mega

IMPLICIT NONE
LOGICAL jeven, symmetry
INTEGER idt, ntheta, nchi, jchi, itheta, npi2
INTEGER parity, nfixed, mega, jinit, irefl
!-----------------------------------------------------------------------
!     On entering:
!     ntheta                             The number of theta angles for
!     finest mesh.
!     nchi                             The number of chi angles for the
!     finest mesh.
!     idt(itheta,jchi)= free             free node
!     idt(itheta,jchi)= cnstrain    Constrained by other nodes.
!     idt(itheta,jchi)= nonxst      non existant node (in the fine mesh
!     but not used.
!     On exit:
!     idt(itheta,jchi)= free             free node
!     idt(itheta,jchi)= cnstrain    Constrained by other nodes.
!     idt(itheta,jchi)= nonxst      non existant node (in the fine mesh
!     but not used.
!     idt(itheta,jchi)= fixed             node removed or fixed at zero
!     [Only the nodes that were original
!     free nodes (idt=free) or constrain
!     (idt=cnstrain) nodes can be fixed
!     zero.]
!-----------------------------------------------------------------------
DIMENSION idt(ntheta,nchi)
nfixed=0
!-----------------------------------------------------------------------
!     I.  No reflection symmetry
!     A.  Even Parity [C*sin(2*n*chi)+D*cos(2*n*chi)]
!     No nodes are required to be zero.
!     Case 1    parity=0 no symmetry (internal z-axis is out of the page)

!     chi=pi/2
!     +
!     chi=pi  ________  chi=0
!
!     +
!     chi=3pi/2

!-----------------------------------------------------------------------
IF(.NOT.symmetry.AND.parity==0)GOTO 55
!-----------------------------------------------------------------------
!     B.  Odd Parity  [C*sin((2*n+1)*chi)+D*cos((2*n+1)*chi)]
!     Nodes at theta=0 must be zero.
!     Case 2    parity=1 no symmetry (internal z-axis is out of the page)

!     chi=pi/2
!     +
!     chi=pi  ________  chi=0
!
!     -
!     chi=3pi/2

!-----------------------------------------------------------------------
IF(.NOT.symmetry.AND.parity==1)THEN
   DO 1 jchi=1,nchi
      IF(idt(1,jchi)==nonxst)GO TO 1
      idt(1,jchi)=fixed
      nfixed=nfixed+1
 1       CONTINUE
   GOTO 55
ENDIF
!     -------------------------------------------------------------------
!     II.  Reflection symmetry
!     irefl is result of reflection about chi=0.
!     beta is result of reflection about chi=pi/2.
!     ---------------------------------------------------------------
jinit=1
IF(jeven)jinit=0
irefl=(-1)**(jinit-mega)
!     --------------------------------------------------------------
!     A.  Even Parity
!     i.  Even refl symm about 0 and pi/2  [cos(2*n*chi)] (A1)
!     No nodes are required to be zero.
!     Case 1'   parity=0 irefl=1  (internal z-axis is out of the page)

!     chi=pi/2
!     +  |  +
!     chi=pi  ____|____  chi=0
!     |
!     +  |  +
!     chi=3pi/2

!     -------------------------------------------------------------------
IF(symmetry.AND.parity==0.AND.irefl==1)GOTO 55
!     -------------------------------------------------------------------
!     ii.  Odd refl symm about 0 and pi/2  [sin(2*n*chi)] (A2)
!     Nodes at theta=0 must be zero.
!     Case 2'   parity=0 irefl=-1 (internal z-axis is out of the page)

!     chi=pi/2
!     -  |  +
!     chi=pi  ____|____  chi=0
!     |
!     +  |  -
!     chi=3pi/2

!     -------------------------------------------------------------------
IF(symmetry.AND.parity==0.AND.irefl==-1)THEN
   DO 2 jchi=1,nchi
      IF(idt(1,jchi)==nonxst)GO TO 2
      idt(1,jchi)=fixed
      nfixed=nfixed+1
 2       CONTINUE
!-----------------------------------------------------------------------
!     Nodes at chi=0 must be zero.
!-----------------------------------------------------------------------
   DO 23 itheta=1,ntheta
      IF(idt(itheta,1)==nonxst)GO TO 23
      idt(itheta,1)=fixed
      nfixed=nfixed+1
 23      CONTINUE
!-----------------------------------------------------------------------
!     Nodes at chi=pi must be zero.
!-----------------------------------------------------------------------
   DO 233 itheta=1,ntheta
      IF(idt(itheta,nchi)==nonxst)GO TO 233
      idt(itheta,nchi)=fixed
      nfixed=nfixed+1
 233     CONTINUE
!-----------------------------------------------------------------------
!     Nodes at chi=pi/2 must be zero.
!-----------------------------------------------------------------------
   npi2=(nchi+1)/2
   DO 24 itheta=1,ntheta
      IF(idt(itheta,npi2)==nonxst)GO TO 24
      idt(itheta,npi2)=fixed
      nfixed=nfixed+1
 24      CONTINUE
   GOTO 55
ENDIF
!-----------------------------------------------------------------------
!     B.  Odd Parity
!     i.  Even about pi/2, odd about 0 [sin((2*n+1)*chi)] (B2)
!     Nodes at theta=0 must be zero.
!     Case 3'   parity=1 irefl=-1 (internal z-axis is out of the page)

!     chi=pi/2
!     +  |  +
!     chi=pi  ____|____  chi=0
!     |
!     -  |  -
!     chi=3pi/2

!     -------------------------------------------------------------------
IF(symmetry.AND.parity==1.AND.irefl==-1)THEN
   DO 3 jchi=1,nchi
      IF(idt(1,jchi)==nonxst)GO TO 3
      idt(1,jchi)=fixed
      nfixed=nfixed+1
 3       CONTINUE
!     -------------------------------------------------------------------
!     Nodes at chi=0 must be zero.
!     -------------------------------------------------------------------
   DO 33 itheta=2,ntheta
      IF(idt(itheta,1)==nonxst)GO TO 33
      idt(itheta,1)=fixed
      nfixed=nfixed+1
 33      CONTINUE
!     -------------------------------------------------------------------
!     Nodes at chi=pi must be zero.
!-----------------------------------------------------------------------
   DO 333 itheta=2,ntheta
      IF(idt(itheta,nchi)==nonxst)GO TO 333
      idt(itheta,nchi)=fixed
      nfixed=nfixed+1
 333     CONTINUE
   GOTO 55
ENDIF
!-----------------------------------------------------------------------
!     ii.  Odd about pi/2, even about 0 [cos((2*n+1)*chi)] (B1)
!     Nodes at theta=0 must be zero.
!     Case 4'   parity=1 irefl=1 (internal z-axis is out of the page)

!     chi=pi/2
!     -  |  +
!     chi=pi  ____|____  chi=0
!     |
!     -  |  +
!     chi=3pi/2

!     -------------------------------------------------------------------
IF(symmetry.AND.parity==1.AND.irefl==1)THEN
   DO 4 jchi=1,nchi
      IF(idt(1,jchi)==nonxst)GO TO 4
      idt(1,jchi)=fixed
      nfixed=nfixed+1
 4       CONTINUE
!     -------------------------------------------------------------------
!     Nodes at chi=pi/2 must be zero.
!     -------------------------------------------------------------------
   npi2=(nchi+1)/2
   DO 44 itheta=2,ntheta
      IF(idt(itheta,npi2)==nonxst)GO TO 44
      idt(itheta,npi2)=fixed
      nfixed=nfixed+1
 44      CONTINUE
   GOTO 55
ENDIF
!-----------------------------------------------------------------------
!     Incorrect parameters
!-----------------------------------------------------------------------
WRITE(Out_unit,*)'***error execution stopped in routine znodes'
WRITE(Out_unit,*)'symmetry=',symmetry,'  parity=',parity,'  jeven=',jeven
STOP 'znodes'
 55   CONTINUE
RETURN
END
