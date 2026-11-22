SUBROUTINE symnod(ntheta, nchi, idt, nid, nodj, beta, nidm, ndisce, parity, jeven, symmetry, jmax, Virtualx,  mega)
USE Numeric_Kinds_Module
USE FileUnits_Module
USE FEM_Module
USE Numbers_Module
!
!     $RCSfile: symnod.f,v $ $Revision: 1.14 $
!     $Date: 89/10/18 14:18:16 $
!     $State: Stable $
!



!     I N P U T    A R G U M E N T S
!     ntheta
!     nchi
!     idt      Identifies each node as free, fixed, constrained, unique, nonexistant
!     nid(NDISCE) Array with the number of terms in the constraint equations
!     nodj
!     beta(ndisce) Constraint Coefficients
!     nidm
!     ndisce    number of displacement constraints.
!     parity
!     jeven
!     symmetry
!     jmax
!     Virtualx
!     mega

IMPLICIT NONE
LOGICAL jeven, symmetry
INTEGER idt, ntheta, nchi, jchi, itheta, nchipi, nchipi2, ndisce, nid, nodj
INTEGER nidm, jmax, jchip, parity,  mega, jinit, irefl, Virtualx
REAL(Kind=WP_Kind) beta, fac
DIMENSION idt(ntheta,nchi), nodj(nidm+1,*), beta(nidm,*), nid(*), Virtualx(ntheta,nchi)
!-----------------------------------------------------------------------
!     On entering:
!     ntheta                             The number of theta angles for
!     finest mesh.
!     nchi                             The number of chi angles for the
!     finest mesh.
!     idt(itheta,jchi)=  free             free node
!     idt(itheta,jchi)= nonxst      non existant node (in the fine mesh
!     but not used.
!     On exit:
!     idt(itheta,jchi)=  free             free node
!     idt(itheta,jchi)= nonxst      non existant node (in the fine mesh
!     but not used.
!     idt(itheta,jchi)= cnstrain    Constrained node.
!     nodj(ndisce)                    Constraining node number.
!     beta(1,ndisce)                     Constraining coefficient for t
!     ndisce-th constraint.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     I.  No reflection symmetry
!     A.  Even Parity
!     nodes at chi+pi are constrained to nodes at chi
!     Beta=1
!     Case 1    parity=0 no symmetry (internal z-axis is out of the page)

!     chi=pi/2
!     +
!     chi=pi  ________  chi=0
!
!     +
!     chi=3pi/2

!-----------------------------------------------------------------------
IF(.NOT.symmetry.AND.parity==0)THEN
   nchipi=(nchi+1)/2
   ndisce=0
   DO 10 jchi=nchipi+1,nchi-1
      DO 1 itheta=1,ntheta
         IF(idt(itheta,jchi)==nonxst)GO TO 1
         idt(itheta,jchi)=cnstrain
         ndisce=ndisce+1
         nid(ndisce)=1
         jchip=jchi-nchipi+1
         nodj(1,ndisce)=Virtualx(itheta,jchip)
         IF(itheta==1)nodj(1,ndisce)=1
         beta(1,ndisce)=one
         nodj(nidm+1,ndisce)=Virtualx(itheta,jchi)
         IF(nodj(1,ndisce)<=0.or.nodj(nidm+1,ndisce)<=0)THEN
            WRITE(Out_unit,*)'error in routine symnod'
            WRITE(Out_unit,*)ndisce, nid(ndisce), nodj(1,ndisce), nodj(nidm+1,ndisce), beta(1,ndisce)
            STOP 'symnod'
         ENDIF
 1          CONTINUE
 10      CONTINUE
   DO 11 itheta=1,ntheta
      IF(idt(itheta,1)==nonxst.AND.idt(itheta,nchipi)== nonxst.AND.idt(itheta,nchi)==nonxst)GOTO 11
      IF(idt(itheta,1)/=nonxst.AND.idt(itheta,nchipi)/=  nonxst)THEN
         idt(itheta,nchipi)=cnstrain
         ndisce=ndisce+1
         nid(ndisce)=1
         beta(1,ndisce)=one
         nodj(1,ndisce)=Virtualx(itheta,1)
         nodj(nidm+1,ndisce)=Virtualx(itheta,nchipi)
         IF(idt(itheta,nchi)/=nonxst)THEN
            idt(itheta,nchi)=cnstrain
            ndisce=ndisce+1
            nid(ndisce)=1
            beta(1,ndisce)=one
            nodj(1,ndisce)=Virtualx(itheta,1)
            nodj(nidm+1,ndisce)=Virtualx(itheta,nchi)
         ENDIF
         GOTO 11
      ELSEIF(idt(itheta,1)/=nonxst.AND.idt(itheta,nchi)/= nonxst)THEN
         idt(itheta,nchi)=cnstrain
         ndisce=ndisce+1
         nid(ndisce)=1
         beta(1,ndisce)=one
         nodj(1,ndisce)=Virtualx(itheta,1)
         nodj(nidm+1,ndisce)=Virtualx(itheta,nchi)
         GOTO 11
      ELSEIF(idt(itheta,nchipi)/=nonxst.AND. idt(itheta,nchi)/=nonxst)THEN
         idt(itheta,nchi)=cnstrain
         ndisce=ndisce+1
         nid(ndisce)=1
         beta(1,ndisce)=one
         nodj(1,ndisce)=Virtualx(itheta,nchipi)
         nodj(nidm+1,ndisce)=Virtualx(itheta,nchi)
      ELSE
         WRITE(Out_unit,*)'The program should never reach', ' this statement'
         STOP 'symnod'
      ENDIF
 11      CONTINUE
   GOTO 99
ENDIF
!-----------------------------------------------------------------------
!     B.  Odd Parity
!     nodes at chi+pi are constrained to nodes at chi
!     Beta=-1
!     Case 2    parity=1 no symmetry (internal z-axis is out of the page)

!     chi=pi/2
!     +
!     chi=pi  ________  chi=0
!
!     -
!     chi=3pi/2

!-----------------------------------------------------------------------
IF(.NOT.symmetry.AND.parity==1)THEN
   nchipi=(nchi+1)/2
   ndisce=0
   DO 20 jchi=nchipi+1,nchi
      DO 2 itheta=1,ntheta
         IF(idt(itheta,jchi)==nonxst)GO TO 2
         idt(itheta,jchi)=cnstrain
         ndisce=ndisce+1
         nid(ndisce)=1
         jchip=jchi-nchipi+1
         nodj(1,ndisce)=Virtualx(itheta,jchip)
         IF(itheta==1)nodj(1,ndisce)=1
         beta(1,ndisce)=-one
         nodj(nidm+1,ndisce)=Virtualx(itheta,jchi)
 2          CONTINUE
 20      CONTINUE
   GO TO 99
ENDIF
!     -------------------------------------------------------------------
!     II.  Reflection symmetry
!     jinit is 0 or 1 as initial j is even or odd.
!     beta is result of reflection about chi=pi/2
!     irefl is result of reflection about chi=0.
!     --------------------------------------------------------------
jinit=1
IF(jeven)jinit=0
irefl=(-1)**(jinit-mega)
!     -----------------------------------------------------------------
!     A.  Even Parity
!     i.  Even reflection symmetry in both reflections (A1)
!     Case 1'   parity=0 irefl=1  (internal z-axis is out of the page)

!     chi=pi/2
!     +  |  +
!     chi=pi  ____|____  chi=0
!     |
!     +  |  +
!     chi=3pi/2

!     -------------------------------------------------------------------
IF(symmetry.AND.parity==0.AND.irefl==1)THEN
   nchipi2=(nchi+1)/2
   ndisce=0
   DO 30 jchi=nchipi2+1,nchi
      DO 3 itheta=1,ntheta
         IF(idt(itheta,jchi)==nonxst)GO TO 3
         idt(itheta,jchi)=cnstrain
         ndisce=ndisce+1
         nid(ndisce)=1
         jchip=2*nchipi2-jchi
         nodj(1,ndisce)=Virtualx(itheta,jchip)
         IF(itheta==1)nodj(1,ndisce)=1
         beta(1,ndisce)=one
         nodj(nidm+1,ndisce)=Virtualx(itheta,jchi)
 3          CONTINUE
 30      CONTINUE
   GO TO 99
ENDIF
!-----------------------------------------------------------------------
!     ii.  Odd reflection symmetry in both reflections (A2)
!     Case 2'   parity=0 irefl=-1 (internal z-axis is out of the page)

!     chi=pi/2
!     -  |  +
!     chi=pi  ____|____  chi=0
!     |
!     +  |  -
!     chi=3pi/2

!-----------------------------------------------------------------------
IF(symmetry.AND.parity==0.AND.irefl==-1)THEN
   nchipi2=(nchi+1)/2
   ndisce=0
   DO 40 jchi=nchipi2+1,nchi-1
      DO 4 itheta=2,ntheta
         IF(idt(itheta,jchi)==nonxst)GO TO 4
         idt(itheta,jchi)=cnstrain
         ndisce=ndisce+1
         nid(ndisce)=1
         jchip=2*nchipi2-jchi
         nodj(1,ndisce)=Virtualx(itheta,jchip)
         IF(itheta==1)nodj(1,ndisce)=1
         beta(1,ndisce)=-one
         nodj(nidm+1,ndisce)=Virtualx(itheta,jchi)
 4          CONTINUE
 40      CONTINUE
   GO TO 99
ENDIF
!-----------------------------------------------------------------------
!     B.  Odd Parity
!     i.  Even reflection symmetry about pi/2, odd about 0 (B2)
!     Case 3'   parity=1 irefl=-1 (internal z-axis is out of the page)

!     chi=pi/2
!     +  |  +
!     chi=pi  ____|____  chi=0
!     |
!     -  |  -
!     chi=3pi/2

!-----------------------------------------------------------------------
IF(symmetry.AND.parity==1.AND.irefl==-1)THEN
   nchipi2=(nchi+1)/2
   ndisce=0
   DO 50 jchi=nchipi2+1,nchi-1
      DO 5 itheta=2,ntheta
         IF(idt(itheta,jchi)==nonxst)GO TO 5
         idt(itheta,jchi)=cnstrain
         ndisce=ndisce+1
         nid(ndisce)=1
         jchip=2*nchipi2-jchi
         nodj(1,ndisce)=Virtualx(itheta,jchip)
         IF(itheta==1)nodj(1,ndisce)=1
         beta(1,ndisce)=one
         nodj(nidm+1,ndisce)=Virtualx(itheta,jchi)
 5          CONTINUE
 50      CONTINUE
   GO TO 99
ENDIF
!-----------------------------------------------------------------------
!     ii. Odd reflection symmetry about pi/2, even about 0 (B1)
!     Case 4'   parity=1 irefl=1 (internal z-axis is out of the page)

!     chi=pi/2
!     -  |  +
!     chi=pi  ____|____  chi=0
!     |
!     -  |  +
!     chi=3pi/2

!-----------------------------------------------------------------------
IF(symmetry.AND.parity==1.AND.irefl==1)THEN
   nchipi2=(nchi+1)/2
   ndisce=0
   DO 60 jchi=nchipi2+1,nchi-1
      DO 6 itheta=2,ntheta
         IF(idt(itheta,jchi)==nonxst)GO TO 6
         idt(itheta,jchi)=cnstrain
         ndisce=ndisce+1
         nid(ndisce)=1
         jchip=2*nchipi2-jchi
         nodj(1,ndisce)=Virtualx(itheta,jchip)
         IF(itheta==1)nodj(1,ndisce)=1
         beta(1,ndisce)=-one
         nodj(nidm+1,ndisce)=Virtualx(itheta,jchi)
 6          CONTINUE
 60      CONTINUE
   GO TO 99
ENDIF
!-----------------------------------------------------------------------
!     Incorrect parameters
!-----------------------------------------------------------------------
WRITE(Out_unit,*)'***error execution stopped in routine symnod'
WRITE(Out_unit,*)'symmetry=',symmetry,'  parity=',parity,'  jeven=',jeven
STOP 'symnod'
 99   CONTINUE
IF(.NOT.symmetry.AND.parity==0.or. symmetry.AND.parity==0.AND.irefl==1)THEN
   fac=1
   DO 98 jchi=2,nchi
      IF(idt(1,jchi)==free)THEN
         idt(1,jchi)=cnstrain
         ndisce=ndisce+1
         nid(ndisce)=1
         nodj(1,ndisce)=Virtualx(1,1)
         beta(1,ndisce)=fac
         nodj(nidm+1,ndisce)=Virtualx(1,jchi)
      ENDIF
 98      CONTINUE
ENDIF
RETURN
END
