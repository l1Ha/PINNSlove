SUBROUTINE eqncn2(ntheta, nchi, idt, nid, nodj, beta, nidm,  ndisce, virtualx, nelem, IntTheta, IntChi)
USE Numeric_Kinds_Module
USE FileUnits_Module
USE FEM_Module
USE Numbers_Module
!
!     $RCSfile: eqncn2.f,v $ $Revision: 1.14 $
!     $Date: 89/10/18 14:18:02 $
!     $State: Stable $
!

!     P U R P O S E   O F    S U B R O U T I N E


!-----------------------------------------------------------------------
!     This routine sets up the constraint equations for the nodal points
!     that are not in common with nodal points from adjacent
!     elements
!-----------------------------------------------------------------------



!     I N P U T    A R G U M E N T S
!
!     ntheta
!     nchi
!     idt      Identifies each node as free, fixed, constrained, unique, nonexistant
!     nid(NDISCE) Array with the number of terms in the constraint equations
!     nodj
!     beta(ndisce) Constraint Coefficients
!     nidm
!     ndisce    number of displacement constraints.
!     IntTheta
!     IntChi
!     virtualx
!     nelem


IMPLICIT NONE
INTEGER ntheta, nchi, idt, nid, nodj, nidm, ndisce, jchi, itheta, i1, i2, i3, j1, j2, j3, node, imap, &
           iel, nelem, ncom, virtualx, idt1, idt2, idt3, jdt1, jdt2, jdt3
INTEGER IntTheta(9,nelem), IntChi(9,nelem)
REAL(Kind=WP_Kind) beta, xval
DIMENSION idt(ntheta, nchi), nid(1), nodj(nidm+1,1), beta(nidm,1),  virtualx(ntheta,nchi), imap(9)
DATA imap/7,1,3,9,4,2,6,8,5/

DO jchi=1,nchi
   DO 21 itheta=1,ntheta
      IF(idt(itheta,jchi)/=UniqueNode)GOTO 21
      node=virtualx(itheta,jchi)
      DO 10 iel=1,nelem
         i1=IntTheta(imap(2),iel)
         i2=IntTheta(imap(6),iel)
         i3=IntTheta(imap(3),iel)
         idt1=idt(i1,jchi)
         idt2=idt(i2,jchi)
         idt3=idt(i3,jchi)
         j1=IntChi(imap(2),iel)
         j2=IntChi(imap(5),iel)
         j3=IntChi(imap(1),iel)
         jdt1=idt(itheta,j1)
         jdt2=idt(itheta,j2)
         jdt3=idt(itheta,j3)
         IF(i1==itheta.AND.j1==jchi)GOTO 10
         IF(i2==itheta.AND.j1==jchi)GOTO 10
         IF(i3==itheta.AND.j1==jchi)GOTO 10
         IF(i1==itheta.AND.j2==jchi)GOTO 10
         IF(i2==itheta.AND.j2==jchi)GOTO 10
         IF(i3==itheta.AND.j2==jchi)GOTO 10
         IF(i1==itheta.AND.j3==jchi)GOTO 10
         IF(i2==itheta.AND.j3==jchi)GOTO 10
         IF(i3==itheta.AND.j3==jchi)GOTO 10
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
         IF(i1==itheta.or.i3==itheta)THEN
            IF(jchi<j1.or.jchi>j3)GOTO 10
            IF(jdt1/=fixed.AND.jdt2/=fixed.AND. jdt3/=fixed)THEN
               ndisce=ndisce+1
               nid(ndisce)=3
               nodj(1,ndisce)=virtualx(itheta,j1)
               nodj(2,ndisce)=virtualx(itheta,j2)
               nodj(3,ndisce)=virtualx(itheta,j3)
               nodj(nidm+1,ndisce)=node
               xval=(jchi-j1)*(jchi-j2)/REAL((j3-j1)*(j3-j2),WP_Kind) -(jchi-j2)*(jchi-j3)/REAL((j1-j2)*(j1-j3),WP_Kind)
               beta(1,ndisce)=half*xval*(xval-one)
               beta(2,ndisce)=-(xval+one)*(xval-one)
               beta(3,ndisce)=half*xval*(xval+one)
            ELSEIF(jdt1==fixed.AND.jdt2/=fixed.AND. jdt3/=fixed)THEN
               ndisce=ndisce+1
               nid(ndisce)=2
               nodj(1,ndisce)=virtualx(itheta,j2)
               nodj(2,ndisce)=virtualx(itheta,j3)
               nodj(nidm+1,ndisce)=node
               xval=(jchi-j1)*(jchi-j2)/REAL((j3-j1)*(j3-j2),WP_Kind) -(jchi-j2)*(jchi-j3)/REAL((j1-j2)*(j1-j3),WP_Kind)
               beta(1,ndisce)=-(xval+one)*(xval-one)
               beta(2,ndisce)=half*xval*(xval+one)
            ELSEIF(jdt1/=fixed.AND.jdt2==fixed.AND.jdt3/=fixed)THEN
               ndisce=ndisce+1
               nid(ndisce)=2
               nodj(1,ndisce)=virtualx(itheta,j1)
               nodj(2,ndisce)=virtualx(itheta,j3)
               nodj(nidm+1,ndisce)=node
               xval=(jchi-j1)*(jchi-j2)/REAL((j3-j1)*(j3-j2),WP_Kind) -(jchi-j2)*(jchi-j3)/REAL((j1-j2)*(j1-j3),WP_Kind)
               beta(1,ndisce)=half*xval*(xval-one)
               beta(2,ndisce)=half*xval*(xval+one)
            ELSEIF(jdt1/=fixed.AND.jdt2/=fixed.AND. jdt3==fixed)THEN
               ndisce=ndisce+1
               nid(ndisce)=2
               nodj(1,ndisce)=virtualx(itheta,j1)
               nodj(2,ndisce)=virtualx(itheta,j2)
               nodj(nidm+1,ndisce)=node
               xval=(jchi-j1)*(jchi-j2)/REAL((j3-j1)*(j3-j2),WP_Kind) -(jchi-j2)*(jchi-j3)/REAL((j1-j2)*(j1-j3),WP_Kind)
               beta(1,ndisce)=half*xval*(xval-one)
               beta(2,ndisce)=-(xval+one)*(xval-one)
            ELSE
               STOP 'eqnc2'
            ENDIF
            GOTO 21
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
         ELSEIF(j1==jchi.or.j3==jchi)THEN
            IF(itheta<i1.or.itheta>i3)GOTO 10
            IF(idt1/=fixed.AND.idt2/=fixed.AND.idt3/=fixed)THEN
               ndisce=ndisce+1
               nid(ndisce)=3
               nodj(1,ndisce)=virtualx(i1,jchi)
               IF(i1==1)nodj(1,ndisce)=virtualx(1,1)
               nodj(2,ndisce)=virtualx(i2,jchi)
               nodj(3,ndisce)=virtualx(i3,jchi)
               nodj(nidm+1,ndisce)=node
               xval=(itheta-i1)*(itheta-i2)/REAL((i3-i1)*(i3-i2),WP_Kind)-(itheta-i2)*(itheta-i3)/ REAL((i1-i2)*(i1-i3),WP_Kind)
               beta(1,ndisce)=half*xval*(xval-one)
               beta(2,ndisce)=-(xval+one)*(xval-one)
               beta(3,ndisce)=half*xval*(xval+one)
            ELSEIF(idt1==fixed.AND.idt2/=fixed.AND.  idt3/=fixed)THEN
               ndisce=ndisce+1
               nid(ndisce)=2
               nodj(1,ndisce)=virtualx(i2,jchi)
               nodj(2,ndisce)=virtualx(i3,jchi)
               nodj(nidm+1,ndisce)=node
               xval=(itheta-i1)*(itheta-i2)/REAL((i3-i1)*(i3-i2),WP_Kind)-(itheta-i2)*(itheta-i3)/ REAL((i1-i2)*(i1-i3),WP_Kind)
               beta(1,ndisce)=-(xval+one)*(xval-one)
               beta(2,ndisce)=half*xval*(xval+one)
            ELSEIF(idt1/=fixed.AND.idt2==fixed.AND. idt3/=fixed)THEN
               ndisce=ndisce+1
               nid(ndisce)=2
               nodj(1,ndisce)=virtualx(i1,jchi)
               IF(i1==1)nodj(1,ndisce)=virtualx(1,1)
               nodj(2,ndisce)=virtualx(i3,jchi)
               nodj(nidm+1,ndisce)=node
               xval=(itheta-i1)*(itheta-i2)/REAL((i3-i1)*(i3-i2),WP_Kind)-(itheta-i2)*(itheta-i3)/REAL((i1-i2)*(i1-i3),WP_Kind)
               beta(1,ndisce)=half*xval*(xval-one)
               beta(2,ndisce)=half*xval*(xval+one)
            ELSEIF(idt1/=fixed.AND.idt2/=fixed.AND. idt3==fixed)THEN
               ndisce=ndisce+1
               nid(ndisce)=2
               nodj(1,ndisce)=virtualx(i1,jchi)
               IF(i1==1)nodj(1,ndisce)=virtualx(1,1)
               nodj(2,ndisce)=virtualx(i2,jchi)
               nodj(nidm+1,ndisce)=node
               xval=(itheta-i1)*(itheta-i2)/REAL((i3-i1)*(i3-i2),WP_Kind)-(itheta-i2)*(itheta-i3)/REAL((i1-i2)*(i1-i3),WP_Kind)
               beta(1,ndisce)=half*xval*(xval-one)
               beta(2,ndisce)=-(xval+one)*(xval-one)
            ELSE
               STOP 'eqncn2'
            ENDIF
            GOTO 21
         ENDIF
 10         CONTINUE
      WRITE(Out_unit,*)'error in routine eqncn2'
      WRITE(Out_unit,*)'itheta,jchi = ',itheta,jchi
      STOP 'eqncn2'
 21      CONTINUE
ENDDO
ncom=0
DO jchi=1,nchi
   DO itheta=1,ntheta
      IF(idt(itheta,jchi)==UniqueNode)THEN
         ncom=ncom+1
         idt(itheta,jchi)=cnstrain
      ENDIF
   ENDDO
ENDDO
RETURN
END
