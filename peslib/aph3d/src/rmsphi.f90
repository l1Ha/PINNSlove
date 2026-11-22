SUBROUTINE rmsphi (ntheta, nchi, phirms, nel, id, th, ch, nphi, numnp, phi, thetaval, chivals)
Use FileUnits_Module
USE Numeric_Kinds_Module

!            P U R P O S E   O F    S U B R O U T I N E

!-----------------------------------------------------------------------
! This routine calculates the probability density.
!
! On entering:
!     ntheta     number of theta angles.
!     nchi       number of chi angles.
!     nel        number of elements.
!     nphi       number of phi functions.
!     numnp      number of nodal points.
!     id         list of nodal numbers.
!     th
!     ch
!     phi        value of the surface functions at each nodal point.
!
! On exit:
!     phirms     probability density on the Virtualx grid.
!                The probability is the square root of the
!                mean of the sum of the surface functions squared.
!-----------------------------------------------------------------------

!            I N P U T    A R G U M E N T S
!   ntheta
!   nchi
!   phirms
!   nel
!   id
!   th
!   ch
!   nphi
!   numnp
!   phi
!   thetaval
!   chivals
LOGICAL little, medium, full
INTEGER ntheta, nchi, nel, numnp, k, iel, ithphi, nphi, jchi, itheta
INTEGER node, i1, i9, j1, j9, ithcall, ithsub, theta, chi
INTEGER id(9,nel), imap(9), th(numnp), ch(numnp)
INTEGER node1, node2, node3, node4, node5, node6, node7, node8, node9
REAL(Kind=WP_Kind) phival, phisum, r, s, gtot, gtheta
REAL(Kind=WP_Kind) gchi, r2, r2p, r2m, s2, s2p, s2m
REAL(Kind=WP_Kind) phirms(ntheta,nchi), phi(numnp,nphi), h(9)
REAL(Kind=WP_Kind) thetaval(ntheta), chivals(nchi)
INTRINSIC sqrt
EXTERNAL funct2, popt, mxoutd, grms
DATA imap / 2, 6, 3, 5, 9, 7, 1, 8, 4 /
!------------------------------------------------------------------
! Determine printing options.
!------------------------------------------------------------------
DATA little/.false./,medium/.false./,full/.false./,   ithcall/0/,ithsub/0/
CALL popt('rmsphi  ',little,medium,full,ithcall,ithsub)
IF(little)THEN
   WRITE(Out_unit,*)'rmsphi: ntheta, nchi, nel, nphi, numnp ', ntheta, nchi, nel, nphi, numnp
ENDIF
! temporary reads
REWIND TempNode_Bin_Unit
DO iel=1,nel
   DO k=1,9
      READ(TempNode_Bin_Unit)id(k,iel), theta, chi
      node=id(k,iel)
      th(node)=theta
      ch(node)=chi
   ENDDO
ENDDO
IF(full)THEN
   WRITE(Out_unit,*)'id: ',id
   WRITE(Out_unit,*)'th: ',th
   WRITE(Out_unit,*)'ch: ',ch
ENDIF
!-----------------------------------------------------------------------
! Initialize phirms to a negative number. Since phirms is nonegative
! this will allow one to make sure that the entire array has been
! defined.
!-----------------------------------------------------------------------
DO jchi=1,nchi
   DO itheta=1,ntheta
      phirms(itheta,jchi)=-1.0d0
   ENDDO
ENDDO
!-----------------------------------------------------------------------
! Evaluate phirms at the nodal points.
!-----------------------------------------------------------------------
DO iel=1,nel
   DO k=1,9
      node=id(k,iel)
      itheta=th(node)
      jchi=ch(node)
      IF(phirms(itheta,jchi)<0.0d0)THEN
         CALL grms(thetaval(itheta), chivals(jchi), gtot, gtheta, gchi)
         phisum=0.0d0
         DO ithphi=1,nphi
            phisum=phisum+(gtot*phi(node,ithphi))**2
         ENDDO
         phirms(itheta,jchi)=sqrt(phisum/nphi)
      ENDIF
   ENDDO
ENDDO
!-----------------------------------------------------------------------
! Interpolate to find phirms at all other points in the fine mesh.
!-----------------------------------------------------------------------
DO iel=1,nel
   node1=id(1,iel)
   node2=id(2,iel)
   node3=id(3,iel)
   node4=id(4,iel)
   node5=id(5,iel)
   node6=id(6,iel)
   node7=id(7,iel)
   node8=id(8,iel)
   node9=id(9,iel)
   i1=th(id(1,iel))
   i9=th(id(9,iel))
   j1=ch(id(1,iel))
   j9=ch(id(9,iel))
   IF(j9-j1/=1)THEN
   DO jchi=j1,j9
      DO itheta=i1,i9
!-----------------------------------------------------------------------
! IF phisum is nonegative the probability density has already been
! calculated hence skip the calculation.
!-----------------------------------------------------------------------
         IF(phirms(itheta,jchi)<0.d0)THEN
            CALL grms(thetaval(itheta), chivals(jchi), gtot, gtheta,  gchi)
            s=  REAL(i1+i9-2*itheta,WP_Kind)/REAL(i9-i1,WP_Kind)
            r= -REAL(j1+j9-2*jchi,WP_Kind)/REAL(j9-j1,WP_Kind)
            r2p=(r*r+r)/2d0
            r2m=(r*r-r)/2d0
            s2p=(s*s+s)/2d0
            s2m=(s*s-s)/2d0
            r2=1.0d0-r*r
            s2=1.0d0-s*s
!-----------------------------------------------------------------------
! Find the interpolation functions (h).
!-----------------------------------------------------------------------
            h(1)=r2p*s2p
            h(2)=r2m*s2p
            h(3)=r2m*s2m
            h(4)=r2p*s2m
            h(5)=r2*s2p
            h(6)=r2m*s2
            h(7)=r2*s2m
            h(8)=r2p*s2
            h(9)=r2*s2
            phisum=0.d0
            DO ithphi=1,nphi
               phival=h(2)*phi(node1,ithphi)
               phival=phival+h(6)*phi(node2,ithphi)
               phival=phival+h(3)*phi(node3,ithphi)
               phival=phival+h(5)*phi(node4,ithphi)
               phival=phival+h(9)*phi(node5,ithphi)
               phival=phival+h(7)*phi(node6,ithphi)
               phival=phival+h(1)*phi(node7,ithphi)
               phival=phival+h(8)*phi(node8,ithphi)
               phival=phival+h(4)*phi(node9,ithphi)
               phisum=phisum+(gtot*phival)**2
            ENDDO
            phirms(itheta,jchi)=sqrt(phisum/nphi)
         ENDIF
       ENDDO
    ENDDO
    ENDIF
ENDDO
!-----------------------------------------------------------------------
! Make sure that the phirms array has been completely defined.
!-----------------------------------------------------------------------
DO jchi=1,nchi
   DO itheta=1,ntheta
      IF(phirms(itheta,jchi)<0.d0)THEN
         WRITE(Out_unit,*) 'error in routine rmsphi: phirms is less than', ' zero'
         STOP 'rmsphi'
      ENDIF
   ENDDO
ENDDO
IF(full)THEN
   WRITE(Out_unit,*)' The phirms array: '
   CALL MxOutD(phirms, ntheta, nchi, 0)
ENDIF
WRITE(PhiRms1_Unit)phirms
RETURN
END
