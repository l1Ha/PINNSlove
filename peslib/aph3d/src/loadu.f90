SUBROUTINE loadu(nnc, jmax, lcor, usys, chii, iskp, jskp, idiv,  dthe, dchi, nelem, rho, sinlist, philist,  &
           phirms, mxnode, mxtheta, mxchi, mxelmnt, ntheta, nchi)
USE Numeric_Kinds_Module
USE FileUnits_Module
USE Numbers_Module
!
!     $RCSfile: loadu.f,v $ $Revision: 1.14 $
!     $Date: 89/10/18 14:18:05 $
!     $State: Stable $
!

!     I N P U T    A R G U M E N T S
!
!     unc
!     jmax
!     lcor
!     usys
!     chii
!     iskp
!     jskp
!     idiv
!     dthe
!     dchi
!     nelem
!     rwho
!     sinlist
!     philist
!     phirms
!     mxnode
!     mxtheta
!     mxchi
!     mxelmnt
!     ntheta
!     nchi


IMPLICIT NONE
LOGICAL little, medium, full
INTEGER ithcall, ithsub
INTEGER jmax, iskp, nelem, mxnode,  mxtheta, mxchi, mxelmnt, ntheta, nchi, iel
INTEGER :: n, ibox, jbox, icor, mj, jcor, j, mi, i
INTEGER :: nnc(3), jskp(3), lcor(2,nelem), idiv(nelem)
REAL(Kind=WP_Kind) :: usys, dang, dthe, rho, phim, chi, theta, thetm
REAL(Kind=WP_Kind) :: dchi(3), chii(3), sinlist(mxelmnt), philist(mxelmnt), phirms(ntheta,nchi)

!     P U R P O S E   O F   S U B R O U T I N E

!     this routine returns the array philist which contains the largest valu
!     of the phirms for each element and array sinlist containing the
!     sin(theta) at the point of large phirms.


DATA ithcall/0/, ithsub/0/
DATA little/.false./, medium/.false./, full/.false./
CALL popt ('loadu   ', little, medium, full, ithcall, ithsub)

DO 165 iel=1,nelem

!     now go through each element and find the largest phirms in each elemen
!     store these phirms in an array indexed by element number. the followin
!     two loops look at all the nodes of the finest possible mesh
!     allowed by ndiv and keep the largest phirms found.

   CALL corner (iel, lcor, nnc, jmax, iskp, jskp, n, idiv, ibox,  jbox, dthe, dchi, dang, mxelmnt, icor, jcor)
   phim=zero
   DO 166 mj=0,2*jbox
      j=jcor+mj
      IF(n==1)THEN
         chi=(j-1)*dchi(n)
      ELSEIF (n==2)THEN
         chi=chii(2)+(j-nnc(1))*dchi(n)
      ELSE
         chi=chii(3)+(j-nnc(2)-nnc(1)+1)*dchi(n)
      ENDIF
      DO 167 mi=0,2*ibox
         i=icor+mi
         theta=(i-1)*dthe
         IF(phirms(i,j)>phim)THEN
            phim=phirms(i,j)
            thetm=theta
         ENDIF
 167        CONTINUE
 166     CONTINUE
   philist(iel) = phim
   IF(philist(iel)<=zero)philist(iel)=1.d-60
   sinlist(iel)=sin(thetm)
   IF(sinlist(iel)<=zero)sinlist(iel)=1.d-20
 165  CONTINUE
IF(full)THEN
   WRITE(Out_unit,*)'(philist(iel),iel=1,nelem) and nelem = ',nelem
!     WRITE(Out_unit,*)(philist(iel),iel=1,nelem)
   WRITE(Out_unit,*)'exiting routine loadu'
ENDIF
RETURN
END
