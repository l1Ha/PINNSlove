SUBROUTINE evndiv(nelem, lcor, nnc, jmax, iskp, jskp, numnod, idiv, dthe, dchi, imax, ndiv, lnod, linkel, &
           numsto, lowel, iwork, idatwr, rho, usys,  nodmax, chii, sinlist, philist, phirms,              &
           mxnode, ntheta, nchi, mxelmnt, work1, work2, iwork33)
USE FileUnits_Module
USE opts_Module
USE tolexpc_Module
USE philis_Module
USE Numbers_Module
!
!     $RCSfile: evndiv.f,v $ $Revision: 1.14 $
!     $Date: 89/10/18 14:18:03 $
!     $State: Stable $
!

!     P U R P O S E   O F    S U B R O U T I N E

!     this SUBROUTINE takes the element link list and the ordered element
!     list made in crudely dividing the mesh and uses them to finely divide
!     the mesh.  in particular, this routine divides the elements
!     starting with those with the largest rms wavefunction.
!     it THEN divides the new mesh in the same fashion until
!     the number of mesh points is greater than nodmax, or all
!     criteria are satisfied.
!     numnod is a running count of how many nodes are in the mesh, nelnew
!     is how many new elements have been added and should be added to nelem
!     before leaving the routine to give the total number of elements.
!     linkel is a link list of elements which shows in what order the elemen
!     were made. it comes from SUBROUTINE divide. lowel is a ordered list of
!     element numbers, with the element with the largest phirms value in low
!     philist(iel) is the largest phirms in element iel.
!
!     I N P U T    A R G U M E N T S
!     nelem
!     lcor
!     nnc
!     jmax
!     iskp
!     jskp
!     numnod
!     idiv
!     dthe
!     dchi
!     imax
!     ndiv
!     lnod
!     linkel
!     numsto
!     lowel
!     iwork
!     idatwr
!     rho
!     usys
!     nodmax
!     chii
!     sinlist
!     philist
!     phirms
!     mxnode
!     ntheta
!     nchi
!     mxelmnt
!     work1
!     work2
!     iwork33


IMPLICIT NONE
LOGICAL didiv
INTEGER :: nelem, jmax, iskp, numnod, imax, ndiv, ibox, numsto, nodmax, mxnode
INTEGER :: ntheta, nchi, mxelmnt, icall, iel, idatwr, nelnew, iloop, ioel, jbox, jdiv, icor, jcor, n
INTEGER :: lcor(2,mxelmnt), nnc(3), idiv(mxelmnt), lnod(2*mxelmnt), linkel(mxelmnt), jskp(3)
INTEGER :: lowel(mxelmnt), iwork(3*mxelmnt), iwork33(2*mxelmnt)
REAL(Kind=WP_Kind) :: rho, tolsav, usys, tol, dang
REAL(Kind=WP_Kind) :: delth, aelchi, aelth, dthe, phimax, phim, tolfac, delchi
REAL(Kind=WP_Kind) :: chii(3), sinlist(mxelmnt), philist(mxelmnt), phirms(ntheta,nchi), work1(mxelmnt), work2(mxelmnt), dchi(3)
DATA icall /0/
!---------------------------------------------------------------------
!     calculate some tolerance parameters.
!---------------------------------------------------------------------
tolsav=(fortyeight*tolfn)**onethird/xkmesh
tol=tolsav/rho
WRITE(Out_unit,*)'tolfn=',tolfn
!
!     reset the link list of elements
!
DO iel=1,nelem
   linkel(iel) = iel+1
ENDDO

nelnew = 0

DO 100 iloop=1,mxelmnt

   didiv = .false.

   DO 90 ioel=1,nelem

      iel = lowel(ioel)
!     WRITE(Out_unit,*)'ioel=',ioel,' iel=',iel

!     IF the element has already been divided the maximum number of times,
!     DO not divide it.

      IF(idiv(iel)>=ndiv) GOTO 90

      CALL corner (iel, lcor, nnc, jmax, iskp, jskp, n, idiv, ibox, jbox, dthe, dchi, dang, mxelmnt, icor, jcor)
!     phimax is largest phirms on whole surface.
      IF(ioel==1) phimax=philist(iel)
!     phim is largest phirms on the current element.  save it because
!     hiphi redefines philist(iel) while adding new elements.  for
!     a given ioel and iel philist(iel) and sinlist(iel)
!     are kept ok until that element is divided.
      phim=philist(iel)
!     calc size of cell allowed by tolerance.
      IF(scheme==1)THEN
         tolfac=tol*(phimax/philist(iel))**onethird
      ELSEIF(scheme==2)THEN
         tolfac=tol*(phimax/philist(iel))**tolexp
      ENDIF
      delchi=tolfac/sinlist(iel)
      delth=2.*tolfac
!     compare with actual size of present cell and decide whether to divide
!     aelchi and aelth are the actual dimensions of the present cell.
      aelchi=2.*jbox*dchi(n)
      aelth=2.*ibox*dthe
      IF(delth>=aelth.AND.delchi>=aelchi) GOTO 90
!     divide the element into 4 new elements

      jdiv = idiv(iel) + 1
      CALL divide(mxelmnt, nelem, jbox, ibox, imax, jdiv, iel, icor, jcor, nelnew, numnod, idiv, lcor, lnod, linkel,  numsto)

      ibox = ibox/2
      jbox = jbox/2

      CALL hiphi(nnc, chii, dthe, dchi, jbox, ibox, iel, lcor, linkel, n, sinlist, philist, phirms,  &
                   mxnode, ntheta, nchi, mxelmnt, lowel, idiv)

      didiv = .true.

!     IF we have the desired number of nodal points, and the next largest
!     phirms is not equal to the current phirms, THEN we can quit.


      IF(2*numnod>=nodmax-75)THEN
         IF(ioel==nelem)THEN
            GOTO 110
         ELSEIF(phim>philist(lowel(ioel+1)))THEN
            GOTO 110
         ENDIF
!     this makes sure numnod never exceeds nodmax.
         IF(numnod>=nodmax-50) GOTO 110
      ENDIF
 90      CONTINUE


   IF(didiv)THEN
      nelem = nelem + nelnew
      nelnew = 0
      WRITE(Out_unit,*)'numnod=',numnod,' nelem=',nelem, ' at END of iloop=', iloop

!     WRITE(Out_unit,*)'philist before cormov 2 =',(philist(i), i=1,nelem)
!     WRITE(Out_unit,*)'linkel=',(linkel(i), i=1,nelem)
      CALL cormov(nelem, idiv, lcor, iwork, linkel, iwork33, sinlist, philist, work1, work2)

!     reset the link list of elements

      DO iel=1,nelem
         linkel(iel) = iel+1
      ENDDO

!     WRITE(Out_unit,*)'philist after cormov 2 =',(philist(i), i=1,nelem)
!     STOP 'evndiv'
      CALL lorder(philist, lowel, nelem, mxnode, mxelmnt)

   ELSE

      WRITE(Out_unit,*)' did not divide. All elements satisfied', ' tolfn or ndiv.'
      GOTO 110
   ENDIF

 100  CONTINUE
WRITE(Out_unit,*)'***error***: did not get enough nodes in evndiv!'
STOP 'evndiv'
110   CONTINUE
nelem=nelem+nelnew
CALL cormov(nelem, idiv, lcor, iwork, linkel, iwork33, sinlist, philist, work1, work2)

WRITE(Out_unit,*)'numnod=',numnod,' nelem=',nelem,' at END of evndiv'
RETURN
END
