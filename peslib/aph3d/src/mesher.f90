SUBROUTINE mesher (nlt, nlc, ndiv, id, th, ch, mxelmnt, numnp,  nelem, chimax, mapem, nodmax, idiv,            &
           numnod, lnod, IntTheta, IntChi, lcor, scr1,  scr2, iscr1, iscr2, iscr3, mxtheta,                    &
           mxchi, mapr, ido, idivo, lcoro, mapnod, sinlist, philist, phirms, mxnode, work1, work2, iwork33,    &
           ithrho, ntheta, nchi)
USE Numeric_Kinds_Module
USE FileUnits_Module
USE delang_Module
USE Numbers_Module
USE Masses_Module
USE rhosur_Module
USE totj_Module
USE rhowrt_Module
!
!     $RCSfile: mesher.f,v $ $Revision: 1.14 $
!     $Date: 89/10/18 14:18:08 $
!     $State: Stable $
!

!     P U R P O S E   O F    S U B R O U T I N E

!-----------------------------------------------------------------------
!     this routine determines the element mesh and nodal DATA.
!     this routine uses nlt, nlc, and ndiv to determine
!     a Virtualx grid that divides the space into the smallest
!     units that will be possible;  it THEN knows the smallest possible
!     increments of theta and chi; IF all the points were turned on
!     there would be rows and columns of points, labeled by i and j
!     respectively; imax is the number of rows in the theta direction and
!     jmax is the number of columns in the chi direction; .
!     mxtheta must be at least as large as imax detnd below.
!     mxchi must be at least as large as jmax before jmax is halved.
!-----------------------------------------------------------------------


!     I N P U T    A R G U M E N T S
!     nlt
!     nlc
!     ndiv
!     id(numnp)   id=-nce if the is constrained and positive otherwise
!     th
!     ch
!     mxelmnt
!     numnp    Number of nodal poin
!     nelem
!     chimax
!     mapem
!     nodmax
!     idiv
!     numnod
!     lnod
!     IntTheta
!     IntChi
!     lcor
!     scr1
!     scr2
!     iscr1
!     iscr2
!     iscr3
!     mxtheta
!     mxchi
!     mapr
!     ido
!     idivo
!     lcoro
!     mapnod
!     sinlist
!     philist
!     phirms
!     mxnode
!     work1
!     work2
!     iwork33
!     ithrho
!     ntheta
!     nchi


IMPLICIT NONE
LOGICAL mapem
INTEGER :: nlt, nlc(3), ndiv, id(9,mxelmnt), mxelmnt, numnp, nelem, nodmax, numnod, IntTheta(9,mxelmnt), IntChi(9,mxelmnt)
INTEGER :: mxchi, mxnode, ithrho, ntheta, nchi, numsto, idatwr
INTEGER :: imax, nltmax, jmax, i, nnc(3), iskp, j, nmax, nrc, iel, jcor, n, lc, lt, icor, jskp1, nod
INTEGER :: ijnode, nmndlp, inode, nodlap, nodend, nodnum, indx, ibit, jmask, in, node, iboth, nelem2, ibox, jbox
INTEGER :: ibcor, jbcor, m, l, k2, jmaxsv, maxelp, ielp, k, mxtheta
INTEGER :: kmap(9), jskp(3), lnod(2*mxelmnt), lcor(2,mxelmnt), nlcrmx(3), idiv(mxelmnt)
INTEGER :: iscr1(mxtheta*mxchi), iscr2(mxtheta*mxchi), iscr3(*), iwork33(mxtheta*mxchi)
INTEGER :: mapr(*), ido(*), idivo(*), lcoro(*), mapnod(*)
REAL(Kind=WP_Kind) :: th(9,mxelmnt), ch(9,mxelmnt), scr1(mxnode), scr2(mxnode), sinlist(mxelmnt)
REAL(Kind=WP_Kind) :: philist(mxelmnt), phirms(ntheta,nchi), work1(mxelmnt), work2(mxelmnt)
REAL(Kind=WP_Kind) :: chii(3), thetaval(mxtheta*mxchi), chivals(mxtheta*mxchi), dang, chi, theta, chimax

DATA kmap /7,8,9,4,5,6,1,2,3/
DATA numsto /30/
idatwr=0
IF(jtot<mega)THEN
   WRITE(Out_unit,*) '***error***: mega is greater than jtot'
   WRITE(Out_unit,*) 'jtot=',jtot,' mega=',mega
ENDIF
!-----------------------------------------------------------------------
!     find the singularity lines
!     note that chii(1) < chii(3) < chii(2) < chimax
!-----------------------------------------------------------------------
chii(1)=zero
chii(2)=pi+atan2(-mass(3),-usys)
chii(3)=atan2(mass(2),-usys)
IF(ithrho==1)WRITE(Out_unit,*) 'chii=',chii
!     after calculating chi angles, switch usys to au.

!     determine the maximum number of elements in each region and
!     THEN find the maximum number of independent nodes

nltmax=nlt*(2**ndiv)
imax=nltmax*2+1
jmax=0
DO i=1,3
   nlcrmx(i)=nlc(i)*(2**ndiv)
   nnc(i)=nlcrmx(i)*2+1
   jmax=jmax+nnc(i)
ENDDO
jmax=jmax-2
IF(symmetry)jmax=(jmax+1)/2

!     compute the smallest angle displacement available in each region

dchi(1)=chii(2)/(nnc(1)-1)
dchi(2)=(chii(3)-chii(2))/(nnc(2)-1)
dchi(3)=(pi-chii(3))/(nnc(3)-1)
dthe=halfpi/(imax-1)

IF(ithrho==1)THEN
   WRITE(Out_unit,*) 'nlt=',nlt,' ndiv=',ndiv,' nlc=',nlc
   WRITE(Out_unit,*) 'nltmax=',nltmax,' imax=',imax,' jmax=',jmax
   WRITE(Out_unit,*) 'nlcrmx=',nlcrmx,' nnc=',nnc
   WRITE(Out_unit,*) 'dthe=',dthe
   WRITE(Out_unit,*) 'dchi=',dchi
ENDIF

!     this completes the Virtualx (fine) mesh.

!     find the number of nodes to skip in making the coarse grid

iskp=(imax-1)/nlt
DO j=1,3
   jskp(j)=(nnc(j)-1)/nlc(j)
ENDDO

!     lnod keeps track of whether a node is "turned on" or not.  IF there
!     is a 1 in the nth bit of lnode(k)THEN that node is
!     being included in the final mesh; each kth part of lnode will be used
!     to represent numsto nodes (numsto=60 on the cray, 30 on the vax); star
!     numbering the nodes at i=j=0 and increase i fastest (theta direction)
!     we go all the way from i=1 to i=imax for each j.

!     make the coarse grid

numnod=0
nmax=3
IF(symmetry)THEN
   nmax=2
   IF(nlc(1)/=nlc(3))THEN
      WRITE(Out_unit,*)'***error***:symmetry=true requires', ' nlc(1)=nlc(3)'
      STOP 'mesher'
   ENDIF
   IF(mod(nlc(2),2)/=0)THEN
      WRITE(Out_unit,*) '***error***:symmetry=true requires', ' nlc(2)=even'
      STOP 'mesher'
   ELSE
      nlc(2)=nlc(2)/2
   ENDIF
ENDIF
lnod=0
nrc=0
iel=0
jcor=1
DO n=1,nmax
   nrc=jcor
   DO lc=1,nlc(n)
      DO lt=1,nlt
         icor=1+(lt-1)*iskp
         iel=iel+1
         lcor(1,iel)=icor
         lcor(2,iel)=jcor
         iscr2(iel)=iel+1
         jskp1=jskp(n)
         CALL nodson (icor, jcor, iskp, jskp1, numnod,  lnod, imax, mxelmnt, numsto)
      ENDDO
      jcor=nrc+lc*jskp(n)
   ENDDO
ENDDO
nelem=iel
!      WRITE(Out_unit,*)'(lnod(i),i=1,50)'
!      WRITE(Out_unit,*) (lnod(i),i=1,50)
idiv=0
!     this completes the coarse grid.

!     rearrange the element's corner and divison arrays so that they corresp
!     with the new element numbers due to having new elements. note that thi
!     effectively renumbers the elements. in this CALL, scr2 is a dummy arra
!     that is just needed to fill out the CALL requirements.

CALL cormov(nelem, idiv, lcor, iscr3, iscr2, iscr3(mxelmnt), sinlist, philist, work1, work2)

!     divide the mesh only IF ndiv gt 0

IF(ndiv>0)THEN
!     READ the current phirms from disc PhiRms2_Unit

   IF(ithrho>1)THEN
      READ(PhiRms2_Unit) phirms
   ELSE
      CALL phint (phirms, ntheta, nchi, thetaval, chivals)
   ENDIF

!     now go through each element and find the largest phirms in each elemen
!     store these phirms in an array philist indexed by element number.

   CALL loadu(nnc, jmax, lcor, usys, chii, iskp, jskp, idiv, dthe, dchi, nelem, rhosurf, sinlist, philist, &
              phirms, mxnode, mxtheta, mxchi, mxelmnt, ntheta, nchi)

!     keep an array that is a list of elements going from largest
!     phirms to smallest.

   CALL lorder(philist, iscr3, nelem, mxnode, mxelmnt)

!     starting with the largest phirms, divide each element until all eleme
!     satisfy the tolerance or until the max no of nodal points is reached.

   CALL evndiv(nelem, lcor, nnc, jmax, iskp, jskp, numnod, idiv, dthe, dchi, imax, ndiv, lnod, iscr2, numsto,         &
              iscr3, iscr3(3*mxelmnt), idatwr, rhosurf, usys, nodmax, chii, sinlist, philist, phirms, mxnode, ntheta, &
              nchi, mxelmnt, work1, work2,iwork33)

ENDIF

!     number the nodes globally.
!     renumber the nodes so that they are in an order
!     usful to the surface function routine; that means going from i-j label
!     to elements of nine nodes and giving each nodal point an unique
!     nodal number (a global number)

nod=0
ijnode=0
nmndlp=0
DO i=1,imax
   inode=0
   nodlap=0
   nodend = 0
   DO j=1,jmax

!     we now find the node's position in the on/off list;
!     IF ion=1 the node has been turned on and we have another node
!     nodnum is the node's position on the grid, indx is which array element
!     in lnod the node is stored in, and ibit is which bit the node is.

      nodnum=imax*(j-1)+i
      indx=(nodnum-1)/numsto+1
      IF(indx>2*mxelmnt)THEN
         WRITE(Out_unit,*) '***error***: indx too big:',indx,2*mxelmnt
         STOP 'mesher'
      ENDIF
!            WRITE(Out_unit,*)'nodnum=',nodnum,' indx=',indx,' numsto=',numsto
      ibit = nodnum - (indx-1)*numsto-1
      IF(ibit>numsto)THEN
         WRITE(Out_unit,*) '***error***:ibit too big: ibit =',ibit
         WRITE(Out_unit,*) 'indx=',indx,' nodnum=',nodnum, ' i=',i,' j=',j
         STOP 'mesher'
      ENDIF
!      WRITE(Out_unit,*)'ibit=',ibit,' (2**ibit)=',2**ibit
      jmask = 2**ibit
!     use ior on vax and or on the cray
!     IF(ior(lnod(indx),jmask)==lnod(indx))THEN
      IF(ior(lnod(indx),jmask)==lnod(indx))THEN
!     IF this test is true, THEN there was a node turned on at this position
!     and the global node number, nod, is incremented and stored; note that
!     we are only looking at the lower half of the mesh currently
!      WRITE(Out_unit,*)'found a node in mesher'
         nod=nod+1
         iscr1(j+(i-1)*jmax)=nod

!     we will have to double the range so we
!     have to keep track of how many nodes there are on this row up to
!     the symmetry point;  we also need to know how many nodes lie on the
!     symmetry column. inode is how many nodes are on this constant i line,
!     nodlap is wether or not the node at the symmetry line is on, and iscr2
!     is the index between the global node number and the position on this r
!     nodend is how many nodes are on the j=1 line (used by
!     symmetry=false only)

         inode=inode+1
         iscr2(inode)=j
         IF(j==jmax)nodlap=1
         IF(j==1)nodend = 1
      ENDIF
   ENDDO

!     double the number of nodes in this row

!     nod is the global node number; inode is the number of nodes in the ith
!     row between chi=0 and chi=chisym; j is the global node number of the
!     in'th node in the lower half plane; node is the global node number of
!     the node in the upper half plane corresponding to the in'th node.
!     nodlap=1 IF the node at chi=chisym is turned on and nodend=1 IF the no
!     at chi=0 IF turned on (both are for the i'th row only)

   IF(symmetry)THEN
!     symmetry=true does a reflection about the symmetry axis
      DO in=1,inode
         j=iscr2(in)
         node=ijnode+2*inode-in
         IF(nodlap==0) node=node+1
         iscr3(2*jmax-j+(i-1)*jmax)=node
      ENDDO
      nod=nod+inode-nodlap
      IF(nodlap==1) nmndlp=nmndlp+1

   ELSEIF(.NOT.symmetry)THEN
!     symmetry=false has a identical node at chi2=chi+pi for every node
      iboth = 0
      IF(nodlap==1.AND.nodend==1)THEN
         nmndlp = nmndlp + 1
         iboth=1
      ENDIF
      DO in=1,inode
         j =  iscr2(in)
         node = ijnode + inode + in - iboth
         iscr3(2*jmax-j+(i-1)*jmax)=node
      ENDDO
      nod = nod + inode - iboth
   ENDIF

   ijnode=nod
ENDDO

IF(nod>mxnode)THEN
   WRITE(Msg_unit,*) '***error***: too many nodes !'
   WRITE(Msg_unit,*) 'numnp=',nod,' mxnode=',mxnode
   WRITE(Out_unit,*) '***error***: too many nodes !'
   WRITE(Out_unit,*) 'numnp=',nod,' mxnode=',mxnode
   STOP 'mesher'
ENDIF

!     store the number of nodes and elements in the doubled mesh

numnp=nod
nelem2=2*nelem

IF(nelem2>mxelmnt)THEN
   WRITE(Out_unit,*) '***error***:too many elements to double'
   STOP 'mesher'
ENDIF
IF(numnp/=(numnod*2-nmndlp))THEN
   WRITE(Out_unit,*) '***error***: numnp/=numnod'
   WRITE(Out_unit,*) 'numnp=',numnp,' numnod=',numnod*2-nmndlp
   WRITE(Out_unit,*) 'nmndlp=',nmndlp
   STOP 'mesher'
ENDIF
WRITE(Out_unit,*) 'numnp=',numnp,' full nelem=',nelem2

!     load the output arrays for the surface function code so that
!     the node information is presented by local node number in each
!     element;

DO iel=1,nelem
   CALL corner (iel, lcor, nnc, jmax, iskp, jskp, n, idiv, ibox,jbox, dthe, dchi, dang, mxelmnt, ibcor, jbcor)
   DO m=0,2
      j=jbcor+m*jbox
      IF(n==1)THEN
         chi=(j-1)*dchi(n)
      ELSEIF (n==2)THEN
         chi=chii(2)+(j-nnc(1))*dchi(n)
      ELSE
         chi=chii(3)+(j-nnc(2)-nnc(1)+1)*dchi(n)
      ENDIF
      DO l=0,2
         i=ibcor+l*ibox
         k=3*m+l+1
         theta=(i-1)*dthe
         id(k,iel)=iscr1(j+(i-1)*jmax)
         th(k,iel)=theta
         ch(k,iel)=chi

!     IF the range is doubled, load the upper half of the arrays

         IF(symmetry)THEN
            k2=kmap(k)
         ELSE
            k2 = k
         ENDIF
         id(k2,iel+nelem)=iscr3(2*jmax-j+(i-1)*jmax)
      ENDDO
   ENDDO
ENDDO
!  double the range of chi IF necessary; note that
!  double doubles nelem and maps the new elements into the old; it does
!  not number the upper half of the nodes globally or double numnp

jmaxsv=jmax
IF(symmetry)THEN
   CALL double (nelem, chimax, th, ch, id, maxelp, idiv, lcor, mapr, mxelmnt)
   jmax=2*jmax-1
   nlc(2)=nlc(2)*2
ELSEIF(.NOT.symmetry)THEN
   WRITE(Out_unit,*)'calling doubpi, nelem=',nelem
   CALL doubpi(nelem, th, ch, maxelp, mxelmnt, mapr)
ENDIF

!     WRITE the mesh out to elmdef (not used by sfunc)

WRITE(Out_unit,*)'nelem=',nelem,' just before elmdef WRITE'
IF(rhosurf>=rhowrit)THEN
   WRITE(ElmDef_Unit) nelem,symmetry,nlt,nlc,ndiv,rhosurf
   WRITE(ElmDef_Unit) ((id(k,iel),k=1,9),iel=1,nelem)
   WRITE(ElmDef_Unit) ((th(k,iel),k=1,9),iel=1,nelem)
   WRITE(ElmDef_Unit) ((ch(k,iel),k=1,9),iel=1,nelem)
ENDIF
DO iel=1,nelem
   ielp=iel
   IF(iel>nelem/2)ielp=iel-nelem/2
   IF(symmetry.AND.(iel==nelem.or.iel==maxelp))THEN
      IF(iel==nelem)THEN
         ielp=maxelp-nelem/2
      ELSE
         ielp=nelem/2
      ENDIF
   ENDIF
   CALL corner(ielp, lcor, nnc, jmax, iskp, jskp, n, idiv, ibox, jbox, dthe, dchi, dang, mxelmnt, ibcor, jbcor)
   DO m=0,2
      j=jbcor+m*jbox
      IF(symmetry.AND.iel>nelem/2)j=jbcor+2*jbox-m*jbox
      IF(iel>nelem/2)THEN
         IF(.NOT.symmetry)j=j+jmaxsv-1
         IF(symmetry)j=jmax-j+1
      ENDIF
      DO l=0,2
         i=ibcor+l*ibox
         k=3*m+l+1
         node=id(k,iel)
         scr1(node)=th(k,iel)
         scr2(node)=ch(k,iel)
         iscr3(i+(j-1)*imax)=node
         IntTheta(k,iel)=i
         IntChi(k,iel)=j
      ENDDO
   ENDDO
ENDDO

!     IF desired the mesh from the last rho value is mapped
!     into the mesh just created for this rho value;  the assumptions
!     that mapel makes require that the difference in rho must be small
!     so that the mesh changes very little

!     IF(mapem)THEN
   CALL mapel (id, lcor, mxelmnt, idiv, rhosurf, nelem, symmetry, ndiv, &
            numnp, maxelp, mega, idatwr, ido, idivo,  lcoro, mapnod)
!     ENDIF

RETURN
ENDSUBROUTINE Mesher
