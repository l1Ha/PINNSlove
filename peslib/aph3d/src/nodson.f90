SUBROUTINE nodson (icor, jcor, iskp, jskp1, n, lnod, imax, mxelmnt, numsto)
USE FileUnits_Module
!
!     $RCSfile: nodson.f,v $ $Revision: 1.14 $
!     $Date: 89/10/18 14:18:12 $
!     $State: Stable $
!


!     P U R P O S E   O F    S U B R O U T I N E
!     this SUBROUTINE takes an element and turns on its 9 node points
!     (icor,jcor) specifies the upper left hand corner of the
!     element and iskp, jskp specify the length of the sides  in
!     units of the maxium desity nodes


!     I N P U T    A R G U M E N T S
!     icor
!     jcor
!     iskp
!     jskp1
!     n
!     lnod
!     imax
!     mxelmnt
!     numsto

IMPLICIT NONE
INTEGER i, ibit, icor, imax, indx, iskp, j, jcor, jmask, jskp1,  k, l, lnodod, mxelmnt, n, nodnum, numsto
INTEGER lnod(2*mxelmnt)

!     lnod has size mxelmnt when numsto=30, mxelmnt when numsto=60

!     loop thru the 9 nodes in the new element finding the row and column nu
!     (i&j) (theta and chi)
DO k=0,2
   j=jcor+k*jskp1/2
   DO l=0,2
      i=icor+l*iskp/2
!     find the global node number
      nodnum=imax*(j-1)+i
!     find which array element this node is in
      indx=(nodnum-1)/numsto+1
      IF(indx>2*mxelmnt)THEN
         WRITE(Out_unit,*) '***error***:indx too big',indx,2*mxelmnt
         STOP 'nodson'
      ENDIF
!     compute which bit in the array element we want to turn on
      ibit = nodnum-(indx-1)*numsto-1
      IF(ibit>numsto)THEN
         WRITE(Out_unit,*) '***error***:ibit too big: ibit =',ibit
         WRITE(Out_unit,*) 'indx=',indx,' nodnum=',nodnum, ' i=',i,' j=',j
         STOP 'nodson'
      ENDIF
!     make the mask needed to turn on the desired bit
      jmask = 2**ibit
      WRITE(Out_unit,*)'nodnum=',nodnum,' indx=',indx,' ibit=',ibit,' jmask=',jmask
!     turn on the desired bit;
!     use ior on the vax; this is an inclusive or hence IF either or both argum
!     have the bit turned on, the resultant has the bit turned on.
!     use or on the cray (bit-by-bit LOGICAL sum)
      lnodod = lnod(indx)
!     lnod(indx) = ior(lnodod,jmask)
      lnod(indx) = ior(lnodod,jmask)
      WRITE(Out_unit,*)'lnodod=',lnodod,' lnod(indx)=',lnod(indx)
!     IF the array element changed, THEN we turned on a new node and must c
      IF(lnod(indx)/=lnodod)n=n+1
      WRITE(Out_unit,*)'n=',n,' in nodson'
   ENDDO
ENDDO
!     WRITE(Out_unit,*)'n=',n,' in nodson at END of loop'
!     WRITE(Out_unit,*)'lnod in nodson',(lnod(i), i=1,50)
RETURN
END
