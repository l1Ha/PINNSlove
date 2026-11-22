SUBROUTINE addres(maxa, mht, neq, nwk, ma)
USE FileUnits_Module
!
!   $RCSfile: addres.f,v $ $Revision: 1.15 $
!   $Date: 89/11/16 11:02:49 $
!   $State: Stable $

!   P U R P O S E  O F   S U B R O U T I N E

!   I N P U T  A R G U M E N T S
!  maxa
!     mht       column heights.
!  neq
!  nwk
!  ma

!
IMPLICIT NONE

!   I N T E G E R S
INTEGER maxa, mht, neq, nwk, ma, i, ithcall, ithsub

!    L O G I C A L S
LOGICAL  little, medium, full

!----------------------------------------------------------------------
!    p r o g r a m
!         to calculate addresses of diagonal elements in banded
!         matrix whose column heights are known.
!
!         mht  = active column heights.
!         maxa = addresses of diagonal elements.
!         neq  = number of equations.
!         nwk  = number of elements below the skyline.
!         ma   = maximum column height.
!
!           x x 0 0 0 0 0 0            1  3
!             x x 0 x 0 0 0               2  5    10
!               x 0 x 0 x 0                  4     9    18
!                 x x x x 0                     6  8 13 17
!                   x x x 0                        7 12 16
!                     x x x                          11 15 21
!                       x x                             14 20
!                         x                                19
!
!         Matrix with skyline             Storage Scheme
!
!   The matrix is assumed to be symmetric and to have zero elements
!   above the skyline.
!
!                   | 1|
!                   | 2|    MAXA gives the location of diagonal elements
!                   | 4|     in the matrix. MAXA must has neq+1 elements
!                   | 6|     where neq is the order of the matrix.
!          MAXA =   | 7|
!                   |11|
!                   |14|
!                   |19|
!                   |22|



DIMENSION maxa(*), mht(*)

!-----------------------------------------------------------------------
!  Determine printing options.
!-----------------------------------------------------------------------
DATA ithcall/0/, ithsub/0/
DATA little/.false./, medium/.false./, full/.false./
CALL popt ('adress  ', little, medium, full, ithcall, ithsub)

maxa(1)=1
maxa(2)=2
ma=0
IF(neq>1)THEN
   DO i=2,neq
      IF(mht(i)>ma) ma=mht(i)
      maxa(i+1)=maxa(i)+mht(i)+1
   ENDDO
ENDIF
nwk=maxa(neq+1)-maxa(1)

IF(medium)THEN
   WRITE(Out_unit,*)'neq,nwk = ',neq,nwk
   WRITE(Out_unit,*)'maxa array'
   WRITE(Out_unit,*)(maxa(i),i=1,neq+1)
   WRITE(Out_unit,*)' exiting addres'
ENDIF

RETURN
ENDSUBROUTINE addres
