SUBROUTINE double (nel, chimax, th, ch, id, maxelp, idiv, lcor,  mapr, mxelmnt)
USE Numeric_Kinds_Module
USE FileUnits_Module
!
!     $RCSfile: double.f,v $ $Revision: 1.14 $
!     $Date: 89/10/18 14:17:59 $
!     $State: Stable $
!

!     P U R P O S E   O F    S U B R O U T I N E

!-----------------------------------------------------------------------
!     this routine doubles the range on chi and the number of elements.
!     it does a reflection about the chisym axis; only good for
!     symmetry=true.
!-----------------------------------------------------------------------

!     I N P U T    A R G U M E N T S
!
!
!     nel
!     chimax
!     th
!     ch
!     id(numnp)   id=-nce if the is constrained and positive otherwise
!     maxelp
!     idiv
!     lcor
!     mapr
!     mxelmnt


IMPLICIT NONE

!     I N T E G E R S
INTEGER id, idiv, iel, ielp, k, lcor, m1, m2, mapr, maxelp, mxelmnt,  nel


!     R E A L S
REAL(Kind=WP_Kind) c1, c2, ch, chimax, cmax, t1, t2, th, tmax,one

!     D I M E N S I O N S

DIMENSION th(9,mxelmnt), ch(9,mxelmnt), id(9,mxelmnt),  idiv(mxelmnt), lcor(2,mxelmnt), mapr(mxelmnt)

!     E X T E R N A L S

EXTERNAL swap


PARAMETER (one=1.d0)

IF(2*nel>mxelmnt)THEN
   WRITE(Out_unit,30) nel, mxelmnt
   STOP 'double'
ENDIF
ielp=nel
cmax=-one
tmax=-one
DO iel=1,nel
   ielp=ielp+1
   mapr(ielp)=iel
   th(1,ielp)=th(7,iel)
   th(2,ielp)=th(8,iel)
   th(3,ielp)=th(9,iel)
   th(4,ielp)=th(4,iel)
   th(5,ielp)=th(5,iel)
   th(6,ielp)=th(6,iel)
   th(7,ielp)=th(1,iel)
   th(8,ielp)=th(2,iel)
   th(9,ielp)=th(3,iel)
   ch(1,ielp)=2*chimax-ch(7,iel)
   ch(2,ielp)=2*chimax-ch(8,iel)
   ch(3,ielp)=2*chimax-ch(9,iel)
   ch(4,ielp)=2*chimax-ch(4,iel)
   ch(5,ielp)=2*chimax-ch(5,iel)
   ch(6,ielp)=2*chimax-ch(6,iel)
   ch(7,ielp)=2*chimax-ch(1,iel)
   ch(8,ielp)=2*chimax-ch(2,iel)
   ch(9,ielp)=2*chimax-ch(3,iel)
   IF(ch(9,ielp)>cmax.or.th(9,ielp)>tmax)THEN
      maxelp=ielp
      tmax=th(9,ielp)+1.d-10
      cmax=ch(9,ielp)+1.d-10
   ENDIF
ENDDO
nel=2*nel
DO k=1,9
   c1=ch(k,maxelp)
   t1=th(k,maxelp)
   c2=ch(k,nel)
   t2=th(k,nel)
   ch(k,nel)=c1
   th(k,nel)=t1
   ch(k,maxelp)=c2
   th(k,maxelp)=t2
ENDDO
CALL swap (nel, maxelp, lcor, idiv, id)
m1=mapr(maxelp)
m2=mapr(nel)
mapr(nel)=m1
mapr(maxelp)=m2
RETURN
!
 30   FORMAT(/,'***error***:2*nel>mxelmnt:',2i10)
END
