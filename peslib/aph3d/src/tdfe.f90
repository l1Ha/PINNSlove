SUBROUTINE tdfe (id, theta, chi, mht, lm, yz, nd, neq, nume, nid, nidm, idi, ind, amatrx, nod, s, xx, beta, maxa, nwk)
USE todim_Module
USE var_Module
Use FileUnits_Module
!
!   $RCSfile: tdfe.f,v $ $Revision: 1.15 $
!   $Date: 89/11/16 11:03:33 $
!   $State: Stable $
!

!            P U R P O S E   O F    S U B R O U T I N E


!            I N P U T    A R G U M E N T S
!    id
!    theta
!    chi
!     mht       column heights.
!    lm
!    yz
!    nd
!    neq
!    nume
!    nid
!    nidm      maximum number of terms in the constraint equation.
!    idi
!    ind
!    amatrx
!    nod
!    s
!    xx
!    beta
!    maxa

IMPLICIT NONE
LOGICAL little, medium, full
INTEGER nwk
INTEGER ithcall, ithsub, id, lm, mht, nd, neq, nume, nid, nidm, idi, ind, nod, maxa(neq+1), n, m, i
REAL(Kind=WP_Kind) theta, chi, yz, amatrx, s, xx, beta(*), b
DIMENSION id(*), theta(*), chi(*), mht(neq+1), lm(18,nume), yz(18,nume), nid(*)
DIMENSION idi(nidm,*), amatrx(nwk), nod(9), s(378), xx(18), b(4,18)
DATA ithcall/0/, ithsub/0/
DATA little/.false./, medium/.false./, full/.false./

CALL popt ('tdfe    ', little, medium, full, ithcall, ithsub)
medium=.true.
!-----------------------------------------------------------------------
! This section reads in the element DATA.
!    yz      yz(2i+1,n) is the theta value for node n.
!            yz(2i+2,n) is the chi value for node n.
!    lm      node numbers.
!-----------------------------------------------------------------------
WRITE(Out_Unit,'("Inside tdfe: ind=",i5)')ind
IF(ind>0) GOTO 270
DO n=1,nume
   READ(Msher_Bin_Unit) m,(nod(i),i=1,9)
   IF(full) WRITE(Out_unit,'(A,/,(10I6))') m,(nod(i),i=1,9)
   IF(m/=n)THEN
      WRITE(Out_unit,*)'error n/=m',n,m
      STOP 'tdfe'
   ENDIF
   !-----------------------------------------------------------------------
   ! Store the coordinates of the nodal points
   !-----------------------------------------------------------------------
   DO i=1,9
      yz(2*i-1,n)=theta(nod(i))
      yz(2*i,n)=chi(nod(i))
   ENDDO
   !-----------------------------------------------------------------------
   ! Make sure that the coordinates of the side nodes are at the center
   ! of the corner nodes.
   !-----------------------------------------------------------------------
   CALL center (yz(1,n))
   !-----------------------------------------------------------------------
   ! Store the node numbers.
   !-----------------------------------------------------------------------
   DO i=1,9
      lm(2*i-1,n)=id(nod(i))
      lm(2*i,n)=0
   ENDDO
   CALL colht (mht, nd, lm(1,n), neq, nid, nidm, idi)
ENDDO
IF(medium)THEN
   WRITE(Out_unit,*)' yz array'
   WRITE(Out_unit,'(10ES11.3)')((yz(i,n),i=1,18),n=1,nume)
   WRITE(Out_unit,*)' lm arary'
   WRITE(Out_unit,'(20I6)')((lm(i,n),i=1,18),n=1,nume)
   WRITE(Out_unit,*)' column heights'
   WRITE(Out_unit,'(20I6)')(mht(i),i=1,neq)
   WRITE(Out_unit,*)' exiting tdfe part 1'
ENDIF
RETURN
!-----------------------------------------------------------------------
! This section calculates Hamiltonian matrix elements.
!-----------------------------------------------------------------------
  270 GOTO (280,330,330), ind
  !270 CONTINUE
  IF(IND==1)THEN
  280 CONTINUE
DO n=1,nume
   CALL modcor(yz(1,n))
   DO i=1,18
      xx(i)=yz(i,n)
   ENDDO
   DO i=1,378
      s(i)=0.d0
   ENDDO
   CALL quads (n, b, s, xx)
   CALL addban (amatrx, maxa, s, lm(1,n), nid, nidm, idi, beta, neq)
ENDDO
IF(little)THEN
   WRITE(Out_unit,*)'exiting tdfe part 2'
ENDIF
RETURN
ENDIF
!-----------------------------------------------------------------------
! This section calculates overlap matrix elements between the finite
! element basis.
!-----------------------------------------------------------------------
330 CONTINUE
DO n=1,nume
   CALL modcor(yz(1,n))
   DO i=1,18
      xx(i)=yz(i,n)
   ENDDO
   CALL quadm (n, b, s, xx)
   CALL addban (amatrx, maxa, s, lm(1,n), nid, nidm, idi, beta, neq)
ENDDO
IF(little)THEN
   WRITE(Out_unit,*)'exiting tdfe part 3'
ENDIF
RETURN
END
