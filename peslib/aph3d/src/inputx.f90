SUBROUTINE inputx(id, theta, chi, numnp, neq)
Use FileUnits_Module
USE eigif_Module
!
!   $RCSfile: inputx.f,v $ $Revision: 1.15 $
!   $Date: 89/11/16 11:02:58 $
!   $State: Stable $
!

!            P U R P O S E   O F    S U B R O U T I N E


!            I N P U T    A R G U M E N T S
!  id
!  theta
!  chi
!  nunmp
!  neq


IMPLICIT NONE
LOGICAL little, medium, full
INTEGER idt, ithcall, ithsub, numnp, neq, n, id(1), idum
REAL(Kind=WP_Kind) theta(1), chi(1)

DATA ithcall/0/, ithsub/0/
DATA little/.false./, medium/.false./, full/.false./

CALL popt ('inputx  ', little, medium, full, ithcall, ithsub)
DO WHILE (N/=numnp)
   READ(Msher_Bin_Unit) n, idt, theta(n), chi(n)
   IF(full) WRITE(Out_unit,*) 'inputx: ',n, idt, theta(n), chi(n)
   IF(n<=0.or.n>numnp)THEN
   WRITE(Out_unit,*)'n,numnp = ', n,numnp
   STOP 'inputx'
   CONTINUE
   ENDIF
   id(n)=idt
ENDDO
neq=0
DO n=1,numnp
   idum=id(n)+3
   IF(IDUM==3)THEN
      neq=neq+1
      id(n)=neq
   ELSEIF(Idum==2.or.Idum==4)THEN
      id(n)=0
   ENDIF
ENDDO

RETURN
END
