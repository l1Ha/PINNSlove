SUBROUTINE repset(pointg, jeven, mega, parity, nsym)
USE FileUnits_Module


!            P U R P O S E   O F    S U B R O U T I N E
!  Determine the irreducible representation in C2v or Cs for the FBR
!   which defines the DVR.


!            I N P U T    A R G U M E N T S
!    pointg
!    jeven
!    mega
!    parity

!            O U T P U T    A R G U M E N T S
!    nsym


IMPLICIT NONE

!            L O G I C A L S

LOGICAL jeven

!            I N T E G E R S

INTEGER  mega, parity, jinit, iarg, irefl, nsym

!            C H A R A C T E R S

CHARACTER(LEN=3) pointg

!            C O M M O N S

!            I N T R I N S I C    F U N C T I O N S
!            E X T E R N A L S

!  ---------------------------------------------------------------------
IF(pointg == 'c2v')THEN
!  ---------------------------------------------------------------------
!     II.  Reflection symmetry (C2v)
!          jinit is 0 or 1 as initial j is even or odd.
!          irefl is result of reflection about chi=0.
!  --------------------------------------------------------------
jinit=1
IF(jeven)jinit=0
!     iarg=(jinit-mega)
!     irefl=(-1)**(iarg)
iarg=ABS((jinit-mega))
irefl=mod(iarg,2)
!  -----------------------------------------------------------------
!          A.  Even Parity
!             i.  Even reflection symmetry in both reflections (A1)
!    Case 1   parity=0 irefl=0  (internal z-axis is out of the page)

!             chi=pi/2
!             +  |  +
!    chi=pi  ____|____  chi=0      A1 representation
!                |                (Whitnell and Light)
!             +  |  +
!            chi=3pi/2

!  -------------------------------------------------------------------
           IF(parity==0.AND.irefl==0)THEN
             nsym = 1
             WRITE(Out_Unit,*)'nsym = 1 (A1):  FBR = {cos[(2n)chi]}'
           ENDIF
!  ---------------------------------------------------------------------
!             ii.  Odd reflection symmetry in both reflections (A2)
!    Case 2   parity=0 irefl=1 (internal z-axis is out of the page)

!             chi=pi/2
!             -  |  +
!    chi=pi  ____|____  chi=0      B1 representation
!                |                (Whitnell and Light)
!             +  |  -
!            chi=3pi/2

!  ---------------------------------------------------------------------
            IF(parity==0.AND.irefl==1)THEN
              nsym = 3
              WRITE(Out_Unit,*)'nsym = 3 (B1):  FBR = {sin[(2n)chi]}'
            ENDIF
!  ---------------------------------------------------------------------
!          B.  Odd Parity
!             i.  Even reflection symmetry about pi/2, odd about 0 (B2)
!    Case 3   parity=1 irefl=1 (internal z-axis is out of the page)

!             chi=pi/2
!             +  |  +
!    chi=pi  ____|____  chi=0      A2 representation
!                |                (Whitnell and Light)
!             -  |  -
!            chi=3pi/2

!  ---------------------------------------------------------------------
            IF(parity==1.AND.irefl==1)THEN
              nsym = 2
              WRITE(Out_Unit,*)'nsym = 2 (A2):  FBR = {sin[(2n+1)chi]}'
            ENDIF
!  ---------------------------------------------------------------------
!             ii. Odd reflection symmetry about pi/2, even about 0 (B1)
!    Case 4   parity=1 irefl=0 (internal z-axis is out of the page)

!             chi=pi/2
!             -  |  +
!    chi=pi  ____|____  chi=0     B2 representation
!                |               (Whitnell and Light)
!             -  |  +
!            chi=3pi/2

!  ---------------------------------------------------------------------
            IF(parity==1.AND.irefl==0)THEN
              nsym = 4
              WRITE(Out_Unit,*)'nsym = 4 (B2):  FBR = {cos[(2n+1)chi]}'
            ENDIF
!  ---------------------------------------------------------------------
ELSE
!  ---------------------------------------------------------------------
!     II.  No reflection symmetry (Cs)
!  ---------------------------------------------------------------------
!             i. even about 0
!    Case 5   parity=0 (internal z-axis is out of the page)

!             chi=pi/2
!                +
!    chi=pi  _________  chi=0     A(prime) representation
!                                   (Whitnell and Light)
!                +
!            chi=3pi/2

!  ---------------------------------------------------------------------
            IF(parity==0)THEN
              nsym = 5
!                   WRITE(Out_Unit,*)'nsym = 5 (A(prime)):',
              WRITE(Out_Unit,*)'nsym = 5 (A(prime) and A(prime)(prime):',&
                    'FBR = {cos[(n)chi]} and {sin[(n)chi]}', &
                    'FBR = {cos[(2n)chi]} and {sin[(2n)chi]}'
            ENDIF
!  ---------------------------------------------------------------------
!             ii. odd about 0
!    Case 6   parity=1 (internal z-axis is out of the page)

!             chi=pi/2
!                +
!    chi=pi  _________  chi=0     A(prime)(prime) representation
!                                   (Whitnell and Light)
!                -
!            chi=3pi/2

!  ---------------------------------------------------------------------
            IF(parity==1)THEN
              nsym = 6
!                   WRITE(Out_Unit,*)'nsym = 6 (A(prime)(prime)):',
              WRITE(Out_Unit,*)'nsym = 5 (A(prime) and A(prime)(prime):', &
                   'FBR = {cos[(n)chi]} and {sin[(n)chi]}',  &
                   'FBR = {cos[(2n+1)chi]} and {sin[(2n+1)chi]}'
            ENDIF
!  ---------------------------------------------------------------------
!  ---------------------------------------------------------------------
ENDIF
!
RETURN
END
