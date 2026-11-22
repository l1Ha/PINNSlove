SUBROUTINE INITAD(NUM, BASEID, ADR, N, ADDRSS)
!
!   $RCSfile: initad.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:14 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE - (VER = 1 OR 2)
!
!         INITIALIZES THE ADDRESS VECTOR.
!
!     INPUT PARAMETERS -
!
!         NUM    = THE NUMBER OF VECTORS IN ONE IDENTIFIER GROUP.
!         BASEID = THE IDENTIFIER OF THE FIRST VECTOR IN THE GROUP - 1.
!         ADR    = THE ADDRESS IN THE W-ARRAY.
!
!     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
!     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
!     THE FUNCTION OF THE ROUTINE.
!
!
! **********************************************************************
!
INTEGER   ADDRSS(1), ADR, BASEID, I, ID, N, NUM
!
IF(NUM == 0) GOTO 9999
DO 10 I = 1, NUM
  ID = BASEID + I
  ADDRSS(ID) = ADR
  ADR = ADR - N
10      CONTINUE
!
9999  RETURN
END
