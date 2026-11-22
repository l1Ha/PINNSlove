SUBROUTINE CMPRSS(POINTR, PMAX)
SAVE
!
!   $RCSfile: cmprss.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:05 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE -
!
!         UPDATES POINTR, I.E. DELETES ELEMENTS POINTING TO FALSE
!         EIGENPAIRS.
!
!     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
!     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
!     THE FUNCTION OF THE ROUTINE.
!
!
! **********************************************************************
!
INTEGER   I, ISAVE, J, LEN, NHITS, PMAX, POINTR(1)
!
LEN = POINTR(PMAX + 1)
IF(LEN <= 0) GOTO 9999
!     ********************************************
!     A FALSE ELEMENT IS <= 0. LOCATE THE FIRST.
!     ********************************************
DO 10 I = 1, LEN
  IF(POINTR(I) <= 0) GOTO 20
10      CONTINUE
GOTO 9999
!
20    NHITS = 0
ISAVE = I
!     ****************************
!     WRITE OVER THE FALSE VALUES.
!     ****************************
DO 40 J = ISAVE, LEN
  IF(.NOT. POINTR(J) <= 0) GOTO 30
    NHITS = NHITS + 1
    GOTO 40
!
30      POINTR(I) = POINTR(J)
  I = I + 1
!
40      CONTINUE
!
!     ******************
!     UPDATE THE LENGTH.
!     ******************
POINTR(PMAX + 1) = LEN - NHITS
!
9999  RETURN
END
