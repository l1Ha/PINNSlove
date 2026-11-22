SUBROUTINE FREEID(ID)
SAVE
!
!   $RCSfile: freeid.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:10 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE -
!
!         RELEASES VEC(ID) FROM W.
!
!     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
!     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
!     THE FUNCTION OF THE ROUTINE.
!
!
! **********************************************************************
!
INTEGER   ACTIVE, ID, MV, MXNEW, MXOLD, MXRST, NIL,  SCPX, SOLCPX, V, WAD, X
COMMON   /STLMID/ NIL, MV, V, MXNEW, MXOLD, MXRST, SCPX, SOLCPX, X
COMMON   /STLMWH/ WAD(2), ACTIVE(2)
!
IF(ID == NIL) GOTO 9999
IF(ID == ACTIVE(1)) ACTIVE(1) = NIL
IF(ID == ACTIVE(2)) ACTIVE(2) = NIL
!
9999  RETURN
END
