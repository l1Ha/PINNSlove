SUBROUTINE REPL
USE Numeric_Kinds_Module
!
!   $RCSfile: repl.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:36 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE -
!
!         UPDATES CERTAIN CONTROL VARIABLES.
!
!     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
!     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
!     THE FUNCTION OF THE ROUTINE.
!
!
! **********************************************************************
!
REAL(Kind=WP_Kind)  ALTMU, MU, NEXTMU, OLDMU
INTEGER   BUFF, CNEG, CPOS, ITERNO, MV, MXNEW, MXOLD, MXRST, N, NIL, OLCPOS, P, REST, RNEW, ROLD, SCPX, SOLCPX, TCONV, V, X
LOGICAL   USEMX, ZERBET
COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW, ROLD, REST, P, USEMX, ZERBET
COMMON   /STLMID/ NIL, MV, V, MXNEW, MXOLD, MXRST, SCPX, SOLCPX, X
COMMON   /STLMMU/ MU, OLDMU, NEXTMU, ALTMU
!
OLDMU = MU
MU = NEXTMU
ROLD = RNEW
OLCPOS = CPOS
!     ***********************
!     INTERCHANGE MX-BUFFERS.
!     ***********************
BUFF = MXNEW
MXNEW = MXOLD
MXOLD = BUFF
SOLCPX = SCPX
!
RETURN
END
