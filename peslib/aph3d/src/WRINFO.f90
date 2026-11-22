SUBROUTINE WRINFO(NOR, LOC, W, IW)
USE Numeric_Kinds_Module
!
!   $RCSfile: wrinfo.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:52 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE -
!
!         DRIVER ROUTINE FOR THE ROUTINES THAT WRITE INFORMATION.
!
!
!     INPUT PARAMETERS -
!
!         NOR = NUMBER OF THE CALLING ROUTINE.
!         LOC = LOCATION IN THAT ROUTINE.
!
!     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
!     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
!     THE FUNCTION OF THE ROUTINE.
!
!
! **********************************************************************
!
REAL(Kind=WP_Kind)     W(1)
INTEGER   ALPHA, BETA, BETAPI, COPT, IW(1), LAMBDA, LOC, MSGLVL, MXREST, NOR, NU, PMAX, POINTR, POPT, S, SCR, NERR, NOUT
COMMON   /STLMPL/ PMAX, POPT, COPT, MXREST
COMMON   /STLMPR/ MSGLVL, NERR  , NOUT
COMMON   /STLMPV/ ALPHA, BETA, BETAPI, LAMBDA, NU, POINTR, S, SCR
!
IF(MSGLVL <= 1) GOTO 9999
!
IF(.NOT. (MSGLVL == 3 .OR. (MSGLVL == 2 .AND. NOR == 20 .AND. (LOC == 1 .OR. LOC == 5)))) GOTO 10
!
  IF(NOR == 20 .OR. NOR == 2 .OR. NOR == 19) CALL WRINF3(NOR, LOC, W(LAMBDA), IW(POINTR))
!
10    IF(MSGLVL == 4) CALL WRINF4(NOR, LOC, W(ALPHA), W(BETA), W(BETAPI),  W(LAMBDA), W(NU), IW(POINTR), W(S), PMAX)
!
9999  RETURN
END
