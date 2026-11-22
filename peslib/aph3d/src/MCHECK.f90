SUBROUTINE MCHECK(M, D, N, DIAGM, MEQI, LEN, LOC)
USE Numeric_Kinds_Module
!
!   $RCSfile: mcheck.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:23 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE - (VER = 1)
!
!         THIS ROUTINE CHECKS THE DIAGONAL ELEMENTS IN M.
!
!     INPUT PARAMETERS -
!
!         D     = POINTER VECTOR TO DIAGONAL ELEMENTS OF M (IF
!                 M IS NOT DIAGONAL).
!         DIAGM = PROFIL <= 2
!         M     = M MATRIX STORED ACCORDING TO THE USERGUIDE.
!         MEQI  = PROFIL == 1
!         N     = DIMENSION OF M MATRIX.
!         X     = VECTOR TO BE MULTIPLIED WITH M.
!
!
!     OUTPUT PARAMETERS -
!
!         LEN  = 0, NO ERROR HAS BEEN DETECTED.
!                7, A NEGATIVE DIAGONAL ELEMENT HAS BEEN DETECTED.
!         LOC  = IS THE LOCATION OF THE NEGATIVE ELEMENT.
!
! **********************************************************************
!
REAL(Kind=WP_Kind)  M(1)
INTEGER   D(1), DI, I, LEN, LOC, N
LOGICAL   DIAGM, MEQI
!
LEN = 0
IF(MEQI) GOTO 8888
!
IF(.NOT. DIAGM) GOTO 30
  DO 20 I = 1, N
    IF(.NOT. M(I) < 0.0D0) GOTO 20
      LEN = 7
      LOC = I
      GOTO 8888
20        CONTINUE
!
  GOTO 8888
!
30    DO 40 I = 1, N
  DI = D(I)
  IF(.NOT. M(DI) < 0.0D0) GOTO 40
    LEN = 7
    LOC = I
    GOTO 8888
40      CONTINUE
!
8888  RETURN
END
