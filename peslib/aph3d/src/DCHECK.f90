SUBROUTINE DCHECK(D, N, LEN, LOC)
SAVE
!
!   $RCSfile: dcheck.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:07 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE - (VER = 1)
!
!         THIS ROUTINE CHECKS THE COLUMN LENGTHS.
!
!     INPUT PARAMETERS -
!
!         D     = POINTER VECTOR TO DIAGONAL ELEMENTS OF K.
!         N     = DIMENSION OF M MATRIX.
!
!
!     OUTPUT PARAMETERS -
!
!         LEN = 0, NO ERROR HAS BEEN DETECTED.
!               6, A TOO SHORT, OR A TOO LONG COLUMN HAS BEEN DETECTED.
!         LOC = THE INDEX OF THE NUMBER OF THE COLUMN.
!
! **********************************************************************
!
INTEGER   D(1), DI, DIM1, I, LEN, LENG, LOC, N
!
DIM1 = 0
DO 10 I = 1, N
  DI = D(I)
  LENG = DI - DIM1
  DIM1 = DI
!
  IF(.NOT. (LENG <= 0 .OR. LENG > I)) GOTO 10
    LEN = 6
    LOC = I
    GOTO 8888
!
10 Continue
!
8888  RETURN
END
