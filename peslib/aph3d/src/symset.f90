SUBROUTINE symset(anmomchi, tdvfb, nchi, debug, sym)
USE FileUnits_Module
USE DVR_Module
USE DVRMasses_Module
USE OneD_Module
USE DVR2_Module
USE Parms_Module
IMPLICIT NONE
!
!     symset calculates the dvr points (in chi) and the dvr
!     j^2 matrix (in antotchi) for the representation
!     designated by sym. anmomchi and tdvfb must be
!     declared EXTERNAL in the calling program
!
INTEGER NChi
CHARACTER(LEN=3) sym
LOGICAL debug
EXTERNAL anmomchi, tdvfb
IF(Nchi>NchiMax)THEN
   WRITE(*,*)"Ntheta or Nchi exceeds maximum allowed values"
   WRITE(*,*)"Nchi=",Nchi,"  NchiMax=",NchiMax
   STOP "SymSet"
ENDIF
!
!     calculate dvr points and transformation matrix
!
CALL tdvfb(tchi, ptchi, whchi, nchi, nchimax, debug)
!
!     calculate the dvr j^2 matrix in chi for later use
!
IF((sym/='e  ').AND.(sym/='ep '))THEN
   CALL anmomchi(tchi, ptchi, nchi, nchimax, antotchi, chij, debug)
ELSE
   CALL anmomchi(tchi, ptchi, nchi, nchimax, antotchi, chij1, debug)
ENDIF
RETURN
END
