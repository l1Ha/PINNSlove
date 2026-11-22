FUNCTION RANFX(dummy)
USE Numeric_Kinds_Module
!
!   $RCSfile: RANFX.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:35 $
!   $State: Stable $
!
!
!     This function uses Alliant fortran SUBROUTINEs.
!
!     This version created 880126 by Jan Linderberg.
!
!     Generates random numbers with a psuedo-normal distribution.
!
INTEGER dummy
REAL(Kind=WP_Kind) avrge, RANFX
REAL(KIND=WP_Kind) rand
LOGICAL first
EXTERNAL rand
DATA first/.true./,avrge/1.5/
SAVE first
!
IF(first)THEN
   RANFX = rand(3)
   first = .false.
ENDIF
!
RANFX = rand(0) + rand(0) + rand(0) - avrge
RETURN
END
