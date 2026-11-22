SUBROUTINE Regons
USE Numeric_Kinds_Module
USE FileUnits_Module
USE Regins_Module
!

!  P U R P O S E   O F    S U B R O U T I N E

!--------------------------------------------------------------------
!  This routine sets up the propagation regions.
!
!  The following schematic may be useful.
!
!
!  Log-derivative              VIVS
!  --------------------- -------------------------
!  |         |           |            |            |
!  start      endaph     SwitchToVIVS       enddelve      finish
!  --------- ------------------------ ------------
!  APH            Delves             Jacobi
!
!  distance -----
!--------------------------------------------------------------------
IMPLICIT NONE
!--------------------------------------------------------------------
!  Find the actual stopping distance.
!--------------------------------------------------------------------

REAL(Kind=WP_Kind) max, min

IF(enddelve < endaph)enddelve=endaph
IF(finish < enddelve .or. finish < endaph)finish=max(endaph, enddelve)
!--------------------------------------------------------------------
!  Propagation Range.
!--------------------------------------------------------------------
 10   WRITE(Out_Unit,1)start, finish
    1 FORMAT(1x, 'The propagation will start at: ', e14.7,  /1x, 'The propagation will END at: ', e14.7)
!--------------------------------------------------------------------
!  Log-derivative propagation region.
!--------------------------------------------------------------------
IF(start < SwitchToVIVS)THEN
   WRITE(Out_Unit,2)start, min(SwitchToVIVS, finish)
 2       FORMAT(1x, 'The Log-derivative method will be used from ', e14.7, '  to  ', e14.7)
ELSE
   WRITE(Out_Unit,*)'The Log-derivative method will not be used.'
ENDIF
!--------------------------------------------------------------------
!  VIVS propagation range.
!--------------------------------------------------------------------
IF(SwitchToVIVS < finish)THEN
   WRITE(Out_Unit,3)max(start, SwitchToVIVS), finish
 3       FORMAT(1x, 'The VIVS method will be used from ', e14.7, '  to  ', e14.7)
ELSE
   WRITE(Out_Unit,*)'The VIVS method will not be used'
ENDIF
!--------------------------------------------------------------------
!  APH region
!--------------------------------------------------------------------
IF(start < endaph)THEN
   WRITE(Out_Unit,4)start, min(endaph, finish)
 4       FORMAT(1x, 'The APH region ranges from ', e14.7, '  to  ', e14.7)
ELSE
   WRITE(Out_Unit,*)'NO propagation in the APH region.'
ENDIF
!--------------------------------------------------------------------
!  Delves region.
!--------------------------------------------------------------------
IF(endaph < min(enddelve, finish))THEN
   WRITE(Out_Unit,5)max(start, endaph), min(enddelve, finish)
 5       FORMAT(1x, 'The Delves region ranges from ',  e14.7, '  to  ', e14.7)
ELSE
   WRITE(Out_Unit,*)'NO propagation in the Delves region.'
ENDIF
!--------------------------------------------------------------------
!  Jacobi region.
!--------------------------------------------------------------------
IF(max(endaph, enddelve) < finish)THEN
   WRITE(Out_Unit,6)max(endaph, enddelve), finish
 6       FORMAT(1x, 'The Jacobi region ranges from ', e14.7, '  to  ', e14.7)
ELSE
   WRITE(Out_Unit,*)'NO propagation in the Jacobi region.'
ENDIF
RETURN
!-----------------------***END-regons***-------------------------------
ENDSUBROUTINE Regons
