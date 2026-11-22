SUBROUTINE aphget(nchnls)
!
!this routine is called by:
!         aph3d
!this routine calls
!         popt
!-----------------------------------------------------------------------

USE FileUnits_Module
USE Narran_Module
USE qlam_Module
USE aphchl_Module
USE totj_Module
USE APHchl_Module
IMPLICIT NONE

LOGICAL little, medium, full
INTEGER totj,iaph,nchnls,ithsub,ithcll
DATA ithcll/0/, ithsub/0/
DATA little/.false./, medium/.false./, full/.false./
CALL popt('aphget  ', little, medium, full, ithcll, ithsub)
!      medium = .true.
IF(medium)WRITE(Out_Unit,3)
 3    FORMAT(//, ' aph quantum numbers', //)

totj=jtot
WRITE(Out_Unit,*)'totj,naph=',totj,naph
DO iaph=1,naph
   tee(iaph)=iaph
   lambda(iaph)=0
ENDDO
!IF(naph/=-1)STOP

IF(naph==-1)THEN !This section of code is never executed
OPEN(Unit=APH_Quant_Bin_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/APH_Quant.bin',form='unformatted', status='unknown')
READ(APH_Quant_Bin_Unit)
READ(APH_Quant_Bin_Unit)totj
IF(totj/=jtot)THEN
   WRITE(Msg_Unit,*)' totj/=jtot: totj, jtot=',totj,jtot
   STOP 'aphget'
ENDIF
READ(APH_Quant_Bin_Unit)naph
IF(naph>nstate)THEN
   WRITE(Msg_Unit,*)' naph > nstate: naph, nstate=',naph,nstate
   STOP 'aphget'
ENDIF
WRITE(Out_Unit,*)'-----------------------------------------------------------------------------'
READ(APH_Quant_Bin_Unit) (tee(iaph),iaph=1,naph)
READ(APH_Quant_Bin_Unit) (lambda(iaph),iaph=1,naph)
CLOSE(Unit=APH_Quant_Bin_Unit)
ENDIF  !End of dummy section of code

IF(medium)THEN
   DO iaph=1, naph
      WRITE(Out_Unit,'(1x, " iaph=", i5, "  t=", i5, "  lambda=", i5)')iaph, tee(iaph), lambda(iaph)
   ENDDO
ENDIF

nchnls=naph
RETURN
ENDSUBROUTINE aphget
