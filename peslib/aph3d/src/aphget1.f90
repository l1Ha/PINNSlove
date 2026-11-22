SUBROUTINE aphget1(nchnls)
USE Numeric_Kinds_Module
USE FileUnits_Module
!
!     $RCSfile: aphget.f,v $ $Revision: 1.14 $
!     $Date: 89/11/28 13:27:09 $
!     $State: Stable $
!
USE nstate_Module
USE qlam_Module
USE aphchl_Module
USE totj_Module
IMPLICIT NONE
INTEGER ithcll, ithsub, iaph,itee, nchnls
CHARACTER(LEN=5) mlabel


!#include <nstate.parms>
!#include <Inc.qlam>
!#include <Inc.aphchl>
!#include <Inc.totj>

LOGICAL little, medium, full
DATA ithcll/0/, ithsub/0/
DATA little/.false./, medium/.false./, full/.false./
CALL popt('aphget  ', little, medium, full, ithcll, ithsub)
IF(medium)WRITE(Out_Unit,3)
3    FORMAT(///, ' aph quantum numbers', ///)
DO 1 iaph=1, naph
    IF(iaph == 1)THEN
    itee=1
    ELSEIF(lambda(iaph-1) == lambda(iaph))THEN
    itee=itee+1
    ELSE
    itee=1
    ENDIF
    tee(iaph)=itee
    IF(medium)WRITE(Out_Unit,2)iaph, tee(iaph), lambda(iaph)
1 CONTINUE
2    FORMAT(/, ' iaph=', i5, '  t=', i5, '  lambda=', i5)
nchnls=naph

OPEN(Unit=APH_Quant_Bin_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/APH_Quant.bin',form='unformatted', status='unknown')
    mlabel='aph  '
    WRITE(APH_Quant_Bin_Unit)mlabel
    WRITE(APH_Quant_Bin_Unit)jtot
    WRITE(APH_Quant_Bin_Unit)naph
    WRITE(APH_Quant_Bin_Unit)(tee(iaph),iaph=1,naph)
    WRITE(APH_Quant_Bin_Unit)(lambda(iaph),iaph=1,naph)
CLOSE(Unit=APH_Quant_Bin_Unit)

RETURN
ENDSUBROUTINE aphget1
