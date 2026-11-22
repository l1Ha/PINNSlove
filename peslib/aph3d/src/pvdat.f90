SUBROUTINE pvdat(rhovals, evalues, nenergy, nrho, ounit)
USE FileUnits_Module
IMPLICIT NONE
INTEGER nenergy, nrho, ounit, ienergy, irho
REAL  rhovals(nrho), evalues(nenergy,nrho)

WRITE(ounit,5)
WRITE(ounit,3)rhovals
WRITE(ounit,9)
WRITE(ounit,10)

DO ienergy = 1, nenergy
   IF(ienergy<=9)THEN
      WRITE(ounit,17)ienergy
   ELSEIF(ienergy<=99)THEN
      WRITE(ounit,27)ienergy
   ELSEIF(ienergy<=999)THEN
      WRITE(ounit,37)ienergy
   ELSE
      WRITE(Out_Unit,*)'Error in routine pvdat'
      STOP 'pvdat'
   ENDIF
   WRITE(ounit,4)(evalues(ienergy,irho),irho=1,nrho)
   WRITE(ounit,9)
ENDDO


  3 FORMAT(1x,f10.5)
  4 FORMAT(1x,e14.7)
  5 FORMAT('Enter independent DATA rho')
 17 FORMAT('Enter dependent DATA E',i1)
 27 FORMAT('Enter dependent DATA E',i2)
 37 FORMAT('Enter dependent DATA E',i3)
  9 FORMAT('END')
 10 FORMAT('/*')

RETURN
ENDSUBROUTINE pvdat
