SUBROUTINE wgdat(rhovals, evalues, nenergy, nrho, ounit)
USE FileUnits_Module
IMPLICIT NONE
INTEGER nenergy, nrho, ounit, ienergy, irho, jumpmin, jump, min, max
REAL rhovals(nrho), evalues(nenergy,nrho)

WRITE(ounit,5)
WRITE(ounit,6)nrho
WRITE(ounit,3)rhovals
DO 1 ienergy = 1, nenergy
   WRITE(ounit,7)
   WRITE(ounit,6)nrho
   WRITE(ounit,4)(evalues(ienergy,irho),irho=1,nrho)
   WRITE(ounit,8)
    1 CONTINUE
CLOSE(Unit=ounit)

OPEN(Unit=ounit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'GraphicsOut/Adiabatic.graph',status='unknown')
jumpmin = min(10,nenergy)
DO jump = jumpmin, nenergy, 10
   DO irho = 1, nrho
      WRITE(ounit,9)rhovals(irho),  (evalues(ienergy,irho),ienergy=max(jump-9,1),jump)
   ENDDO
   IF(jump/=nenergy)WRITE(ounit,*)'&'
ENDDO
CLOSE(ounit)

OPEN(Unit=ounit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'GraphicsOut/Adiabatic.csv',status='unknown')
DO irho=1,nrho
   WRITE(ounit,811)rhovals(irho),(evalues(ienergy,irho),ienergy=1,nenergy)
ENDDO

  811 FORMAT(1x,e14.7,10000(",",e14.7))
    3 FORMAT(1x,f10.5)
    4 FORMAT(1x,e14.7)
    5 FORMAT(' $readx',5x,'distance in bohr')
    6 FORMAT(1x,i5)
    7 FORMAT(' $ready',5x,'energy in ev.')
    8 FORMAT(' $CONTINUE')
    9 FORMAT(11e15.7)
RETURN
ENDSUBROUTINE wgdat
