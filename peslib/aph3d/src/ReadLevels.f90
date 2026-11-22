SUBROUTINE ReadLevels(neig,nrho,elevels,rhovals)
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
INTEGER neig, nrho
CHARACTER(LEN=5) readx, ready
CHARACTER(LEN=8) cont
LOGICAL there
INTEGER irho, nr, ienergy, iunit, ounit, neng
REAL(Kind=WP_Kind) rhovals(nrho), elevels(neig,nrho)

WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Reading energy levels'
irho = 1
iunit=21
ounit=91

INQUIRE (file=OutDIR(1:LEN(TRIM(OutDIR)))//'GraphicsOut/Sfelevl_All.rbw', exist=there)
IF(there)THEN
   OPEN(Unit=iunit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'GraphicsOut/Sfelevl_All.rbw',status='old')
ELSE
   WRITE(*,*)'sfelevl must exist'
   WRITE(Out_Unit,*)'sfelevl must exist'
   WRITE(Out_Unit,*)'Stopping in ReadLevels1'
   STOP 'ReadLevels1'
ENDIF

DO irho=1,NRho
   READ(iunit,8)readx
   IF(readx/='readx')THEN
      WRITE(*,*)'readx is not correct'
      WRITE(Out_Unit,*)'readx is not correct'
      WRITE(Out_Unit,*)'Stopping in ReadLevels2'
      STOP 'ReadLevels2'
   ENDIF

   READ(iunit,9)nr
   IF(nr/=1)THEN
      WRITE(*,*)'nr/=1; nr=',nr
      WRITE(Out_Unit,*)'nr/=1; nr=',nr
      WRITE(Out_Unit,*)'Stopping in ReadLevels3'
      STOP 'ReadLevels3'
   ENDIF

   READ(iunit,3)rhovals(irho)
   IF(irho==1.or.irho==nrho)THEN
      WRITE(Out_Unit,*)'irho,rhovals(irho)',irho,rhovals(irho)
   ENDIF
   READ(iunit,8)ready

   IF(ready/='ready')THEN
      WRITE(*,*)'ready is not correct'
      WRITE(Out_Unit,*)'ready is not correct'
      WRITE(Out_Unit,*)'Stopping in ReadLevels4'
      STOP 'ReadLevels4'
   ENDIF

   READ(iunit,*)neng
   IF(neig/=neng)THEN
      WRITE(*,*)'neig/=neng',neig,neng
      WRITE(Out_Unit,*)'neig/=neng','neig=',neig,'  neng=',neng
      WRITE(Out_Unit,*)'Stopping in ReadLevels5'
      !STOP 'ReadLevels5'
   ENDIF

   READ(iunit,*)(elevels(ienergy,irho),ienergy=1,MIN(neig,neng))
   READ(iunit,4)cont
ENDDO

WRITE(Out_Unit,*)'RhoVals  ELevels(eV)'
DO IRho=1,NRho
   WRITE(Out_Unit,'(11es15.6)')RhoVals(Irho),(Elevels(IEnergy,IRho),IEnergy=1,10)
ENDDO

WRITE(Out_Unit,*)'completed reading energy levels'
WRITE(Out_Unit,*)

8 FORMAT(2x,a5)
9 FORMAT(i6)
3 FORMAT(1x,f10.5)
5 FORMAT(1x,e14.7)
4 FORMAT(2x,a8)

RETURN
ENDSUBROUTINE ReadLevels
