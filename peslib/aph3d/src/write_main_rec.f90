subroutine WRITE_main_rec(distance)
!subroutine WRITE_main_rec(unit,distance,sector)
USE FileNames_Module
USE Masses_Module
USE convrsns_Module
USE recovery_Module
USE jtot_Module
USE parity_Module
USE jacobi_Module
USE quanta_Module
USE quantum_numbers_Module
USE nrho_Module
!USE symmetry_Module
IMPLICIT NONE
logical debug
integer sector,unit, jtotal_max
CHARACTER(LEN=40) dummy
REAL(Kind=WP_Kind) dat,distance

debug=.true.
if(unit/=Out_Unit.AND.unit/=msg_unit)THEN
   if(debug) WRITE(Msg_Unit,*) 'Writing main.rec'
   OPEN(unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//work_dir//mainrec_name)
   rewind(unit)
ENDIF
WRITE(unit,'(a7)') '!System'
CALL Write_real(unit,'a mass.......................>',mass(1)/amutoau)
CALL Write_real(unit,'b mass.......................>',mass(2)/amutoau)
CALL Write_real(unit,'c mass.......................>',mass(3)/amutoau)
CALL Write_int(unit,'Starting Total J.............>',jtot)
CALL Write_int(unit,'Ending Total J...............>',jtotal_max)
CALL Write_int(unit,'parity.......................>',parity)
!CALL Write_log(unit,'symmetry.....................>',symmetry)
CALL Write_int(unit,'llambda max .................>',llambda_max)
CALL Write_real(unit,'Starting rho.................>',distance)
CALL Write_real(unit,'Finishing rho................>',rho_finish)
CALL Write_int(unit,'Number of rhos ..............>',nrho-sector+1)
CALL Write_int(unit,'First sector ................>',sector)
CALL Write_real(unit,'rho boundary.................>',rho_boundary)
!
! WRITE down the method for BC
!
dummy = 'Delves                    '
if(jacobi) dummy = 'Jacobi               '
CALL Write_str(unit,'BC: Delves or Jacobi.........>',dummy)
OPEN(Unit=Method_data_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//'method.data')
read(Method_data_Unit,*) dummy
close(Method_data_Unit)
CALL Write_str(unit,'Method.......................>',dummy)
CALL Write_int(unit,'nstep........................>',1)
CALL Write_int(unit,'ksqrtmax.....................>',1)
CALL Write_str(unit,'Energy data..................>',energies_data)
CALL Write_str(unit,'Channels data................>',channels_data)
WRITE(unit,'(a3)') 'end'
call flush(unit)
if(unit/=Out_Unit.AND.unit/=msg_unit)THEN
   close(unit)
ENDIF

return
end
