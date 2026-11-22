SUBROUTINE Quant_Out (n, kchan, elect, mvib, jrot, lorb, energy, xksq, location)
USE FileUnits_Module
USE Numeric_Kinds_Module
USE Masses_Module
USE TotalEng_Module
IMPLICIT NONE
CHARACTER(LEN=*) location
INTEGER n, i, kchan(n), elect(n), mvib(n), jrot(n), lorb(n)
REAL(Kind=WP_Kind) energy(n), xksq(n)
WRITE(Out_Unit,'(A)')'begin Quant_Out: Location='//TRIM(Location)
WRITE(Out_Unit,*)'Quantum Numbers'
WRITE(Out_Unit,'(8A)')'    i','  kchan',' elect','  mvib','  jrot','  lorb','          energy','            xksq'
DO i = 1, n
   WRITE(Out_Unit,'(6I6,2ES16.7)')i, kchan(i), elect(i), mvib(i), jrot(i), lorb(i), energy(i), xksq(i)
   IF(ABS(xksq(i)-usys2*(Etot-energy(i)))/ABS(xksq(i))>1.d-4)THEN
      WRITE(Out_Unit,'(16x,A,3ES16.7)')'ERROR in Quant_Out: ', xksq(i), usys2*(Etot-energy(i)), ABS(xksq(i)-usys2*(Etot-energy(i)))/ABS(xksq(i))
      !STOP 'Quant_Out'
   ENDIF
ENDDO

WRITE(Out_Unit,'(2A)')'finish Quant_Out: Location=',TRIM(Location)
RETURN
ENDSUBROUTINE Quant_Out
