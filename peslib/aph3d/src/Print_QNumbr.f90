SUBROUTINE Print_QNumbr(n, QNumbr, energy, xksq, wvec, NQuant, Location)
USE FileUnits_Asymptotic_Module
USE Numeric_Kinds_Module
USE Convrsns_Module
USE QAsy_Numbers_Module, ONLY: QLabels
IMPLICIT NONE
CHARACTER(LEN=*) location
INTEGER n, NQuant, i, j, QNumbr(NQuant,n)
CHARACTER(LEN=21), PARAMETER:: ProcName='Print_QNumbr'
CHARACTER(LEN=6) Print_Flag
REAL(Kind=WP_Kind) energy(n), xksq(n), wvec(n)

CALL PoptAsy(ProcName, Print_Flag)
WRITE(Out_Unit,*)'begin Print_QNumbr: Location='//TRIM(Location)
WRITE(Out_Unit,*)'Quantum Numbers (Jacobi)'
WRITE(Out_Unit,'(A4,4A16,9A6)')'i','energy(Ha)','energy(eV)','xksq','wvec', (TRIM(QLabels(j)),j=1,NQuant)
DO i = 1, n
   WRITE(Out_Unit,'(I4,4ES16.7,9(I5,1x))')i, energy(i), energy(i)*autoev, xksq(i), wvec(i), (QNumbr(j,i),j=1,NQuant)
ENDDO

WRITE(Out_Unit,*)'finish Print_QNumbr:'
RETURN
ENDSUBROUTINE Print_QNumbr
