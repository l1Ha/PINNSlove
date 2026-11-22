SUBROUTINE psi_graph (xvals, psi, nxval, nchanl, filename)
USE Numeric_Kinds_Module
USE FileUnits_Module
USE DiatomicPot_Module, ONLY: DiatomicPot
IMPLICIT NONE
CHARACTER(LEN=*) :: filename
INTEGER :: k, j, nxval, nchanl, ntime, ifirst, ilast, itime, ngroup
PARAMETER (ngroup = 5)

REAL(Kind=WP_Kind) :: xvals (nxval), psi (nxval, nchanl)
OPEN(Unit=Psi_Graph_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//filename, status = 'unknown', form = 'formatted')
WRITE(Psi_Graph_Unit,'(7(A13,","))')" ","xvals ","Wave 1","Wave 2","Wave 3","Wave 4","Wave 5"
ntime = (nchanl - 1) / ngroup + 1
DO itime = 1, ntime
   ifirst = (itime-1) * ngroup + 1
   ilast = min (ifirst + ngroup - 1, nchanl)
   DO j = 1, nxval
      WRITE(Psi_Graph_Unit, '(1x,A13,",",11(1pe13.5,","))')" ", xvals (j), (psi (j, k), k = ifirst, ilast)
   ENDDO
   IF(itime/=ntime) WRITE(Psi_Graph_Unit, * ) '&'
ENDDO
CLOSE(unit=EPsi_Graph_Unit)

CLOSE(unit=Sol_Bin_Unit)
OPEN(Unit=Sol_Bin_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Solution.bin', status = 'unknown', form = 'unformatted')
WRITE(Sol_Bin_Unit) xvals
WRITE(Sol_Bin_Unit) psi
CLOSE(Sol_Bin_Unit)

RETURN
ENDSUBROUTINE psi_graph
