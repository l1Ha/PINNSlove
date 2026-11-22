SUBROUTINE epsi_graph (xvals, psi, nxval, nchanl, evals, filename)
USE Numeric_Kinds_Module
USE FileUnits_Module
USE DiatomicPot_Module, ONLY: DiatomicPot
IMPLICIT NONE
CHARACTER(LEN=*) :: filename
INTEGER :: k, j, nxval, nchanl, ntime, ifirst, ilast, itime
INTEGER, PARAMETER :: ngroup = 5
REAL(Kind=WP_Kind) :: xvals (nxval), psi (nxval, nchanl), evals (nchanl), psi_max, psi_abs, ak1, ak2

REAL(Kind=WP_Kind), ALLOCATABLE::ak (:)

ALLOCATE (ak (nchanl) )
ak2 = abs (evals (nchanl) - evals (1) ) / 10
IF(nchanl>1)THEN
   ak1 = abs (evals (2) - evals (1) ) / 3.d0
ELSE
   ak1 = ak2
ENDIF
ak (1) = max (ak1, ak2)
IF(nchanl>1)THEN
   ak1 = abs (evals (nchanl) - evals (nchanl - 1) ) / 3.d0
ELSE
   ak1 = ak2
ENDIF
ak (nchanl) = max (ak1, ak2)
DO k = 2, nchanl - 1
ak1 = min (abs (evals (k) - evals (k - 1) ), abs (evals (k + 1) - evals (k) ) ) / 3.d0
ak (k) = max (ak1, ak2)

ENDDO
DO k = 1, nchanl
psi_max = 0.d0
DO j = 1, nxval
psi_abs = abs (psi (j, k) )
IF(psi_max<psi_abs) psi_max = psi_abs
ENDDO
ak (k) = ak (k) / psi_max

ENDDO
OPEN(Unit=EPsi_Graph_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//filename, status = 'unknown', form = 'formatted')
ntime = (nchanl - 1) / ngroup + 1
WRITE(EPsi_Graph_Unit,'(7(A15,","))')" ","xvals ","Shifted Wave 1","Shifted Wave 2","Shifted Wave 3","Shifted 4","Shifted 5"
DO itime = 1, ntime
ifirst = (itime-1) * ngroup + 1
ilast = min (ifirst + ngroup - 1, nchanl)
DO j = 1, nxval
   WRITE(EPsi_Graph_Unit,'(1x,A15,",",11(1pe15.7,","))') " ",xvals (j), (evals (k) + ak (k) * psi (j, k),k = ifirst, ilast)
ENDDO
IF(itime/=ntime) WRITE(EPsi_Graph_Unit, * ) '&'
ENDDO
CLOSE(unit=EPsi_Graph_Unit)

OPEN(Unit=Psi_Bin_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Solution.bin', status = 'unknown', form = 'unformatted')
WRITE(Psi_Bin_Unit) xvals
WRITE(Psi_Bin_Unit) psi
CLOSE(Psi_Bin_Unit)

DEALLOCATE (ak)
RETURN
ENDSUBROUTINE epsi_graph
