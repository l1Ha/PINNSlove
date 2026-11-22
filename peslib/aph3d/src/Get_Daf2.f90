SUBROUTINE get_daf2 (nxval, xval, delta, mval, dafx, hermit, ndaf, nd_daf, toe, usys2, kmax, Identifier)
USE FileUnits_Module
USE DAFSigma_Module
USE DiatomicPot_Module, ONLY: DiatomicPot
!TO get t_chir1
!-----------------------------------------------------------------------
! This routine was written by G. A. Parker
! IF you find an error or have an improvement please send a messge to
! Parker@ou.edu
!-----------------------------------------------------------------------
IMPLICIT NONE
CHARACTER(LEN=65) :: filename
CHARACTER(LEN=15) :: Identifier
INTEGER :: mval, jhalf, ndaf, kmax, idaf, karg, nd_daf, nxval
REAL(Kind=WP_Kind) :: dafx (0:nd_daf, 0:ndaf), hermit (0:mval + ndaf + 1), &
 delta, usys2, aterm, bterm, cterm, sigma_ratio, toe (0: &
 nd_daf), xval (nxval)
WRITE(Out_Unit, * )
WRITE(Out_Unit, * ) 'Begin(Get_DAF2)'
!-----------------------------------------------------------------------
! This section gets a guess for the best sigma. One could refine
! this guess but is probably not very useful.
!-----------------------------------------------------------------------
aterm = 10.4893D0
bterm = - 3.5095D0
cterm = - 0.4744D0 - mval
sigma_ratio = ( - bterm + sqrt (bterm**2 - 4.d0 * aterm * cterm) ) / (2.d0 * aterm)
sigma = sigma_ratio * delta
jhalf = nd_daf / 2 + 1
kmax = (jhalf - 1)
WRITE(Out_Unit, * ) 'sigma_ratio=', sigma_ratio

WRITE(Out_Unit, * ) 'sigma,delta=', sigma, delta
!-----------------------------------------------------------------------
! calculate the daf deltas.
!-----------------------------------------------------------------------
kmax = max (kmax, nd_daf)
dafx = 0.d0

CALL DAF(sigma, delta, mval, dafx, kmax, hermit, ndaf, nd_daf)
!-----------------------------------------------------------------------
! output the daf deltas for subsequent plotting.
!-----------------------------------------------------------------------
DO idaf = 0, ndaf
IF(idaf==0) filename = 'GraphicsOut/dafx(0)_graph'//TRIM(Identifier)//TRIM(DiatomicPot)//'.csv'
IF(idaf==1) filename = 'GraphicsOut/dafx(1)_graph'//TRIM(Identifier)//TRIM(DiatomicPot)//'.csv'
IF(idaf==2) filename = 'GraphicsOut/dafx(2)_graph'//TRIM(Identifier)//TRIM(DiatomicPot)//'.csv'
WRITE(Out_Unit, * ) 'Creating dafx(idaf)_graph.csv for idaf=', idaf
!CALL graf (delta, delta, dafx (0, idaf), nd_daf, filename, 0, idaf)

ENDDO
WRITE(Out_Unit, * ) 'usys2,kmax,nd_daf=', usys2, kmax, nd_daf
DO karg = 0, max (kmax, nd_daf - 1)
IF(karg<=kmax)THEN
!!THIS IS THE DIFFERENCE
   toe (karg) = - dafx (karg, 1) / usys2
ELSE
   toe (karg) = 0.d0
ENDIF
ENDDO

   99 FORMAT(1x,2ES13.5)
WRITE(Out_Unit, * ) 'END(Get_DAF2)'

WRITE(Out_Unit, * )
RETURN
ENDSUBROUTINE Get_DAF2
