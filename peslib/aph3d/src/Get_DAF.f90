SUBROUTINE Get_DAF(nxval, xval, delta, mval, dafx, hermit, ndaf, nd_daf, toe, usys2, kmax, Identifier)
USE Numeric_Kinds_Module
USE FileUnits_Module
USE DAFSigma_Module
USE Numbers_Module
USE DiatomicPot_Module, ONLY: DiatomicPot
!-----------------------------------------------------------------------
! This routine was written by G. A. Parker
! If you find an error or have an improvement please send a messge to
! parker@ou.edu 
!-----------------------------------------------------------------------
IMPLICIT NONE
CHARACTER(LEN=15) :: Identifier
CHARACTER(LEN=65) :: filename
INTEGER :: mval, jhalf, ndaf, kmax, idaf, karg, nd_daf, nxval
REAL(KIND=WP_Kind) :: dafx(0:nd_daf,0:ndaf), hermit(0:mval+ndaf+1), toe(0:nd_daf)
REAL(KIND=WP_Kind) :: delta, usys2, aterm, bterm, cterm, sigma_ratio, xval (nxval)
WRITE(Out_Unit, * )
WRITE(Out_Unit, * ) 'Begin(Get_Daf)'
!----------------------------------------------------------------------
! This section gets a guess for the best sigma. One could refine
! this guess but is probably not very useful.
!-----------------------------------------------------------------------
aterm=10.4893_WP_Kind
bterm=-3.5095_WP_Kind
cterm=-0.4744_WP_Kind-mval
sigma_ratio=(-bterm+sqrt(bterm**2-Four*aterm*cterm))/(Two*aterm)
sigma=sigma_ratio*delta
jhalf=nd_daf/2+1
kmax=(jhalf-1)
WRITE(Out_Unit, * ) 'sigma_ratio=', sigma_ratio
WRITE(Out_Unit,*)'sigma=',sigma,'  delta=',delta

!----------------------------------------------------------------------- 
! calculate the daf deltas.
!-----------------------------------------------------------------------
kmax=MAX(kmax,nd_daf)
dafx=Zero
CALL DAF(sigma, delta, mval, dafx, kmax, hermit, ndaf, nd_daf)

!----------------------------------------------------------------------- 
! output the daf deltas for subsequent plotting.
!-----------------------------------------------------------------------
DO idaf = 0, ndaf
  IF(idaf==0) filename = 'GraphicsOut/dafx(0)_graph'//TRIM(Identifier)//TRIM(DiatomicPot)//'.csv'
  IF(idaf==1) filename = 'GraphicsOut/dafx(1)_graph'//TRIM(Identifier)//TRIM(DiatomicPot)//'.csv'
  IF(idaf==2) filename = 'GraphicsOut/dafx(2)_graph'//TRIM(Identifier)//TRIM(DiatomicPot)//'.csv'
  WRITE(Out_Unit, * ) 'Creating dafx(idaf)_graph.csv for idaf=', idaf
  CALL Graf_DAF(delta, delta, dafx (0, idaf), nd_daf, filename, 0, idaf)
ENDDO

!WRITE(DAF_Unit,*)
!WRITE(DAF_Unit,*)'Mu=',usys2/Two,'  kmax=',kmax,'  nd_daf=',nd_daf
DO karg=0,MAX(kmax,nd_daf-1)
IF(karg<=kmax)THEN
      toe(karg)=-dafx(karg,2)/usys2
   ELSE
      toe(karg)=Zero
   ENDIF
ENDDO

WRITE(Out_Unit, * ) 'END(Get_Daf)'

WRITE(Out_Unit, * )
RETURN
ENDSUBROUTINE Get_DAF
