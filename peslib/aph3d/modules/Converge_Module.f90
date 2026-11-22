MODULE Converge_Module
   USE Numeric_Kinds_Module
   USE Parms_Module, Only: MaxDim
   IMPLICIT NONE
   SAVE
   INTEGER(KIND=IW_Kind) :: ngood=1            ! Number of good eigenvectors desired.
   INTEGER(KIND=IW_Kind) :: nave1=1            ! Range for eigenav.  nave1 < nave2
   INTEGER(KIND=IW_Kind) :: nave2=1            ! Range for eigenav.  nave1 < nave2
   INTEGER(KIND=IW_Kind) :: nfreq(0:MaxDim)=-1 ! nfreq(0:mxmega) is used MxMega < MaxDim
   REAL(Kind=WP_Kind)    :: eigmin=1.0D-03     ! Minimum eigenvalue
   REAL(Kind=WP_Kind)    :: ovrerr=1.0D-04     ! Overlap completeness error
   NAMELIST/convrg/ eigmin, ovrerr, ngood, nfreq, nave1, nave2
END MODULE Converge_Module
