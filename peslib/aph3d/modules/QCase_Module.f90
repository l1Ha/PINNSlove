MODULE QCase_Module
   USE Numeric_Kinds_Module
   USE Parms_Module
   IMPLICIT NONE
   SAVE
   LOGICAL :: qcase=  .True.
   INTEGER(KIND=IW_Kind) ngch                   !
   INTEGER(KIND=IW_Kind) ngth                   !
   REAL(Kind=WP_Kind) :: qswitch= 2.1D0                !
   REAL(Kind=WP_Kind), ALLOCATABLE:: xch(:)  ! Chi abscissas
   REAL(Kind=WP_Kind), ALLOCATABLE:: wch(:)  ! Chi quadrature weights
   REAL(Kind=WP_Kind), ALLOCATABLE:: xth(:)  ! Theta abscissas
   REAL(Kind=WP_Kind), ALLOCATABLE:: wth(:)  ! Theta quadrature weights
   NAMELIST/qcasep/ qcase, qswitch
END MODULE QCase_Module
