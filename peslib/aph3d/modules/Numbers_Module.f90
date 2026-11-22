MODULE Numbers_Module
USE Numeric_Kinds_Module
IMPLICIT NONE
SAVE
! Numbers  ===============================================================================
INTEGER(KIND=IW_Kind), PARAMETER :: ione      = 1
INTEGER(KIND=IW_Kind), PARAMETER :: itwo      = 2
INTEGER(KIND=IW_Kind), PARAMETER :: ithree    = 3
INTEGER(KIND=IW_Kind), PARAMETER :: ifour     = 4

REAL(Kind=WP_Kind), PARAMETER :: zero   =    0.0_WP_kind  !  0
REAL(Kind=WP_Kind), PARAMETER :: one    =    1.0_WP_kind  !  1
REAL(Kind=WP_Kind), PARAMETER :: two    =    2.0_WP_kind  !  2
REAL(Kind=WP_Kind), PARAMETER :: three  =    3.0_WP_kind  !  3
REAL(Kind=WP_Kind), PARAMETER :: four   =    4.0_WP_kind  !  4
REAL(Kind=WP_Kind), PARAMETER :: five   =    5.0_WP_kind  !  5
REAL(Kind=WP_Kind), PARAMETER :: six    =    6.0_WP_kind  !  6
REAL(Kind=WP_Kind), PARAMETER :: seven  =    7.0_WP_kind  !  7
REAL(Kind=WP_Kind), PARAMETER :: eight  =    8.0_WP_kind  !  8
REAL(Kind=WP_Kind), PARAMETER :: nine   =    9.0_WP_kind  !  9
REAL(Kind=WP_Kind), PARAMETER :: ten    =   10.0_WP_kind  ! 10
REAL(Kind=WP_Kind), PARAMETER :: eleven =   11.0_WP_kind  ! 11
REAL(Kind=WP_Kind), PARAMETER :: twelve =   12.0_WP_kind  ! 12
REAL(Kind=WP_Kind), PARAMETER :: thirteen = 13.0_WP_kind  ! 13
REAL(Kind=WP_Kind), PARAMETER :: fourteen = 14.0_WP_kind  ! 14
REAL(Kind=WP_Kind), PARAMETER :: fifteen  = 15.0_WP_kind  ! 15
REAL(Kind=WP_Kind), PARAMETER :: Sixteen  = 16.0_WP_kind  ! 16
REAL(Kind=WP_Kind), PARAMETER :: thirty   = 30.0_WP_kind  ! 30
REAL(Kind=WP_Kind), PARAMETER :: fortyfive =45.0_WP_kind  ! 45
REAL(Kind=WP_Kind), PARAMETER :: fortyeight=48.0_WP_kind  ! 48

! Fractions
REAL(KIND=WP_Kind), PARAMETER :: quarter     = 0.25_WP_kind              !  1/4
REAL(Kind=WP_Kind), PARAMETER :: onethird    = 1.00_WP_kind/3.0_WP_kind  !  1/3
REAL(Kind=WP_Kind), PARAMETER :: half        = 0.50_WP_kind              !  1/2
REAL(Kind=WP_Kind), PARAMETER :: ThreeHalves = 1.50_WP_kind              !  3/2

! Common Reals
REAL(Kind=WP_Kind), PARAMETER :: rt2=SQRT(two)       ! Sqrt(2)
REAL(Kind=WP_Kind), PARAMETER :: rt3=SQRT(three)     ! Sqrt(3)
REAL(Kind=WP_Kind), PARAMETER :: Pi=four*ATAN(one)   ! Pi
REAL(Kind=WP_Kind), PARAMETER :: FourthPi = Pi/Four  ! One fourth of pi.
REAL(Kind=WP_Kind), PARAMETER :: HalfPi   = Pi/Two   ! One half of pi.
REAL(Kind=WP_Kind), PARAMETER :: TwoPi    = Two*Pi   ! Two times pi.

! Numerics ===============================================================================
REAL(Kind=WP_Kind), PARAMETER :: WP_Small  = TINY(one)    ! smallest number
REAL(Kind=WP_Kind), PARAMETER :: WP_Large  = HUGE(one)    ! largest number
REAL(Kind=WP_Kind), PARAMETER :: WP_Epsilon= EPSILON(one) ! smallest difference
INTEGER, PARAMETER :: WP_Prec   = PRECISION(one)          ! number of digits of precision
INTEGER, PARAMETER :: WP_minexp = MINEXPONENT(one)        ! smallest exponent
INTEGER, PARAMETER :: WP_maxexp = MAXEXPONENT(one)        ! largest exponent
END MODULE Numbers_Module
