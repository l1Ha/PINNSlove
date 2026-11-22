MODULE GaussB_Module
   USE Numeric_Kinds_Module
   USE Parms_Module
   USE Pow_Module
   IMPLICIT NONE
   SAVE
   LOGICAL :: old_way = .True.
   INTEGER(KIND=IW_Kind) :: noscil(narran)    ! Number of harmonic oscillators.
   INTEGER(KIND=IW_Kind) :: nhermt(narran)   = 25          ! Number of Gauss_Hermite quadrature points.
   INTEGER(KIND=IW_Kind) :: nglegn(narran)   = 46          ! Number of Gauss_Legendre quadrature points.
   INTEGER(KIND=IW_Kind) :: nlegndre(narran) = 6           ! Number of Legendre polynomials.
   INTEGER(KIND=IW_Kind) :: intwt(narran)    = 1           !
   REAL(Kind=WP_Kind) :: weau(narran)        = 2.00534D-02 ! vibrational constant in Hartree's.
   REAL(Kind=WP_Kind) :: ralpha(narran)      = 3.166D0     ! ralpha*alpha is the
   REAL(Kind=WP_Kind) :: re(narran)          = 1.40112d0   ! equilibrium positions for diatomic fragments.
   REAL(Kind=WP_Kind) :: rx(narran)          = 1.085       ! rx*re is the position of harmonic oscillator basis. This
                                                           ! variable should usually be between 1.0 and 1.2 and 1.1 is
                                                           ! usually the best.
   REAL(Kind=WP_Kind) :: balpha(narran)      = 1.0D0       !
   REAL(Kind=WP_Kind) :: calpha(narran)      = 0.92D0      !
   REAL(Kind=WP_Kind) :: dalpha(narran)      = 0.03D0      !
   REAL(Kind=WP_Kind) :: anharm(narran)      = 0.8D0       !
   REAL(Kind=WP_Kind) :: wexeau(narran)      = 5.52847D-04 ! anharmonicity constant in Hartree's.
   REAL(Kind=WP_Kind) :: delta(narran)       = 0.01D0      !
   REAL(Kind=WP_Kind) :: zeta(narran)        =1.d0         !
   NAMELIST/gauss/ rx, re, noscil, nhermt, npow, zeta, delta, nglegn, weau, wexeau, ralpha,       &
                         balpha, calpha, dalpha, anharm, intwt, old_way, nlegndre  !   nlegndre(iarran)=jmax(0,iarran)
END MODULE GaussB_Module
