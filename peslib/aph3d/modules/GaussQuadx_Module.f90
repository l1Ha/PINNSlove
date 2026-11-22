MODULE GaussQuadx_Module
!-----------------------------------------------------------------------
!  noscil    Number of oscillators in each arrangement channel.
!  nhermt    Number of Gauss_Hermite quadrature points.
!  nlegndre  Number of Legendre polynomials in each arrangement channel.
!  nglegn    Number of Gauss_Legendre quadrature points.
!  alpha     Harmonic oscillator exponential factor.
!  ralpha    Usually set to 0.95 to make the harmonic oscillator cover
!            a larger range. This will improve the rate of convergence
!            for the upper vibrational levels.
!  re        Vibrational equilibrium distance.
!  rx        Usually set to 1.1 which shifts the harmonic basis to
!            larger distances. This improves the rate of convergence
!            for upper vibrational states of anharmonic potentials.
!-----------------------------------------------------------------------
   USE Narran_Module
   USE Numeric_Kinds_Module
   USE pow_Module
   ! PORT FIX (2026-10-02): this legacy module used to declare its own copies of
   ! noscil/nhermt/nglegn/nlegndre/ralpha/re/rx/weau, which were NEVER filled by
   ! the $gauss namelist (that binds GaussB_Module).  Routines that import only
   ! GaussQuadx (e.g. input.f90) therefore read uninitialised memory -> NaN
   ! (wzero=NaN in the Jacobi/Delves basis stage).  Share GaussB's variables.
   USE GaussB_Module, ONLY: noscil, nhermt, nglegn, nlegndre, ralpha, re, rx, weau
   SAVE
   REAL(Kind=WP_Kind) alpha(narran)
   REAL(Kind=WP_Kind) ,ALLOCATABLE :: xpth(:,:),wpth(:,:)
!      NAMELIST/gaussq/ alpha, ralpha, re, rx, xpth, wpth, weau
END MODULE GaussQuadx_Module
