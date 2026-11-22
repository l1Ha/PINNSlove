MODULE Physical_Constants_Module
   USE Numeric_Kinds_Module
   USE Prefix_Module
   REAL(Kind=WP_Kind), PARAMETER:: hbar      = 1.05457266d-34     ! hbar with units erg-sec
   REAL(Kind=WP_Kind), PARAMETER:: hartree   = 4.3597482d-18      ! conversion factor from erg to Hartree
   REAL(Kind=WP_Kind), PARAMETER:: time_conv = hbar/hartree/femto ! time-unit conversion
   REAL(Kind=WP_Kind), PARAMETER:: mhbar=1.d0*time_conv           ! Lifetime and Time-Delay prefactor
END MODULE Physical_Constants_Module
