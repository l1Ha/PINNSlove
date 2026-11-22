MODULE Convrsns_Module
   USE Numeric_Kinds_Module
   IMPLICIT NONE
   SAVE
   REAL(Kind=WP_Kind), PARAMETER :: autoev=27.2113961d0   ! Multiply by this constant to convert from Hartree's to electron volts.
   REAL(Kind=WP_Kind) amutoau  ! Multiply by this constant to convert from atomic mass units to electron mass units (atomic units).
   REAL(Kind=WP_Kind) angtoau  ! Multiply by this constant to convert from angstroms to Hartree's.
   REAL(Kind=WP_Kind) radtodeg ! Multiply by this constant to convert from radians to degrees.
   REAL(Kind=WP_Kind) cmm1toau ! Multiply by this constant to convert from inverse centimeters to Hartree's.
END MODULE Convrsns_Module
