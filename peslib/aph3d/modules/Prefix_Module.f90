MODULE Prefix_Module
   USE Numeric_Kinds_Module
   IMPLICIT NONE
   SAVE
   REAL(Kind=WP_Kind), PARAMETER:: exa  =1.d+18  ! Prefix constant
   REAL(Kind=WP_Kind), PARAMETER:: peta =1.d+15  ! Prefix constant
   REAL(Kind=WP_Kind), PARAMETER:: tera =1.d+12  ! Prefix constant
   REAL(Kind=WP_Kind), PARAMETER:: giga =1.d+9   ! Prefix constant
   REAL(Kind=WP_Kind), PARAMETER:: megaf=1.d+6   ! Prefix constant
   REAL(Kind=WP_Kind), PARAMETER:: kilo =1.d+3   ! Prefix constant
   REAL(Kind=WP_Kind), PARAMETER:: hecto=1.d+2   ! Prefix constant
   REAL(Kind=WP_Kind), PARAMETER:: deca =1.d+1   ! Prefix constant
   REAL(Kind=WP_Kind), PARAMETER:: deci =1.d-1   ! Prefix constant
   REAL(Kind=WP_Kind), PARAMETER:: centi=1.d-2   ! Prefix constant
   REAL(Kind=WP_Kind), PARAMETER:: milli=1.d-3   ! Prefix constant
   REAL(Kind=WP_Kind), PARAMETER:: micro=1.d-6   ! Prefix constant
   REAL(Kind=WP_Kind), PARAMETER:: nano =1.d-9   ! Prefix constant
   REAL(Kind=WP_Kind), PARAMETER:: pico =1.d-12  ! Prefix constant
   REAL(Kind=WP_Kind), PARAMETER:: femto=1.d-15  ! Prefix constant
   REAL(Kind=WP_Kind), PARAMETER:: atto =1.d-18  ! Prefix constant
END MODULE Prefix_Module
