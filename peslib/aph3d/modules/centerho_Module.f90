MODULE centerho_Module
!---------------------------------------------------------------------
! Inc.centerho
!
!   rhocnt
!---------------------------------------------------------------------
USE Numeric_Kinds_Module
   IMPLICIT NONE
   SAVE
   REAL(Kind=WP_Kind)rhocnt
   COMMON / centerho / rhocnt
END MODULE centerho_Module
