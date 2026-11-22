MODULE Numeric_Kinds_Module
!=========================================================================================
! This module defines types and constants to be used throughout other codes. Here is 
! defined the real type wp_kind (working precision) that defines the real type of all
! declared real variables. 
!=========================================================================================
IMPLICIT NONE
SAVE
! Define compiler-dependent type constants as parameters
INTEGER, PARAMETER :: I1_kind = SELECTED_INT_KIND(1)  !  1 byte INTEGER
INTEGER, PARAMETER :: I2_kind = SELECTED_INT_KIND(4)  !  4 byte INTEGER
INTEGER, PARAMETER :: I4_kind = SELECTED_INT_KIND(7)  !  8 byte INTEGER
INTEGER, PARAMETER :: I8_kind = SELECTED_INT_KIND(14) ! 16 byte INTEGER
INTEGER, PARAMETER :: SP_kind = SELECTED_REAL_KIND(p=6,r=30)    ! Single Precision
INTEGER, PARAMETER :: DP_kind = SELECTED_REAL_KIND(p=14,r=300)  ! Double Precision
INTEGER, PARAMETER :: QP_kind = SELECTED_REAL_KIND(p=30,r=4000) ! Quadruple Precision
! note above real kind types also work for complex declarations

! Select real kind to be used throughout the rest of code
INTEGER, PARAMETER :: WP_Kind = DP_kind
INTEGER, PARAMETER :: IW_Kind = I4_kind

END MODULE Numeric_Kinds_Module