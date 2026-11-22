SUBROUTINE NumericInformation
USE Numeric_Kinds_Module
USE FileUnits_Common_Module
IMPLICIT NONE
!     Author: Gregory A. Parker, Department of Physics and Astronomy, University of Oklahoma
!     Prints numeric kind information for integers and reals
INTEGER(KIND=I1_Kind) I1             ! 1 byte INTEGER
INTEGER(KIND=I2_Kind) I2             ! 2 byte INTEGER
INTEGER(KIND=I4_Kind) I4             ! 4 byte INTEGER
INTEGER(KIND=I8_Kind) I8             ! 8 byte INTEGER
REAL(KIND=SP_Kind) S              ! Single Precsion Floating Number
REAL(KIND=DP_Kind) D              ! Double Precsion Floating Number
REAL(KIND=QP_Kind) Q              ! Quadruple Precsion Floating Number

! 1 byte INTEGER
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Integer I1 Results'
WRITE(Out_Unit,*)'For I1_Kind=',I1_Kind, ' the computer results are:'
WRITE(Out_Unit,*)'BIT_SIZE=',BIT_SIZE(I1)
WRITE(Out_Unit,*)'DIGITS=',DIGITS(I1)
WRITE(Out_Unit,*)'HUGE',HUGE(I1)
WRITE(Out_Unit,*)'RADIX',RADIX(I1)
WRITE(Out_Unit,*)'RANGE',RANGE(I1)

! 2 byte INTEGER
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Integer I2 Results'
WRITE(Out_Unit,*)'For I2_Kind=',I2_Kind, ' the computer results are:'
WRITE(Out_Unit,*)'BIT_SIZE=',BIT_SIZE(I2)
WRITE(Out_Unit,*)'DIGITS=',DIGITS(I2)
WRITE(Out_Unit,*)'HUGE',HUGE(I2)
WRITE(Out_Unit,*)'RADIX',RADIX(I2)
WRITE(Out_Unit,*)'RANGE',RANGE(I2)

! 4 byte INTEGER
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Integer I4 Results'
WRITE(Out_Unit,*)'For I4_Kind=',I4_Kind, ' the computer results are:'
WRITE(Out_Unit,*)'BIT_SIZE=',BIT_SIZE(I4)
WRITE(Out_Unit,*)'DIGITS=',DIGITS(I4)
WRITE(Out_Unit,*)'HUGE',HUGE(I4)
WRITE(Out_Unit,*)'RADIX',RADIX(I4)
WRITE(Out_Unit,*)'RANGE',RANGE(I4)

! 8 byte INTEGER
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Integer I8 Results'
WRITE(Out_Unit,*)'For I8_Kind=',I8_Kind, ' the computer results are:'
WRITE(Out_Unit,*)'BIT_SIZE=',BIT_SIZE(I8)
WRITE(Out_Unit,*)'DIGITS=',DIGITS(I8)
WRITE(Out_Unit,*)'HUGE',HUGE(I8)
WRITE(Out_Unit,*)'RADIX',RADIX(I8)
WRITE(Out_Unit,*)'RANGE',RANGE(I8)

! Single Precsion Floating Number
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Single Precision Results'
WRITE(Out_Unit,*)'For SP_Kind=',SP_Kind, ' the computer results are:'
WRITE(Out_Unit,*)'DIGITS=',DIGITS(S)
WRITE(Out_Unit,*)'MAXEXPONENT',MAXEXPONENT(S)
WRITE(Out_Unit,*)'MINEXPONENT',MINEXPONENT(S)
WRITE(Out_Unit,*)'EPSILON',EPSILON(S)
WRITE(Out_Unit,*)'TINY',TINY(S)
WRITE(Out_Unit,*)'HUGE',HUGE(S)
WRITE(Out_Unit,*)'PRECISION',PRECISION(S)
WRITE(Out_Unit,*)'RADIX',RADIX(S)
WRITE(Out_Unit,*)'RANGE',RANGE(S)

! Double Precsion Floating Number
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Double Precision Results'
WRITE(Out_Unit,*)'For DP_Kind=',DP_Kind, ' the computer results are:'
WRITE(Out_Unit,*)'DIGITS=',DIGITS(D)
WRITE(Out_Unit,*)'MAXEXPONENT',MAXEXPONENT(D)
WRITE(Out_Unit,*)'MINEXPONENT',MINEXPONENT(D)
WRITE(Out_Unit,*)'EPSILON',EPSILON(D)
WRITE(Out_Unit,*)'TINY',TINY(D)
WRITE(Out_Unit,*)'HUGE',HUGE(D)
WRITE(Out_Unit,*)'PRECISION',PRECISION(D)
WRITE(Out_Unit,*)'RADIX',RADIX(D)
WRITE(Out_Unit,*)'RANGE',RANGE(D)

! Quadruple Precsion Floating Number
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Quadruple Precision Results'
WRITE(Out_Unit,*)'For QP_Kind=',QP_Kind, ' the computer results are:'
WRITE(Out_Unit,*)'DIGITS=',DIGITS(Q)
WRITE(Out_Unit,*)'MAXEXPONENT',MAXEXPONENT(Q)
WRITE(Out_Unit,*)'MINEXPONENT',MINEXPONENT(Q)
WRITE(Out_Unit,*)'EPSILON',EPSILON(Q)
WRITE(Out_Unit,*)'TINY',TINY(Q)
WRITE(Out_Unit,*)'HUGE',HUGE(Q)
WRITE(Out_Unit,*)'PRECISION',PRECISION(Q)
WRITE(Out_Unit,*)'RADIX',RADIX(Q)
WRITE(Out_Unit,*)'RANGE',RANGE(Q)

END SUBROUTINE NumericInformation