MODULE Quantb_Module
   USE Numeric_Kinds_Module
   USE Parms_Module
   IMPLICIT NONE
   SAVE 
   INTEGER(KIND=IW_Kind) jminval, jmaxval
   INTEGER(KIND=IW_Kind) numinval, numaxval
   INTEGER(KIND=IW_Kind) minvib(narran)                    ! minimum vibrational quantum number for each arrangement channel.
   INTEGER(KIND=IW_Kind) maxvib(narran)                    ! maxvib    maximum vibrational quantum number for each arrangement channel.
   INTEGER(KIND=IW_Kind) jskip(narran)                     ! 2 for channel 1 if symmetry.  Otherwise, 1.
   INTEGER(KIND=IW_Kind) jmin(0:MxVib, narran, 0:MxVib)  ! minimum rotational state for each arrangement channel, each vibrational state
   INTEGER(KIND=IW_Kind) jmax(0:MxVib, narran)            ! maximum rotational state for each arrangement channel and each vibrational
   INTEGER(KIND=IW_Kind), ALLOCATABLE:: chanl(:,:)         ! chanl(mxbasis, 0:mxmega)  channel quantum number for each basis function.
   INTEGER(KIND=IW_Kind), ALLOCATABLE::  nvib(:)           ! nvib(mxbasis) vibrational quantum number for each channel.
   INTEGER(KIND=IW_Kind), ALLOCATABLE::  jrot(:,:)         ! jrot(mxbasis,0:mxmega)    rotational quantum number for each basis fcn.
   NAMELIST/quantum/ minvib, maxvib, jmin, jmax
END MODULE Quantb_Module
