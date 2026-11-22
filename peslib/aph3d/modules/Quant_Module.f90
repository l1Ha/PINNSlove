MODULE Quant_Module
   USE Narran_Module
   USE Parms_Module
   SAVE
   INTEGER(KIND=IW_Kind) minvib(narran)         ! Minimum vibrational state for each channel.
   INTEGER(KIND=IW_Kind) maxvib(narran)         ! Maximum vibrational state for each channel.
   INTEGER(KIND=IW_Kind) jmin(0:MxVib,narran)  ! Minimum rotational state for each vibrational state in each channel.
   INTEGER(KIND=IW_Kind) jmax(0:MxVib,narran)  ! Maximum rotational state for each vibrational state in each channel.
   NAMELIST/quant/ minvib, maxvib, jmin, jmax
END MODULE Quant_Module
