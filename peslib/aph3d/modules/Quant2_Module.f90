MODULE Quant2_Module
   USE Narran_Module
   USE Parms_Module
   SAVE
! minvib(narran)         Minimum vibrational state for each channel.
! maxvib(narran)         Maximum vibrational state for each channel.
! jmin(0:MxVib,narran) Minimum rotational state for each vibrational state in each channel.
! jmax(0:MxVib,narran) Maximum rotational state for each vibrational state in each channel.
INTEGER minvib(narran), maxvib(narran), jmin(0:MxVib,narran), jmax(0:MxVib,narran)
END MODULE Quant2_Module
