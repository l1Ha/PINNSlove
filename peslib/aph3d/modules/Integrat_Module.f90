MODULE Integrat_Module
   USE Narran_Module
   SAVE
   LOGICAL :: integrat(narran)=.True.      ! integrat  Set true for each arrangement channel that is accessible.
   NAMELIST/IntChanl/ integrat
END MODULE Integrat_Module
