MODULE Storage_Module
   SAVE
   LOGICAL :: irdindep = .false.  ! If true reads energy independent information from disk.
   LOGICAL :: iwrindep = .true.   ! If true WRITEs energy independent information to disk.
   NAMELIST/storage/ irdindep, iwrindep
END MODULE Storage_Module
