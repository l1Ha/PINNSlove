MODULE MemAdd_Module
   USE Numeric_Kinds_Module
   SAVE
   INTEGER(KIND=IW_Kind) memlast     ! last storage location
   INTEGER(KIND=IW_Kind) memloc(100) ! memory storage locations.
   INTEGER(KIND=IW_Kind) nmemory     ! number of memory storage locations used.
   REAL(Kind=WP_Kind) precious    ! variable used to see array length have be exceeded.
   NAMELIST/memadd/ memlast, memloc, nmemory, precious
END MODULE MemAdd_Module
