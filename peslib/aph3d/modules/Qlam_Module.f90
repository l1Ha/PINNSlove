MODULE Qlam_Module
   USE Narran_Module
   USE NState_Module
   SAVE
!  INTEGER,ALLOCATABLE :: lambda(:), tee(:)
   INTEGER(KIND=IW_Kind) lambda(nstate)   ! APH lambda quantum number
   INTEGER(KIND=IW_Kind) tee(nstate)      ! APH tee quantum number
   NAMELIST/qlam/ lambda, tee
END MODULE Qlam_Module
