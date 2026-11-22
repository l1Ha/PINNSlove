MODULE DVR_Module
   USE Numeric_Kinds_Module
   IMPLICIT NONE
   SAVE
   INTEGER(KIND=IW_Kind) ndvptmax ! maximum number of DVR points
   INTEGER(KIND=IW_Kind), PARAMETER:: nthetmax=100              ! maximum number of theta points
   INTEGER(KIND=IW_Kind), PARAMETER:: nchimax=200               ! maximum number of chi points
   INTEGER(KIND=IW_Kind), PARAMETER:: nmax=nchimax              ! nchimax=maximum number of dvr points kept at each theta
   INTEGER(KIND=IW_Kind), PARAMETER:: nt2max=8400               ! maximum size of final diagonalization
   INTEGER(KIND=IW_Kind), PARAMETER:: nsfnmax=8400              ! max number of surface functions that can be saved.
END MODULE DVR_Module
