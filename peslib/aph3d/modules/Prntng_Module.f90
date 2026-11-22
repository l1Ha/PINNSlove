MODULE Prntng_Module
   USE Numeric_Kinds_Module
   IMPLICIT NONE
   SAVE
   LOGICAL Little, Medium, Full
   INTEGER(KIND=IW_Kind) IthCll, IthSub
   INTEGER(KIND=IW_Kind) iprt(200)   ! level of debug WRITEs: 0=none, 1=minimum, 2=moderate, 3=maximum, 4=combination (must read in ite
   INTEGER(KIND=IW_Kind) iter(3,200) ! how many times to loop thru SUBROUTINE i with the jth debug level turned on (j=1 to 3)
   INTEGER(KIND=IW_Kind) nsubs       ! number of SUBROUTINEs which will have their debug WRITEs turned on
   INTEGER(KIND=IW_Kind) iprtall     !
   !CHARACTER(LEN=8) ::   option='all     '
   CHARACTER(LEN=8) ::   option='selected'
END MODULE Prntng_Module
