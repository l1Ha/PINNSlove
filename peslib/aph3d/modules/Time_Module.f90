MODULE Time_Module
   USE Numeric_Kinds_Module
   IMPLICIT NONE
   SAVE  
   REAL(Kind=WP_Kind), PARAMETER :: timeframe = 1000000
   REAL(Kind=WP_Kind) SecndsLast, SecndsNow, SecndsElapsed, t2, t1
END MODULE Time_Module
