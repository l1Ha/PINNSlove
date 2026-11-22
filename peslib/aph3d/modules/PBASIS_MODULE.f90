MODULE PBASIS_MODULE
   IMPLICIT NONE
   SAVE
   INTEGER nlam(0:1000),lammin,lammax,megamin,maxmega,coupmega

   NAMELIST/basis/nlam,lammin,lammax,megamin,maxmega,coupmega
END MODULE PBASIS_MODULE
