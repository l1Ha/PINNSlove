MODULE Elemnt_Module
   USE Narran_Module
   IMPLICIT NONE
   SAVE
   LOGICAL:: npar=.TRUE.     ! TRUE if parity decoupling is to be used.
   NAMELIST/elemnt/ npar
END MODULE Elemnt_Module
