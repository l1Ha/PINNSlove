MODULE EtaChn_Module
   USE Numeric_Kinds_Module
   USE Narran_Module
   SAVE
   REAL(Kind=WP_Kind) eta(narran)  ! Exponent coefficient for each arrangement channel. For coulombic systems only.
   NAMELIST/etachn/ eta
END MODULE EtaChn_Module
