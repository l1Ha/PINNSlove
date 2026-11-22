MODULE Dip_Module
USE Numeric_Kinds_Module
SAVE
   LOGICAL lphoto
   INTEGER(KIND=IW_Kind) jf
   INTEGER(KIND=IW_Kind) mf
   INTEGER(KIND=IW_Kind) lamf
   INTEGER(KIND=IW_Kind) prf
   INTEGER(KIND=IW_Kind) jd
   INTEGER(KIND=IW_Kind) md
   INTEGER(KIND=IW_Kind) lamd
   INTEGER(KIND=IW_Kind) prd
   INTEGER(KIND=IW_Kind) ji
   INTEGER(KIND=IW_Kind) mi
   INTEGER(KIND=IW_Kind) lami
   INTEGER(KIND=IW_Kind) pri
   INTEGER(KIND=IW_Kind) nsfuncu
   INTEGER(KIND=IW_Kind) nsfuncl
   INTEGER(KIND=IW_Kind) iaxis
   REAL(Kind=WP_Kind) enddip
   NAMELIST/dip/lphoto,jf,mf,lamf,prf,jd,md,lamd,prd, ji,mi,lami,pri,nsfuncu,nsfuncl,iaxis,enddip
END MODULE Dip_Module
