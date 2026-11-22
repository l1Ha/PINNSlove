MODULE Engpro_Module
   USE Numeric_Kinds_Module
   USE TotalEng_Module
   SAVE
   INTEGER(KIND=IW_Kind) numengs  ! Number of different energies at which scattering calculations will be preformed.
   INTEGER(KIND=IW_Kind) kenergy  ! Current energy index.
   REAL(Kind=WP_Kind) emin     ! Minimum scattering energy.
   REAL(Kind=WP_Kind) ejump    ! Incremental scattering energy.
   REAL(Kind=WP_Kind) engmax   ! Maximum scattering energy.
   NAMELIST/engpro/ emin, ejump, engmax, numengs, kenergy
END MODULE Engpro_Module
