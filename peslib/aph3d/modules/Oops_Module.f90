MODULE Oops_Module
   USE Numeric_Kinds_Module
   SAVE
   LOGICAL ioops            ! Set to true when angular vibrational functions are to be used.
   INTEGER jrott            !
   REAL(Kind=WP_Kind) rhoj  ! Radial coordiate where the angular vibrational functions are defined.
   REAL(Kind=WP_Kind) xktot !
   REAL(Kind=WP_Kind) xkold
   NAMELIST/oops/ rhoj, xktot, jrott, ioops
END MODULE Oops_Module
