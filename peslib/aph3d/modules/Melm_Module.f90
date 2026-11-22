MODULE Melm_Module
   SAVE
   CHARACTER(LEN=5) melmnt   !  melmnt='space'  For space-fixed calculations.
                              !  melmnt='probf'  For projectile body-fixed calculations.
                              !  melmnt='molbf'  For molecular body-fixed calculations.
   NAMELIST/melm/ melmnt
END MODULE Melm_Module
