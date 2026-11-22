MODULE Iopt_Module
   SAVE
   LOGICAL ibody   ! TRUE when using body-fixed coordinates.
   LOGICAL ispace  ! TRUE when using space-fixed coordinates.
   NAMELIST/iopt/ ibody, ispace
END MODULE Iopt_Module
