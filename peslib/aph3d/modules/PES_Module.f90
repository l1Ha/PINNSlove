MODULE PES_MODULE
   IMPLICIT NONE
   SAVE
   CHARACTER(LEN=25) :: PES_Name='You must define PES_Name '   ! A unique name for the potential energy surface being used.
   CHARACTER(LEN=25) :: PesNaml ='                         '   ! A unique name for the potential energy surface being used.
   CHARACTER(LEN=25) :: PesNamu ='                         '   ! A unique name for the potential energy surface being used.
   CHARACTER(LEN=25) :: PesNamd ='                         '   ! A unique name for the potential energy surface being used.
   NAMELIST/PotentialNames/ PES_Name, PesNaml, PesNamu, PesNamd
END MODULE PES_MODULE
