MODULE Qnumsp_Module
   USE Numeric_Kinds_Module
   USE Narran_Module
   SAVE
   INTEGER(KIND=IW_Kind), ALLOCATABLE:: numax(:)  ! Maximum vibrational quantum number for each rotational state.
   INTEGER(KIND=IW_Kind), ALLOCATABLE:: numin(:)  ! Minimum vibrational quantum number for each rotational state.
!  NAMELIST/qnumsp/ numax, numin
END MODULE Qnumsp_Module
