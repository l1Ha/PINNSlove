MODULE Approx_Module
   USE Numeric_Kinds_Module
   SAVE
   LOGICAL:: CenSud=.False.      ! censud  Set to true if centrifugal sudden approximation is desired.
   LOGICAL:: EngSud=.False.      ! engsud  Set to true if energy sudden approximation is desired.
   LOGICAL:: IosSud=.False.      ! iossud  Set to true if infinite order sudden approximation is desired.
   LOGICAL:: Pseudo=.False.      ! pseudo  If set to true in conjunction with one of the above approximations
                                 !         then the full coupled-channel equations is solved with the
                                 !         appropriate coupling set to zero. This option is for debugging purposes only.
   INTEGER(KIND=IW_Kind):: Lbar=0.  ! lbar    Average orbital quantum number used in the centrifugal sudden approximation.
   INTEGER(KIND=IW_Kind):: Jbar=0   ! jbar    Average rotational quantum number used in the energy sudden approximation.
   INTEGER(KIND=IW_Kind):: Mbar=0   ! mbar    Average projection quantum number used in the centrifugal or energy sudden approximations
   INTEGER(KIND=IW_Kind):: ICase=1  ! icase=1 Exact diagonal jtot(jtot+1)+j(j+1)-2omega*omega is used.
                                 ! icase=2 Diagonal elements are lorb*(lorb+1).
                                 ! icase=3 Diagonal elements are lbar*(lbar+1).
   NAMELIST/approx/ censud, engsud, iossud, pseudo, icase, lbar,jbar, mbar
END MODULE Approx_Module
