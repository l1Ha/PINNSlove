MODULE Qall_Module
   USE parms_Module
   USE NState_Module
   USE Numeric_Kinds_Module
SAVE
   INTEGER(KIND=IW_Kind) mvib3(nstate)  ! vibrational quantum numbers
   INTEGER(KIND=IW_Kind) jrot3(nstate)  ! rotor quantum numbers
   INTEGER(KIND=IW_Kind) lorb3(nstate)  ! orbital quantum numbers
   INTEGER(KIND=IW_Kind) mega3(nstate)  ! projection quantum numbers
   INTEGER(KIND=IW_Kind) kchan(nstate)  ! arrangement channel labels, kchan(i)=1,2,3, i=1,nhc
   INTEGER(KIND=IW_Kind) nuj3(nstate)   ! state indices
   INTEGER(KIND=IW_Kind) elect3(nstate) ! electronic quantum numbers
   REAL(KIND=WP_Kind)    xksq3(nstate)  ! two times the reduced mass times the quantity of the total energy minus the internal energy
   REAL(KIND=WP_Kind)    Energy3_Delves(nstate) !
   REAL(KIND=WP_Kind)    xksq3_Delves(nstate)
   REAL(KIND=WP_Kind)    Energy3_Jacobi(nstate) !
   REAL(KIND=WP_Kind)    xksq3_Jacobi(nstate)
   REAL(KIND=WP_Kind)    Energy3_APH(nstate)
   REAL(KIND=WP_Kind)    xksq3_APH(nstate)
END MODULE Qall_Module
