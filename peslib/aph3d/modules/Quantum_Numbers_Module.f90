MODULE Quantum_Numbers_Module
   USE Numeric_Kinds_Module
   USE maxnstate_Module
   IMPLICIT NONE
   SAVE
   LOGICAL symmetry                                !
   INTEGER(KIND=IW_Kind) j_min                        !
   INTEGER(KIND=IW_Kind) j_max                        !
   INTEGER(KIND=IW_Kind) nvib(maxnstate)              ! vibrational quantum numbers
   INTEGER(KIND=IW_Kind) jrot(maxnstate)              ! rotor quantum numbers
   INTEGER(KIND=IW_Kind) lorb(maxnstate)              ! array of orbital quantum numbers
   INTEGER(KIND=IW_Kind) mega(maxnstate)              ! projection quantum numbers
   INTEGER(KIND=IW_Kind) chanl(maxnstate)             ! INTEGER array of arrangement channel labels, chanl(i)=1,2,3, i=1,nhc
   INTEGER(KIND=IW_Kind) nuj(maxnstate)               ! quantum indices
   INTEGER(KIND=IW_Kind) elect(maxnstate)             ! Electronic quantum numbers
   INTEGER(KIND=IW_Kind) tee(maxnstate)               !
   INTEGER(KIND=IW_Kind) lambda(maxnstate)            !
   REAL(Kind=WP_Kind) xksq(maxnstate)              ! two times the reduced mass times the quantity of the total energy minus the int
   REAL(Kind=WP_Kind) energies(maxnstate)          !
   REAL(Kind=WP_Kind) wvec(maxnstate)              !
   REAL(Kind=WP_Kind), ALLOCATABLE:: UnitMat(:,:)  !
   REAL(Kind=WP_Kind), ALLOCATABLE:: LorbMat(:,:)  !
   REAL(Kind=WP_Kind), ALLOCATABLE:: JrotMat(:,:)  !
END MODULE Quantum_Numbers_Module
