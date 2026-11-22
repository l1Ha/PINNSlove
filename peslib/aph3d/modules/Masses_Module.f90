MODULE Masses_Module
   USE Numeric_Kinds_Module
   IMPLICIT NONE
   SAVE
   INTEGER NAtoms
   REAL(Kind=WP_Kind):: mass(3)      ! Array of Masses
   REAL(Kind=WP_Kind):: xmass(3)     ! Masses of the three atoms (A,B,C)
   REAL(Kind=WP_Kind):: Mtot         ! Total mass of the system
   REAL(Kind=WP_Kind):: USys         ! Hyperspherical reduced mass.
   REAL(Kind=WP_Kind):: Usys2        ! Two times the hyperspherical reduced mass.
   REAL(Kind=WP_Kind):: UABC         ! UABC=USys
   REAL(Kind=WP_Kind):: UABC2        ! UABC2=USys2
   REAL(Kind=WP_Kind):: UBC          ! Reduced mass of BC molecule
   REAL(Kind=WP_Kind):: UAB          ! Reduced mass of AB molecule
   REAL(Kind=WP_Kind):: UAC          ! Reduced mass of AC molecule
   REAL(Kind=WP_Kind):: U_Diatom(3)  ! Diatomic Reduced Masses
   CHARACTER(LEN=3) :: AtomicSymbol(3) = ['   ','   ','   ']  ! PORT: init blanks (AtomicMasses branch leaves them unset)
   CHARACTER(LEN=100) Notes(3)
   INTEGER AtomicNumber(3), MassNumber(3)
   REAL(KIND=WP_Kind) AtomicMass(3), Abundance(3), AtomicWeight(3)
   
   NAMELIST/Atoms/ AtomicSymbol, MassNumber
   NAMELIST/AtomicMasses/ mass
END MODULE Masses_Module
