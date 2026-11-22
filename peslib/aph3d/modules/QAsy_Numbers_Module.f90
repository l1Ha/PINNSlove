MODULE QAsy_Numbers_Module
   USE Numeric_Kinds_Module
!   USE Masses_Module!, ONLY::AtomicSymbol, Notes, AtomicNumber, MassNumber, Abundance, AtomicWeight, AtomicMass, Usys2, NAtoms
   IMPLICIT NONE
   SAVE

   !INTEGER(Kind=IW_Kind) NAtoms
   !CHARACTER(LEN=3),      ALLOCATABLE:: AtomicSymbol(:)  ! Atomic symbol
   !CHARACTER(LEN=100),    ALLOCATABLE:: Notes(:)         ! Notes about the atoms
   !INTEGER(Kind=IW_Kind), ALLOCATABLE:: AtomicNumber(:)  ! Atomic number
   !INTEGER(Kind=IW_Kind), ALLOCATABLE:: MassNumber(:)    ! Mass number
   !REAL(Kind=WP_Kind),    ALLOCATABLE:: Abundance(:)     ! Isotopic abundance
   !REAL(Kind=WP_Kind),    ALLOCATABLE:: AtomicWeight(:)  ! Atomic Weight
   !REAL(Kind=WP_Kind),    ALLOCATABLE:: AtomicMass(:)    ! Atomic mass
   !REAL(Kind=WP_Kind) :: USys2                           ! two times the reduced mass
   
   
   CHARACTER(LEN=5),      ALLOCATABLE:: QLabels(:)       ! Quantum labesl
   CHARACTER(LEN=9),      ALLOCATABLE:: ConLabels(:)     ! Conserved labels
   
   INTEGER(Kind=IW_Kind), PARAMETER:: NQuant=5
   INTEGER(Kind=IW_Kind), PARAMETER:: NConserved=4
   INTEGER(Kind=IW_Kind), ALLOCATABLE:: QNumbr(:,:)      ! chanl, elec, nvib, jrot, lorb
   INTEGER(Kind=IW_Kind), ALLOCATABLE:: Conserved(:)     ! Conserved quantum numbers
   REAL(Kind=WP_Kind),    ALLOCATABLE:: xksq(:)          ! two times the reduced mass times the quantity of the total energy minus the int
   REAL(Kind=WP_Kind),    ALLOCATABLE:: AsymptEnergies(:)! Asymptotic internal energies
   REAL(Kind=WP_Kind),    ALLOCATABLE:: wvec(:)          ! wavevector
   REAL(Kind=WP_Kind),    ALLOCATABLE:: UnitMat(:,:)     ! unit matrix
   REAL(Kind=WP_Kind),    ALLOCATABLE:: LorbMat(:,:)     ! lorb matrix
   REAL(Kind=WP_Kind),    ALLOCATABLE:: JrotMat(:,:)     ! jrot matrix
END MODULE QAsy_Numbers_Module
