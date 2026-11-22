MODULE Energy_Module
   USE Numeric_Kinds_Module
   USE Convrsns_Module
   IMPLICIT NONE
   SAVE
   INTEGER            :: NEnergy     = 101          ! Number of energies
   LOGICAL            :: eV_Input    = .False.      ! If True then input energies are in eV
   REAL(Kind=WP_Kind) :: Efirst_eV   = .3D0         ! First energy=.3eV
   REAL(Kind=WP_Kind) :: DeltaEng_eV = .01D0        ! Energy grid spacing=.01eV
   REAL(Kind=WP_Kind) :: efirst      = .3/autoev    ! First energy=.3eV converted to Ha
   REAL(Kind=WP_Kind) :: deltaeng    = .01/autoev   ! Energy grid spacing=.01 converted to Ha
   ! PORT (2026-10-03): reference energy (Ha) of the asymptotic channel
   ! thresholds relative to the zero used in 'efirst'.  The open/closed
   ! decision and the asymptotic k^2 use (efirst - (AsymptEnergies - eref)),
   ! i.e. efirst is the collision energy measured from the channel at eref.
   ! eref = 0 keeps the historical behaviour (AsymptEnergies already relative).
   REAL(Kind=WP_Kind) :: eref        = 0.D0
   LOGICAL Asymptotic_Test
   NAMELIST/TotEnergy/ efirst, deltaeng, nenergy, eV_Input, Efirst_eV, DeltaEng_eV, eref
END MODULE Energy_Module
