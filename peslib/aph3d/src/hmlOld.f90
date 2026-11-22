SUBROUTINE hmlOld (nrho, rho_val, ntheta, tht_val, nchi, nchir, &
 chi_val, t_rho, t_theta, t_chi, t_chir, t_chi1, t_chir1, hamil, &
 v_pot, nirrep, irrep, n1, nsize, jtot, parity, usys2, a, b, c, d, &
 rhom2, tblock)
USE FileUnits_Module
IMPLICIT NONE
LOGICAL :: debug
INTEGER :: nirrep
INTEGER :: nrho, ntheta, nchi, nchir, n1 (nirrep), irho, itheta, &
 ichi, irrep, i, j, nsize, Indxr_New, jrho, jtheta, jchi, nbegin, &
 nlast, jrrep, jtot, lam, lammin, lammax, parity

REAL(KIND=WP_KinD) :: t_rho (nrho, nrho), t_theta (ntheta, ntheta), t_chi ( &
 nchi, nchi), t_chir (nchir, nchir), v_pot (nrho, ntheta, nchi), &
 rho_val (nrho), tht_val (ntheta), chi_val (nchi + 1), t_chi1 ( &
 nchi, nchi), t_chir1 (nchir, nchir), hamil (nsize, nsize), &
 usys2, a (ntheta), b (ntheta), c (ntheta), d (ntheta), rhom2 ( &
 nrho), afac, afactor, cfac, cfactor, v_eff, v0, v1, v2, v3, &
 tblock (nchi, nchir)
WRITE(Out_Unit, * )
WRITE(Out_Unit, * ) 'Begin(hmlold)'
debug = .False.
!--------------------------------------------------------------------
! Determine Minimum and Maximum lambda projection.
!--------------------------------------------------------------------
IF (parity==0)THEN
   lammin = 0
ELSE
   lammin = 1
ENDIF

lammax = jtot
!--------------------------------------------------------------------
! WRITE input PARAMETERs
!--------------------------------------------------------------------
IF (debug)THEN
   WRITE(Out_Unit,*) 'nrho=  ', nrho
   WRITE(Out_Unit,*) 'ntheta=', ntheta
   WRITE(Out_Unit,*) 'nchi=  ', nchi,   '   nchir= ', nchir
   WRITE(Out_Unit,*) 'nsize= ', nsize
   WRITE(Out_Unit,*) 'nirrep=', nirrep, '   irrep= ', irrep
   WRITE(Out_Unit,*) 'n1=    ', n1
   WRITE(Out_Unit,*) 'jtot=  ', jtot
   WRITE(Out_Unit,*) 'parity=', parity
   WRITE(Out_Unit,*) 'lammin=', lammin, '   lammax=', lammax
ENDIF

!--------------------------------------------------------------------
! Store
!--------------------------------------------------------------------
nbegin = 1
DO jrrep = 1, irrep - 1
   nbegin = nbegin + n1 (jrrep)
ENDDO
nlast = nbegin + n1 (irrep) - 1

IF (debug)THEN
   WRITE(Out_Unit,  * ) 'nbegin', nbegin, '   nlast=', nlast
ENDIF

DO ichi = nbegin, nlast
   DO jchi = nbegin, nlast
      t_chir (ichi + 1 - nbegin, jchi + 1 - nbegin) = t_chi (ichi, jchi)
      t_chir1 (ichi + 1 - nbegin, jchi + 1 - nbegin) = t_chi (ichi, jchi)
   ENDDO
ENDDO

IF (debug)THEN
   WRITE(Out_Unit, * ) 't_chir'
   CALL MxOut(t_chir, nchir, nchir)
   WRITE(Out_Unit, * ) 't_chir1'
   CALL MxOut(t_chir1, nchir, nchir)
ENDIF
!--------------------------------------------------------------------
! Initialize the Hamiltonian to zero
!--------------------------------------------------------------------

hamil = 0.d0
!--------------------------------------------------------------------
! Calculate diagonal Lambda blocks.
!--------------------------------------------------------------------
DO lam = lammin, lammax
   IF (debug)THEN
      WRITE(Out_Unit, * ) 'lam=', lam
   ENDIF
!--------------------------------------------------------------------
! Rho Kinetic energy term.
!--------------------------------------------------------------------
   DO irho = 1, nrho
      DO jrho = 1, nrho
         DO itheta = 1, ntheta
            DO ichi = 1, nchir
               i = Indxr_New (lam, lammin, nrho, ntheta, nchir, irho, itheta, ichi)
               j = Indxr_New (lam, lammin, nrho, ntheta, nchir, jrho, itheta, ichi)
               hamil (i, j) = t_rho (irho, jrho)
            ENDDO
         ENDDO
      ENDDO
   ENDDO

   IF (debug)THEN
      WRITE(Out_Unit, * ) 'After t_rho'
      CALL MxOut(hamil, nsize, nsize)
   ENDIF

!--------------------------------------------------------------------
! Chi Kinetic energy term.
!--------------------------------------------------------------------
   IF(nchir>1)THEN
   DO irho = 1, nrho
      DO itheta = 1, ntheta
         DO ichi = 1, nchir
            DO jchi = 1, nchir
               i = Indxr_New (lam, lammin, nrho, ntheta, nchir, irho, itheta, ichi)
               j = Indxr_New (lam, lammin, nrho, ntheta, nchir, irho, itheta, jchi)
               hamil (i, j) = t_chir (ichi, jchi) * 2.d0 * b (itheta) * rhom2 ( irho) + hamil (i, j)
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   ENDIF

   IF (debug)THEN
      WRITE(Out_Unit, * ) 'After t_chir'
      CALL MxOut(hamil, nsize, nsize)
   ENDIF

!---------------------------------------------------------------------
! Theta Kinetic energy term.
!---------------------------------------------------------------------
   IF(ntheta>1)THEN
   DO irho = 1, nrho
      DO itheta = 1, ntheta
         DO jtheta = 1, ntheta
            DO ichi = 1, nchir
               i = Indxr_New (lam, lammin, nrho, ntheta, nchir, irho, itheta, ichi)
               j = Indxr_New (lam, lammin, nrho, ntheta, nchir, irho, jtheta, ichi)
               !hamil (i, j) = hamil (i, j) + t_theta (itheta, jtheta) * rhom2 (irho)
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   ENDIF

IF (debug)THEN
   WRITE(Out_Unit, * ) 'After t_theta'
   CALL MxOut(hamil, nsize, nsize)
ENDIF

!---------------------------------------------------------------------
! Add Effective potential
!---------------------------------------------------------------------
   IF (debug)THEN
      WRITE(Out_Unit, * ) 'jtot=', jtot
      WRITE(Out_Unit, * ) 'lam=', lam
   ENDIF
   WRITE(Out_Unit,  * ) 'irrep=', irrep, '   nchir=', nchir, '   nbegin=', nbegin
   DO irho = 1, nrho
      DO itheta = 1, ntheta
         DO ichi = 1, nchir
            i = Indxr_New (lam, lammin, nrho, ntheta, nchir, irho, itheta, ichi)
            v0 = 0.d0
            v0 = v_pot (irho, itheta, ichi)
            hamil (i, i) = hamil (i, i) + v0
            v1 = jtot * (jtot + 1) * (a (itheta) + b (itheta) ) / 2.d0 * rhom2 (irho) / usys2
            v1 = 0.d0
            hamil (i, i) = hamil (i, i) + v1
            v2 = 15.d0 / usys2 / 4.d0 * rhom2 (irho)
            v2 = 0.d0
            hamil (i, i) = hamil (i, i) + v2
            v3 = (2.d0 * c (itheta) - (a (itheta) + b (itheta) ) ) * lam**2 * rhom2 (irho) / usys2
            v3 = 0.d0
            hamil (i, i) = hamil (i, i) + v3
            v_eff = v0 + v1 + v2 + v3
         ENDDO
      ENDDO
   ENDDO

   IF (debug)THEN
      WRITE(Out_Unit, * ) 'After V_eff'
      CALL MxOut(hamil, nsize, nsize)
   ENDIF

ENDDO

!-------------------------------------------------------------------
! construct Coriolis coupling terms
!-------------------------------------------------------------------
DO lam = lammin, lammax - 1
   cfactor = cfac (jtot, lam, lam + 1, parity)
   DO irho = 1, nrho
      DO itheta = 1, ntheta
         DO ichi = 1, nchir
            DO jchi = 1, nchir
               i = Indxr_New (lam, lammin, nrho, ntheta, nchir, irho, itheta, ichi)
               j = Indxr_New (lam + 1, lammin, nrho, ntheta, nchir, irho, itheta, jchi)
               hamil (i, j) = hamil (i, j) - rhom2 (irho) / usys2 * d (itheta) * t_chir1 (ichi, jchi) * cfactor
               hamil (j, i) = hamil (i, j)
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

IF (debug)THEN
   WRITE(Out_Unit, * ) 'After Coriolis'
   CALL MxOut(hamil, nsize, nsize)
ENDIF

!-------------------------------------------------------------------
! construct Asymmetric Top coupling term for lambda=lambda'=1
! This is a diagonal contribution for only lambda=1.
!-------------------------------------------------------------------
IF (lammax>=1)THEN
   cfactor = cfac (jtot, 1, 1, parity)
   DO irho = 1, nrho
      DO itheta = 1, ntheta
         DO ichi = 1, nchir
            i = Indxr_New (1, lammin, nrho, ntheta, nchir, irho, itheta, ichi)
            hamil (i, i) = hamil (i, i) - rhom2 (irho) / usys2 * (a (itheta) - b (itheta) ) / 2.d0 * cfactor
         ENDDO
      ENDDO
   ENDDO
ENDIF

IF (debug)THEN
   WRITE(Out_Unit, * ) 'After Asymmetric Top'
   CALL MxOut(hamil, nsize, nsize)
ENDIF

!-------------------------------------------------------------------
! Construct Asymmetric Top Coupling terms.
!-------------------------------------------------------------------
DO lam = lammin, lammax - 2
   afactor = afac (jtot, lam, lam + 2, parity)
   DO irho = 1, nrho
      DO itheta = 1, ntheta
         DO ichi = 1, nchir
            i = Indxr_New (lam, lammin, nrho, ntheta, nchir, irho, itheta, ichi)
            j = Indxr_New (lam + 2, lammin, nrho, ntheta, nchir, irho, itheta, ichi)
            hamil (i, j) = hamil (i, j) + (a (itheta) - b (itheta) ) / 2.d0 * rhom2 (irho) * afactor / usys2
            hamil (j, i) = hamil (i, j)
         ENDDO
      ENDDO
   ENDDO
ENDDO

IF (debug)THEN
   WRITE(Out_Unit, * ) 'After Asymmetric Top (Lambda=1)'
   CALL MxOut(hamil, nsize, nsize)
ENDIF

CALL sym_check (hamil, nsize)
WRITE(Out_Unit, * ) 'END(hmlold)'

WRITE(Out_Unit, * )
RETURN
ENDSUBROUTINE hmlOld
