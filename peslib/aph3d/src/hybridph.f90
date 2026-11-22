SUBROUTINE hybridph (nt, rmat, eye11, w, vecnew, tstore, vecold, xsq, rbegin, rmax)
!

!             P U R P O S E    O F    S U B R O U T I N E

! Integrates second-order coupled differential equations across a
! region.
! written by G. A. Parker

!             I N P U T   A R G U M E N T S

!   nt        Number of coupled second-order differential equations.
!   rmat      Wigner R-matrix at the begining of the region.
!   eye11     Temporary storage array.
!   w         Temporary storage array.
!   vecnew    Temporary storage array.
!   tstore    Temporary storage array.
!   vecold    Temporary storage array.
!   xsq       Temporary storage array.
!   rbegin    Propagation distance at the beginning of the region.
!   rmax      Propagation distance at the end of the region.
!   nsteps     The number of steps in the log-derivative region.
!   rbegin    Beginning distance for the log-derivative region.
!   SwitchToVIVS    Ending distance for the log-derivative region and
!             the starting distance for the vivs region.
!   rmax      Ending distance for the vivs region.
!             if SwitchToVIVS<rmax then the integration goes from
!             rbegin to SwitchToVIVS and rmax is set equal to SwitchToVIVS.

!             O U T P U T     A R G U M E N T S

!   rmat    Wigner R-matrix at the end of the region.

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  USE FileUnits_Module
  USE nstate_Module
  USE parms_Module
  USE Arrch_Module
  USE melm_Module
  USE Propag_Module
  USE rhovalue_Module
  USE region_Module
  USE Oops_Module
  USE VFunc_Module
  USE tlogd_Module
  USE Storage_Module
  USE dip_Module
  USE Boundary_Module
  USE BasisRegions_Module
  IMPLICIT NONE
  SAVE
  LOGICAL isym, ifalse, little, medium, full
  ! PORT (2026-10-02): state for the Delves-region sector grid
  LOGICAL delv_grid_done
  REAL(Kind=WP_Kind) delv_rbegin, delv_rmax
  DATA delv_grid_done/.false./, delv_rbegin/0.d0/, delv_rmax/0.d0/
  INTEGER nt, n, n1, n2, isect, nrho,i, j, nconsidr, iset, ithsub, ithcll
  REAL(Kind=WP_Kind) rmat, eye11, w, vecnew, deltamin, deltamax
  REAL(Kind=WP_Kind) tstore, vecold, xsq, rbegin, rmax
  REAL(Kind=WP_Kind) rright, rleft, rhov1, olast, rsave, xk
  REAL(Kind=WP_Kind) rhov2, ovrtol, rhovals, okeep
!-----------------------------------------------------------------------
! the following should have length of at least nchan where nchan
! is the number of channels to be integrated.
!-----------------------------------------------------------------------
  DIMENSION xsq(nstate), xk(nstate)
!-----------------------------------------------------------------------
! arrays dimensioned as matrices (nchan*nchan).
!-----------------------------------------------------------------------
  DIMENSION eye11(*), w(*), rmat(*), vecnew(nt,nt), tstore(*), vecold(*)   ! PORT: assumed-size
  ! PORT (2026-10-02): optional dump of the log-derivative (R) matrix in the
  ! APH sector basis, used to validate the 2x2 inner propagator against the
  ! single-surface chain.  Enabled by APH3D_RDUMP=1 (window APH3D_RDUMP_LO/HI).
  LOGICAL rdump_on, rdump_first
  REAL(Kind=WP_Kind) rdump_lo, rdump_hi
  CHARACTER(LEN=64) rdump_env
  INTEGER rdump_len, rdump_stat
  DATA rdump_first/.true./
  EXTERNAL popt, upsiln, sector, matgenph, mxout, overlap, trnsfm

!-----------------------------------------------------------------------
! set default data
!-----------------------------------------------------------------------

  DATA isym /.false./, ifalse /.false./

!-----------------------------------------------------------------------
! Obtain print options
!-----------------------------------------------------------------------
  DATA ithcll/0/,ithsub/0/, isect/0/
  DATA little/.false./,medium/.false./,full/.false./
  CALL popt('hybrid  ', little, medium, full, ithcll, ithsub)

  IF(rdump_first)THEN
     rdump_first = .false.
     rdump_on = .false.; rdump_lo = 25.d0; rdump_hi = 30.d0
     CALL get_environment_variable('APH3D_RDUMP', rdump_env, rdump_len, rdump_stat)
     IF(rdump_stat==0 .AND. rdump_len>0)THEN
        IF(rdump_env(1:1)=='1' .OR. rdump_env(1:1)=='T' .OR. rdump_env(1:1)=='t') rdump_on=.true.
     ENDIF
     CALL get_environment_variable('APH3D_RDUMP_LO', rdump_env, rdump_len, rdump_stat)
     IF(rdump_stat==0 .AND. rdump_len>0) READ(rdump_env(1:rdump_len),*) rdump_lo
     CALL get_environment_variable('APH3D_RDUMP_HI', rdump_env, rdump_len, rdump_stat)
     IF(rdump_stat==0 .AND. rdump_len>0) READ(rdump_env(1:rdump_len),*) rdump_hi
  ENDIF
  okeep = 1.0d0 - ovrtol
!-----------------------------------------------------------------------
! Write or read data for subsequent energies.
!-----------------------------------------------------------------------

  nrho = 1
!-----------------------------------------------------------------------
! Write out intgration parameters.
!-----------------------------------------------------------------------

  WRITE(Out_Unit,40) nsteps, rbegin, SwitchToVIVS, rmax

!-----------------------------------------------------------------------
! calculate lower and upper bounds for energy-labeled
! channels which are to be propagated
!-----------------------------------------------------------------------

  n = nt
  n1 = 1
  n2 = nt
  WRITE(Out_Unit,110) nt, n, n1, n2
  olast = 0.0d0
  rleft = rbegin
  rright = rmax
  rhov1 = rbegin
  !isect = 0
  olast = 0.0d0

!-----------------------------------------------------------------------
! Calculate the vibration-rotational basis.
!-----------------------------------------------------------------------

  IF(iregion /= 'aph    ')THEN
     CALL upsiln (vecold, tstore, xsq, eye11, rhov1, chinuj)
  ENDIF

!-----------------------------------------------------------------------
! if one is outside of the jacobi region then determine the
! the sector length.
!-----------------------------------------------------------------------
10 WRITE(Out_Unit,*)'Starting hybridph'
  IF(iregion /= 'jacobi ')THEN
     !CALL sector1 (rhov1, rhov2, olast, okeep, rleft, rright, rmax, isect, nrho)
     !CALL sector (rhov1, rhov2, rleft, rright, rmax, isect, nrho)
     ! PORT (2026-10-02): the Delves region has no Ovr_Unit sector records and
     ! the APH-region arrays Basis_Dist/Sector_Boundaries do not cover it, so
     ! build a Delves sector grid for [rbegin, rmax] and step through it with
     ! an internal counter (previously this path read past the arrays and
     ! crashed at line 129).
     IF(iregion == 'delves ')THEN
        IF(.NOT.delv_grid_done)THEN
           delv_rbegin = MIN(rbegin, rmax)
           delv_rmax   = MAX(rbegin, rmax)
           CALL Delves_Sector_Grid(delv_rbegin, delv_rmax)
           delv_grid_done = .TRUE.
           isect = 0
           WRITE(Msg_Unit,'("Delves sector grid: ",I6," sectors, ",ES12.4," -> ",ES12.4)') &
                 NSectors, delv_rbegin, delv_rmax
           WRITE(Out_Unit,'("Delves sector grid: ",I6," sectors, ",ES12.4," -> ",ES12.4)') &
                 NSectors, delv_rbegin, delv_rmax
        ENDIF
        isect = isect + 1
        IF(isect > NSectors) GOTO 30
     ELSE
        READ(Ovr_Unit)ISect,Rhov1,rhov2
        IF(ISect > NSectors) GOTO 30
     ENDIF
     Rhov1=Basis_Dist(Isect)
     IF(ISect/=NSectors)Then
        Rhov2=Basis_Dist(Isect+1)
     ELSE
        Rhov2=Basis_Dist(Isect)
     ENDIF
     Rleft=Sector_Boundaries(Isect)
     Rright=Sector_Boundaries(Isect+1)
     Rmax=Sector_Boundaries(Nsectors+1)
     Nrho=Nsectors
     ! PORT (2026-10-02): the energy-independent storage (frst_unit) is
     ! consumed by the per-energy path and must contain the APH-region records
     ! only; do not append Delves-region records.
     IF(iregion /= 'delves ') WRITE(frst_unit) rhov1, rhov2, rleft, rright, rmax, isect, nrho
     IF(Isect==1)THEN
       WRITE(Msg_Unit,*)"Nrho=",nrho
       WRITE(Msg_Unit,'(A5,5A14)')"ISect", "rleft", "rhov1", "rright", "rhov2", "rmax"
     ENDIF
     WRITE(Msg_Unit,'(I5,5ES14.6)')isect, rleft, rhov1, rright, rhov2, rmax
     IF(lphoto)THEN
        WRITE(fph_unit) rhov1, rhov2, rleft, rright, rmax, isect, nrho
     ENDIF
  ENDIF
!-----------------------------------------------------------------------
! integrate from rleft to rright using the log-derivative integrator.
!-----------------------------------------------------------------------

  rsave = rright
  IF(rleft < SwitchToVIVS)THEN
     IF(rright > SwitchToVIVS) rright = SwitchToVIVS
     IF(little) WRITE(Out_Unit,100) rleft,rright,nsteps
     CALL matgenph (tstore, rmat, w, nt, rleft, rright, nsteps, xsq)
     IF(full)THEN
        WRITE(Out_Unit,90) rright
        CALL MxOut(rmat, nt, nt)
     ENDIF
     IF(rright < SwitchToVIVS) GOTO 20
  ENDIF

  IF(rleft >= rmax) GOTO 20
!-----------------------------------------------------------------------
! If we are not in the Jacobi region, sector adiabatic basis functions
! are used.  Therefore, at the end of each sector we need to obtain
! the overlap matrix and transform to the new basis.
!-----------------------------------------------------------------------

20 IF(iregion /= 'jacobi ')THEN
     IF(isect == nrho-1) GOTO 30
        WRITE(Out_Unit,*)
        WRITE(Out_Unit,*)'Calling overlap'
        CALL overlap (rhov1, rhov2, vecnew, nt, vecold, tstore, xsq, xk, eye11, olast)
        WRITE(Out_Unit,*)'Now back in hybridph'
        WRITE(Out_Unit,*)'rhov1=',rhov1,rhov2,nt
        IF(iregion /= 'delves ')THEN
           WRITE(frst_unit)rhov1, rhov2, nt
           DO j = 1,nt
              WRITE(frst_unit)(vecnew(i,j),i = 1,nt)
           ENDDO
        ENDIF
!cccc            WRITE(frst_unit)(xsq(i),i = 1,nt), (xk(i),i = 1,nt)
     IF(rdump_on .AND. rhov2 >= rdump_lo .AND. rhov2 <= rdump_hi)THEN
        OPEN(Unit=93, File='rdump.txt', Status='UNKNOWN', Position='APPEND')
        WRITE(93,'(A,ES22.14,A,I6)') '# R matrix at rho=', rhov2, '  nt=', nt
        DO i = 1, nt
           WRITE(93,'(200ES22.14)') (rmat((i-1)*nt+j), j = 1, nt)
        ENDDO
        CLOSE(93)
     ENDIF
     rhov1 = rhov2

!-----------------------------------------------------------------------
! End of sector. Go to the beginning and propagate across another
! sector.
!-----------------------------------------------------------------------

     GOTO 10
  ENDIF
30 CONTINUE
  WRITE(Out_Unit,*)'Leaving hybridph'
  RETURN

40 FORMAT(//,'nsteps=',i5,5x,'rstart=',f10.5,5x,'SwitchToVIVS=',f10.5,5x, 'rmax=',f10.5)
90 FORMAT(//,'r-matrix at r=',f10.5)
100 FORMAT(//,'rstart=',f10.5,5x,'rend=',f10.5,5x,'nsteps=',i5,5x)
110 FORMAT(//,5x,'nt=',i4,5x,'n=',i4,5x,'n1=',i4,5x,'n2=',i4,5x,'eshift=',e14.7)
  END
