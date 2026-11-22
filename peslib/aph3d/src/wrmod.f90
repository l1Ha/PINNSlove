SUBROUTINE wrmod (energy, phi, idi, numnp, neq, nfreq, nq, xnode, beta, nid, nidm, ndisce, id, nlancz, nel)
USE Numeric_Kinds_Module
USE Narran_Module
USE Masses_Module
USE STST_Module
USE potzero_Module
USE rhosur_Module
USE rhowrt_Module
USE totj_Module
USE convrsns_Module
Use FileUnits_Module
USE Mshdat_Module,only:nave1_FEM,nave2_FEM
IMPLICIT NONE

!
!   $RCSfile: wrmod.f,v $ $Revision: 1.15 $
!   $Date: 89/11/16 11:03:36 $
!   $State: Stable $
!

!            P U R P O S E   O F    S U B R O U T I N E

!-----------------------------------------------------------------------
! This routine WRITEs the surface functions on disk.
!-----------------------------------------------------------------------

!            I N P U T    A R G U M E N T S
!    energy
!    phi
!    idi
!    numnp
!    neq
!    nfreq
!    nq
!    xnode
!    beta
!    nid
!    nidm      maximum number of terms in the constraint equation.
!    ndisce
!    id
!    nlancz
!    nel
LOGICAL little, medium, full
INTEGER ithcall, ithsub
INTEGER numnp, neq, nfreq, ndisce, nidm, nel
INTEGER nq, jtheqn, irho, kthnode, ithmode, nlancz, ndum, nqloc, i, j, nlambda, ione
INTEGER idi(nidm,ndisce), nid(ndisce), id(numnp)
REAL(Kind=WP_Kind) energy(nfreq), phi(numnp), xnode(numnp,nq), beta(nidm,ndisce)
REAL(Kind=WP_Kind) zero, two, eigenav
EXTERNAL popt, condis

PARAMETER (zero=0.d0, two=2.d0, ione=1)

!-----------------------------------------------------------------------
!  Determine printing options.
!-----------------------------------------------------------------------
DATA ithcall/0/, ithsub/0/
DATA little/.false./, medium/.false./, full/.false./

CALL popt ('wrmod   ', little, medium, full, ithcall, ithsub)
!-----------------------------------------------------------------------
! WRITE the surface functions and energies to disk for use by
! AphDel_Old or AphDel_New in the propagator.  Note that only the unconstrained part
! of the functions is written along with the information needed by
! condis to reconstruct the whole function.  This saves approx.
! a factor of two in disk space.
!-----------------------------------------------------------------------
irho=1
nlambda=1
IF(rhosurf>=rhowrit)THEN
  REWIND Temp2_FEM_Unit
  WRITE(Nmodes_FEM_Bin_Unit) jtot, mega, parity, symmetry, numnp, mass(1), mass(2), mass(3), &
    rhosurf, nel, nfreq, irho, nlambda, ndisce, nidm, neq, betaw, seq, chncof
  WRITE(Nmodes_FEM_Bin_Unit) (energy(ithmode),ithmode = 1,nfreq)
  WRITE(Nmodes_FEM_Bin_Unit) (nid(i), i=1,ndisce)
  WRITE(Nmodes_FEM_Bin_Unit) ((idi(i,j), i=1,nidm), j=1,ndisce)
  WRITE(Nmodes_FEM_Bin_Unit) ((beta(i,j), i=1,nidm), j=1,ndisce)
  WRITE(Nmodes_FEM_Bin_Unit) (id(i), i=1,numnp)
ENDIF
!-----------------------------------------------------------------------
!  Calculate the energy in Hartree atomic units.
!-----------------------------------------------------------------------
DO ithmode = 1,nfreq
   energy(ithmode) = energy(ithmode)/(two*rhosurf**2*usys)
   IF(rhosurf>=rhowrit)THEN
      READ(Temp2_FEM_Unit) (phi(jtheqn),jtheqn = 1,neq)
      WRITE(Nmodes_FEM_Bin_Unit) ithmode,(phi(jtheqn),jtheqn = 1,neq)
   ENDIF
ENDDO
!---------------------------------------------------------------------
!  Calculate the average of the first nave2_FEM-nave1_FEM+1 energies in eV.
!  This needs to be changed to READ nav from DATA.
!---------------------------------------------------------------------
eigenav=0.d0
DO ithmode=nave1_FEM,nave2_FEM
   eigenav=eigenav+energy(ithmode)
ENDDO
eigenav = autoev*eigenav/(nave2_FEM-nave1_FEM+1)
!-----------------------------------------------------------------------
! WRITE the calculated eigenvalues in electron volts.
!-----------------------------------------------------------------------
WRITE(Out_unit,*) ' The calculated eigenvalues in eV at rho=: ', rhosurf
DO Ithmode=1,nfreq
   WRITE(Out_unit,'(1x,4d17.9)') autoev*energy(ithmode) !(autoev*energy(ithmode),ithmode = 1,nfreq) !GregParker TEMPMOD
ENDDO
WRITE(Out_Unit,*)' rho, nave1_FEM, ', 'nave2, eigenav=',rhosurf, nave1_FEM, nave2_FEM, eigenav
WRITE(Msg_Unit,*)' The calculated eigenvalues in eV at rho=: ', rhosurf," eigenav=",eigenav
WRITE(Msg_Unit,*) (autoev*energy(ithmode),ithmode = 1,MIN(12,nfreq))

!-----------------------------------------------------------------------
! WRITE the calculated eigenvalues in electron volts to a seperate file.
!-----------------------------------------------------------------------
WRITE(Wrmod_Unit,5)
WRITE(Wrmod_Unit,6)ione
WRITE(Wrmod_Unit,3)rhosurf
WRITE(Wrmod_Unit,7)
WRITE(Wrmod_Unit,6)nfreq
WRITE(Wrmod_Unit,4)(autoev*(energy(ithmode)-vzero),ithmode=1,nfreq)
WRITE(Wrmod_Unit,8)

!-----------------------------------------------------------------------
! Calculate the surface function at all nodal points.
! WRITE the surface functions onto disk for subsequent rhos.
!-----------------------------------------------------------------------
!     IF(nlancz==0)THEN
  nqloc=nq
!    ELSE
!       nqloc=nfreq
!     ENDIF
!
!  -----------------------------------------------------------------
!  WRITE the last set of fcns calculated on Nzromg_Unit for use at next mega.
!  IF mega=megamin WRITE the functions on ZroMeg_Unit for use at next rho.
!  use Lancos80_Unit and Lancos81_Unit instead of ZroMeg_Unit and Nzromg_Unit IF lanczos is used.
!  -------------------------------------------------------------------
IF(mega==megamin)THEN
   REWIND ZroMeg_Unit
   WRITE(ZroMeg_Unit) nfreq,nqloc,numnp
   IF(nlancz==1)THEN
      REWIND Lancos80_Unit
      WRITE(Lancos80_Unit) nfreq,nqloc,numnp
   ENDIF
ENDIF

REWIND Nzromg_Unit
WRITE(Nzromg_Unit) nfreq,nqloc,numnp
IF(nlancz==1)THEN
   REWIND Lancos81_Unit
   WRITE(Lancos81_Unit) nfreq,nqloc,numnp
ENDIF
!
IF(nlancz==1)THEN
   !-----------------------------------------------------------------------
   ! Calculate the Lanczos eigenvectors at all nodal points.
   ! WRITE the Lanczos vectors onto disk for subsequent rhos.
   !-----------------------------------------------------------------------
   REWIND Lancos_Phi_Unit
   !     DO 95 ithmode = 1,nfreq
   DO ithmode = 1,nqloc
      READ(Lancos_Phi_Unit) ndum, (phi(jtheqn),jtheqn = 1,neq)
      CALL condis (nid, idi, beta, phi, nidm, ndisce, neq)
      DO kthnode = 1,numnp
         IF(id(kthnode)<0)THEN
            xnode(kthnode, ithmode) = phi(neq-id(kthnode))
         ELSEIF(id(kthnode)>0)THEN
            xnode(kthnode, ithmode) = phi(id(kthnode))
         ELSEIF(id(kthnode)==0)THEN
            xnode(kthnode, ithmode) = zero
         ENDIF
      ENDDO
      IF(mega==megamin)THEN
         WRITE(Lancos80_Unit) ithmode,(xnode(kthnode, ithmode), kthnode=1,numnp)
         WRITE(Lancos81_Unit) ithmode,(xnode(kthnode, ithmode), kthnode=1,numnp)
      ELSE
         WRITE(Lancos81_Unit) ithmode,(xnode(kthnode, ithmode), kthnode=1,numnp)
      ENDIF
   ENDDO
ENDIF
!-----------------------------------------------------------------------
! Determine the Surface function at each nodal point for sspace.
!-----------------------------------------------------------------------
REWIND Temp2_FEM_Unit
DO ithmode = 1,nqloc
   READ(Temp2_FEM_Unit) (phi(jtheqn),jtheqn = 1,neq)
   CALL condis (nid, idi, beta, phi, nidm, ndisce, neq)
   DO kthnode = 1,numnp
      IF(id(kthnode)<0)THEN
         xnode(kthnode, ithmode) = phi(neq-id(kthnode))
      ELSEIF(id(kthnode)>0)THEN
         xnode(kthnode, ithmode) = phi(id(kthnode))
      ELSEIF(id(kthnode)==0)THEN
         xnode(kthnode, ithmode) = zero
      ENDIF
   ENDDO
   IF(mega==megamin)THEN
      WRITE(ZroMeg_Unit) ithmode,(xnode(kthnode, ithmode), kthnode=1,numnp)
      WRITE(Nzromg_Unit) ithmode,(xnode(kthnode, ithmode), kthnode=1,numnp)
   ELSE
      WRITE(Nzromg_Unit) ithmode,(xnode(kthnode, ithmode),  kthnode=1,numnp)
   ENDIF
ENDDO
!
!-----------------------------------------------------------------------
! IF desired print out the surface functions.
!-----------------------------------------------------------------------
IF(full)THEN
   DO ithmode = 1,nqloc
      WRITE(Out_unit,*)' Surface function for ithmode = ',ithmode
      WRITE(Out_unit,*) (kthnode, xnode(kthnode, ithmode), kthnode = 1, numnp)
   ENDDO
ENDIF
RETURN

3 FORMAT(1x,f10.5)
4 FORMAT(1x,e14.7)
5 FORMAT(/,'$readx',5x,'distance in bohr')
6 FORMAT(1x,i5)
7 FORMAT(/,'$ready',5x,'energy in ev.')
8 FORMAT(/,'$CONTINUE')
END
