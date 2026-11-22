SUBROUTINE readfrst(ithrho, naph, rhoend, id, th, thetaval, ch, chivals, &
    phi, energy, nel, numnp, ntheta, nchi, mxnode, mxelem, lambda, lammin)
USE rhochk_Module
USE rhowrt_Module
USE restar_Module
USE fuzzy_Module
USE FileUnits_Module
USE Numbers_Module
!
!   $RCSfile: readfrst.f,v $ $Revision: 1.15 $
!   $Date: 89/11/16 11:03:25 $
!   $State: Stable $
!
!            P U R P O S E   O F    S U B R O U T I N E

!-----------------------------------------------------------------------
! This routine restarts the code.
! It is currently written assuming that there was a clean STOP
! to the code.  That is, you ran to 25 bohr and now want to
! run to 50 bohr.
!-----------------------------------------------------------------------
! input variables:
! ithrho  the number of rho values previously done.
! rhoend  the last rho value in the previous run.
!-----------------------------------------------------------------------
! Written by B.J. Archer
!-----------------------------------------------------------------------

!            I N P U T    A R G U M E N T S
!   ithrho
!   naph
!   rhoend
!   id
!   th
!   thetaval
!   ch
!   chivals
!   phi
!   energy
!   nel
!   numnp
!   ntheta
!   nchi
!   mxnode
!   mxelem
!   lambda
!   lammin

IMPLICIT NONE
LOGICAL logc
INTEGER i, ithrho, naphrd, naph, j, id, nelem, numnp, nlcrd, nel, ntheta, nchi, lambda, lammin, mxnode, mxelem, numnprd, th, ch
REAL(Kind=WP_Kind) rhocent, rhoend, rho1, rho2, thetaval, chivals, phi, energy
DIMENSION id(*), th(*), thetaval(*), ch(*), chivals(*), nlcrd(3), phi(*), energy(*)
EXTERNAL rdwrsurf
!-----------------------------------------------------------------------
!  This section is for restarting after ithrho rho values have
! been run and more are desired.
! The following files must be saved from the previous run. Note
! that you have to start at the next rho value after the last rho of
! the previous run.
! UNIT  FILE
! Pmatrx_FEM_Unit     Pmatrx_FEM
! 19    Ovrlp_FEM
! 22    elmdef
! Nmodes_FEM_Bin_Unit    nmodes
! Surf33_Unit      surf33
! sf34_unit        surf34
! OldMsh_Unit      oldmsh
! Zromsh_Bin_Unit  zromsh
! ZroMeg_Unit      zromeg
! Nzromg_Unit      nzromg
! PhiRms2_Unit     rmsphi2
! PhiRms1_Unit     rmsphi1
!-----------------------------------------------------------------------
! First reset Pmatrx to the END of the file.
!-----------------------------------------------------------------------
IF(rhoend<=zero)THEN
  WRITE(Out_unit,*)'rhoend must be set to restart'
  STOP 'readfrst'
ENDIF
lrestar=.true.
DO i=1,ithrho
  READ(Pmatrx_FEM_Unit)rhocent,naphrd
  READ(Pmatrx_FEM_Unit)
  DO j=1,3*naph
    READ(Pmatrx_FEM_Unit)
  ENDDO
ENDDO

IF(ABS(rhocent-rhoend)>fuzz)THEN
  WRITE(Out_unit,*)'Pmatrx not lined up, rhocent=',rhocent,' rhoend=', rhoend
  STOP 'readfr'
ENDIF

IF(naph/=naphrd)THEN
  WRITE(Out_unit,*)'wrong naph, naphrd=',naphrd,' naph=',naph
  STOP 'readfr'
ENDIF
!-----------------------------------------------------------------------
! READ to the END of the overlap file.
!-----------------------------------------------------------------------
DO 20 i=2,ithrho
  READ(OVR_FEM_unit)
  READ(OVR_FEM_unit)naphrd,rho1,rho2
  DO 25 j=1,naph
    READ(OVR_FEM_unit)
 25     CONTINUE
 20   CONTINUE
IF(naphrd/=naph)THEN
  WRITE(Out_unit,*)'wrong naph: naphrd=',naphrd,' naph=',naph,' Ovrlp'
  STOP 'readfr'
ENDIF
IF(ABS(rho2-rhoend)>fuzz)THEN
  WRITE(Out_unit,*)'Ovrlp not lined up, rho2=',rho2,' rhoend=',rhoend
  STOP 'readfr'
ENDIF
!-----------------------------------------------------------------------
! Find which surf file contains the correct surface function;
! the surface function for rhoend must be put in surf34.
!-----------------------------------------------------------------------
  CALL rdwrsurf(id,th,thetaval,ch,chivals,phi,energy,nelem,numnprd,naph,ntheta,nchi,34,1,mxnode,mxelem,lambda,lammin)
IF(ABS(rhor-rhoend)>fuzz)THEN
  CALL rdwrsurf(id,th,thetaval,ch,chivals,phi,energy,nelem,numnprd,naph,ntheta,nchi,33,1,mxnode,mxelem,lambda,lammin)
  IF(ABS(rhor-rhoend)<fuzz)THEN
    CALL rdwrsurf(id,th,thetaval,ch,chivals,phi,energy,nelem,numnprd,naph,ntheta,nchi,34,2,mxnode,mxelem,lambda,lammin)
  ELSE
    WRITE(Out_unit,*)'can not find the surface functions'
    STOP 'readfr'
  ENDIF
ENDIF
IF(numnprd/=numnp)THEN
  WRITE(Out_unit,*)'numnp not same, numnprd=',numnprd,' numnp=',numnp
  STOP 'readfr'
ENDIF
IF(nel/=nelem)THEN
  WRITE(Out_unit,*)'nelem not same, nel=',nel,' nelem=',nelem
  STOP 'readfr'
ENDIF
!-----------------------------------------------------------------------
! IF rhoend is greater/equal to rhowrit, THEN we need to reset the
! nmodes and elmdef files so they are in position to recieve the next
! rho value.
!-----------------------------------------------------------------------
IF(rhoend>=rhowrit.or.ABS(rhoend-rhowrit)<fuzz)THEN
  DO 40 i=1,ithrho
    READ(Nmodes_FEM_Bin_Unit)numnprd,numnprd,numnprd,logc,numnprd,rhocent,rhocent,rhocent,rho1,numnprd,naphrd
    IF(naph/=naphrd)THEN
      WRITE(Out_unit,*)'nmodes, wrong naph, naphrd=',naphrd,' naph=',naph
      STOP 'readfr'
    ENDIF
    READ(Nmodes_FEM_Bin_Unit)
    READ(Nmodes_FEM_Bin_Unit)
    READ(Nmodes_FEM_Bin_Unit)
    READ(Nmodes_FEM_Bin_Unit)
    READ(Nmodes_FEM_Bin_Unit)
    DO j=1,naphrd
      READ(Nmodes_FEM_Bin_Unit)
    END DO
    READ(ElmDef_Unit)numnprd,numnprd,numnprd,nlcrd,numnprd,rho2
    READ(ElmDef_Unit)
    READ(ElmDef_Unit)
    READ(ElmDef_Unit)
    IF(ABS(rho1-rho2)>fuzz)THEN
      WRITE(Out_unit,*)'elmdef and nmodes DO not agree',rho1,rho2
      STOP 'readfr'
    ENDIF
    IF(ABS(rho1-rhoend)<=fuzz)GOTO 50
 40     CONTINUE
  WRITE(Out_unit,*)'did not find rhoend in nmodes/elmdef', rho1
  STOP 'readfr'
 50     CONTINUE
ENDIF
!--------------------------------------------------------------------
RETURN
!----------------------***END-readfrst-------------------------------
END
