SUBROUTINE readbas(rho, nbasis, megaval, nang, phinow, ThetAPH, chiaph, weight, vect, neigmin, eigen, nbasiss, IthRho)
USE Parms_Module
USE Convrsns_Module
USE Converge_Module
USE AzBzCz_Module
USE FileUnits_Module
USE Boundary_Module
USE TotJ_Module
USE Numbers_Module
!
!   $RCSfile: readbas.f,v $ $Revision: 1.13 $
!   $Date: 89/10/09 11:19:20 $
!   $State: Stable $
!
!            P U R P O S E   O F    S U B R O U T I N E
!   This routine calculates the surface functions using an analytic
!   basis set.


!            I N P U T    A R G U M E N T S
!  rho      current hyperradius.
!  mega     projection of the total angular momentum on the body-fixed
!           z-axis.
!  nbasis   the total number of raw primitive basis functions.
!  nbasiss  the number of symmetry-adapted primitive basis fcns.
!           IF symmetry=.false., nbasiss=nbasis.

!            O U T P U T    A R G U M E N T S
!  phinow   primitive basis at the quadrature points.
!  vect     matrix of coeffs of surface fcns expanded in prim. basis.


IMPLICIT NONE
LOGICAL there
INTEGER megaval, nbasis, nang, ibig, ibasis, iang, jbasis, neigmin, nbasiss, i, nangx, IthRho
REAL(Kind=WP_Kind) rho, big
REAL(Kind=WP_Kind) eigen(neigmin), vect(nbasiss, neigmin), phinow(nang, nbasis)
REAL(Kind=WP_Kind) ThetAPH(nang), chiaph(nang), weight(nang)
EXTERNAL popt, mxout
!-----------------------------------------------------------------------
! Determine printing options.
!-----------------------------------------------------------------------
LOGICAL little, medium, full
INTEGER ithcall, ithsub
DATA ithcall /0/, ithsub /0/
DATA little /.false./, medium /.false./, full /.false./
CALL popt ('readbas ', little, medium, full, ithcall, ithsub)
Full=.False.
IF(Mega/=MegaVal)THEN
   WRITE(Out_Unit,*)"Mega/=MegaVal:",Mega,MegaVal
   WRITE(Msg_Unit,*)"Mega/=MegaVal:",Mega,MegaVal
   STOP "Stopping in ReadBas"
ENDIF
!  --------------------------------------------------------------
!  READ vect(nbasiss,neigmin) on unit Vfrst_Unit ='vector.first'
!  --------------------------------------------------------------
IF(IthRho==IthStart_ABM)THEN
   INQUIRE(file=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Vector_First.bin',exist=there)
   IF(there)THEN
      OPEN(Unit=Vfrst_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Vector_First.bin',form='unformatted', status='old')
      READ(Vfrst_Unit) nbasiss, neigmin, rho
      DO jbasis=1,neigmin
         READ(Vfrst_Unit)(vect(ibasis,jbasis), ibasis=1,nbasiss)
      ENDDO
      CLOSE(Unit=Vfrst_Unit,status='keep')
   ELSE
      WRITE(Msg_Unit,*)'Error: Vector_First.bin does not exist'
      STOP 'readbas'
   ENDIF
ELSE
   INQUIRE(file=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Vector_ABM.bin',exist=there)
   IF(there)THEN
      OPEN(Unit=Vec_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Vector_ABM.bin',form='unformatted', status='old')
      READ(Vec_Unit) nbasiss, neigmin, rho
      DO jbasis=1,neigmin
         READ(Vec_Unit)(vect(ibasis,jbasis), ibasis=1,nbasiss)
      ENDDO
      CLOSE(Unit=Vec_Unit,status='keep')
   ELSE
      WRITE(Msg_Unit,*)'Error: Vector_ABM.bin does not exist'
      STOP 'readbas'
   ENDIF
ENDIF

!  --------------------------------------------------------------
! UNIT=Pfrst_Unit=phiprim.first is an unformatted file for storing the primitive basis.
!-----------------------------------------------------------------------
IF(IthRho==IthStart_ABM)THEN
   INQUIRE(file=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Phiprim_First.bin',exist=there)
   IF(there)THEN
      IF(mega==megamin)THEN
         OPEN(Unit=Pfrst_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Phiprim_First.bin', form='unformatted',status='old')
         REWIND Pfrst_Unit
      ENDIF
      READ(Pfrst_Unit) mega, nangx, nbasis
      WRITE(Out_Unit,*)' mega, nang, nbasis = ', mega, nang, nbasis
      DO jbasis=1, nbasis
         READ(Pfrst_Unit) (phinow(iang,jbasis), iang=1, nang) !XXXX
      ENDDO
      IF(mega==megamax)THEN
         CLOSE(Unit=Pfrst_Unit,status='keep')
      ENDIF
   ELSE
      WRITE(Msg_Unit,*)'Error: phiprim.first does not exist'
      STOP 'readbas'
   ENDIF
ELSE
   INQUIRE(file=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/PhiPrim_ABM.bin',exist=there)
   IF(there)THEN
      OPEN(Unit=Phip_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/PhiPrim_ABM.bin', form='unformatted',status='old')
      REWIND Phip_Unit
      READ(Phip_Unit) mega, nangx, nbasis
      WRITE(Out_Unit,*)' mega, nang, nbasis = ', mega, nang, nbasis
      DO jbasis=1, nbasis
         READ(Phip_Unit) (phinow(iang,jbasis), iang=1, nang)
      ENDDO
      !CLOSE(Unit=Phip_Unit,status='keep')
   ELSE
      WRITE(Msg_Unit,*)'Error: PhiPrim_ABM.bin does not exist'
      STOP 'readbas'
   ENDIF
ENDIF

IF(IthRho==IthStart_ABM)THEN
   INQUIRE(file=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Weights_First.bin',exist=there)
   WRITE(Out_Unit,*)'nang=',nang
   IF(there)THEN
      OPEN(Unit=wfrst_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Weights_First.bin', form='unformatted', status='old')
      READ(wfrst_Unit)(ThetAPH(iang),iang=1,nang)
      READ(wfrst_Unit)(chiaph(iang),iang=1,nang)
      READ(wfrst_Unit)(weight(iang),iang=1,nang)
      READ(wfrst_Unit)(eigen(i),i=1,neigmin)
      READ(wfrst_Unit)az, bz, cz
      WRITE(Out_Unit,*)'ABM Eigenenergies'
      IF(medium) WRITE(Out_Unit,*)eigen
      CLOSE(Unit=wfrst_Unit,status='keep')
   ELSE
      WRITE(Msg_Unit,*)'Error: weights.first does not exist'
      STOP 'readbas'
   ENDIF
ELSE
   INQUIRE(file=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Weights.bin',exist=there)
   WRITE(Out_Unit,*)'nang=',nang
   IF(there)THEN
      OPEN(Unit=Weight_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Weights.bin', form='unformatted', status='old')
      READ(Weight_Unit)(ThetAPH(iang),iang=1,nang)
      READ(Weight_Unit)(chiaph(iang),iang=1,nang)
      READ(Weight_Unit)(weight(iang),iang=1,nang)
      READ(Weight_Unit)(eigen(i),i=1,neigmin)
      READ(Weight_Unit)az, bz, cz
      WRITE(Out_Unit,*)'ABM Eigenenergies'
      IF(medium) WRITE(Out_Unit,*)eigen
      CLOSE(Unit=Weight_Unit,status='keep')
   ELSE
      WRITE(Msg_Unit,*)'Error: weights.first does not exist'
      STOP 'readbas'
   ENDIF
ENDIF

IF(medium)THEN
   WRITE(Out_Unit,*)'routine readbas'
   WRITE(Out_Unit,'(A,Es15.7,3I5)')'rho, mega = ',rho, mega
   WRITE(Out_Unit,'(A,3I5,A,ES15.7)')'nang, nbasiss, neigmin = ', nang,nbasiss,neigmin,' overlap eigmin = ', eigmin
ENDIF

IF(full)THEN
   WRITE(Out_Unit,*)'nbasiss,neigmin=',nbasiss,neigmin
   DO jbasis = 1, neigmin
      ibig=0
      big = 0.0d0
      DO ibasis = 1, nbasiss
         IF(ABS(vect(ibasis,jbasis))>big)THEN
             ibig = ibasis
             big = vect(ibasis,jbasis)
         ENDIF
      ENDDO
      WRITE(Out_Unit,*)'jbasis,ibig,big',jbasis,ibig,big
      WRITE(Msg_Unit,*)'jbasis,ibig,big',jbasis,ibig,big
   ENDDO
   WRITE(Out_Unit,*)'eigenvectors'
   CALL MxOut(vect, nbasiss, nbasiss)
   WRITE(Out_Unit,*)'primitive basis'
   CALL MxOut(phinow, nang, nbasis)
ENDIF

RETURN
ENDSUBROUTINE readbas
