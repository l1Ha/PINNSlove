SUBROUTINE DVR_to_ABM
USE FileUnits_Module
USE Opts_Module
USE Parms_Module
USE Masses_Module
USE Numbers_Module
USE ChiAng_Module
USE Convrsns_Module
USE Converge_Module
USE TotJ_Module
USE Quantb_Module
USE Gaussb_Module
USE AzBzCz_Module
USE Integopt_Module
USE InputFile_Module
USE Narran_Module
USE Boundary_Module
!
!            P U R P O S E   O F    S U B R O U T I N E
!     This program calculates the surface functions using
!     analytic basis functions. It also calculates overlap
!     matrix elements between the surfaces functions at the
!     previous rho with the surface functions at the current
!     rho. THEN potential energy interaction matrix elements
!     are calculated.

!           V A R I A B L E S
!     today     CHARACTER variable containing the current date.
!     hour      CHARACTER variable containing the current time.
!     naph      number of surface functions to be calculated.
!     jtot      total angular momentum.
!     mega      Projection of the total angular momentum on
!               the body-fixed z-axis.
!     megamin   minimum value of mega.
!     megamax   maximum value of mega.
!     parity    parity of the surface functions.
!     amatrx   amatrx and s are scratch matrices of
!              dimension mxbasis x mxbasis.
!     Phi_ABM   present primitive basis at present quadrature points.
!     Phi_DVR  previous rho prim basis at  present quadrature points.
!     vect     contains the coeffs of the surf fcns in the prim. basis


IMPLICIT NONE
LOGICAL there

!            C H A R A C T E R    S T R I N G S
CHARACTER(LEN=10) hour
CHARACTER(LEN=8) today
CHARACTER(LEN=5) curzone
INTEGER dtvalues(8)
CHARACTER(LEN=200) filename, prefix

!            I N T E G E R S
INTEGER nang, maxn, maxl, lnblnk

!            R E A L S
REAL(Kind=WP_Kind) Rho
REAL(Kind=WP_Kind) rholast,  cpu1, cpu2, overhead, tovr

!            D I M E N S I O N S
INTEGER, ALLOCATABLE:: nbasis(:), neigmin(:), nbasiss(:), nangch(:)
REAL(Kind=WP_Kind), ALLOCATABLE:: Phi_DVR(:,:), Phi_ABM(:,:),ThetAPH(:), chiaph(:), s(:,:)
REAL(Kind=WP_Kind), ALLOCATABLE:: vect(:,:), amatrx(:,:), eigen(:,:), weight(:)

!            E X T E R N A L S
EXTERNAL constant, cputime
!-----------------------------------------------------------------------
! UNIT=Ovr_Unit Ovrlp is an unformatted file for storing the overlap
!        matrix elements.
!-----------------------------------------------------------------------
WRITE(Out_Unit,*)'Called    DVR_To_ABM'
CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/BasisForScattDVR_To_ABM.txt',Form='FORMATTED')
OPEN(Unit=Ovr_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Ovrlp_DVR_to_ABM.bin', form='unformatted', status='unknown')

!-------------------------------------------------------------
! Determine the time and todays date.
!-------------------------------------------------------------

CALL Date_And_Time(today, hour, curzone, dtvalues)
CALL Print_Date_And_Time(Today, Hour, CurZone, DtValues, Out_Unit)
!-----------------------------------------------------------------------
! Establish all conversion factors and physical constants
!-----------------------------------------------------------------------
CALL constant

!-----------------------------------------------------------------------
! Determine overhead time for CALL timing routine.
!-----------------------------------------------------------------------
CALL cputime(cpu1)
CALL cputime(cpu2)
overhead = cpu2 - cpu1
!-----------------------------------------------------------------------
!  READ setbasis.first which contains basis used by abm at its first rho.
!-----------------------------------------------------------------------

INQUIRE(file=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/SetBasis_First.bin',exist=there)
ALLOCATE(NAngch(Narran+1))
IF(there)THEN
   OPEN(Unit=SBFrst_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/SetBasis_First.bin',form='unformatted', status='old')
   READ(SBFrst_Unit)MxAng,MxBasis,mxbasiss,MxMega
   WRITE(Out_Unit,*)'MxAng,MxBasis,mxbasiss,MxMega=',MxAng,MxBasis,mxbasiss,MxMega
   ALLOCATE(nbasis(0:MxMega))
   ALLOCATE(neigmin(0:MxMega))
   ALLOCATE(nbasiss(0:MxMega))
   READ(SBFrst_Unit)nang,nbasis,maxn,maxl,nangch,nbasiss,neigmin
   WRITE(Out_Unit,'(1x,A,5I5)')'MxAng,NAng,MxBasis,mxbasiss,MxMega=',MxAng,NAng,MxBasis,mxbasiss,MxMega
   ALLOCATE(Phi_DVR(Nang, MxBasis))
   ALLOCATE(Phi_ABM(Nang, MxBasis))
   ALLOCATE(s(MxBasis,MxBasis))
   ALLOCATE(vect(MxBasis, MxBasis))
   ALLOCATE(amatrx(MxBasis,MxBasis))
   ALLOCATE(ThetAPH(Nang))
   ALLOCATE(chiaph(Nang))
   ALLOCATE(weight(Nang))
   ALLOCATE(eigen(MxBasis, 0:MxMega))
   ALLOCATE(chanl(MxBasis, 0:MxMega), nvib(MxBasis), jrot(MxBasis,0:MxMega))
   READ(SBFrst_Unit)minvib, jmin, maxvib, jmax, chanl, nvib, jrot, jskip
   READ(SBFrst_Unit)weau, ralpha, re, rx, noscil, nhermt, nglegn, nlegndre, balpha, calpha,  dalpha, anharm, wexeau, intwt
   READ(SBFrst_Unit)jtot, parity, mega, megamin, megamax,  symmetry, jeven
   READ(SBFrst_Unit)rhoabm,scheme
   CLOSE(Unit=SBFrst_Unit)
ELSE
   WRITE(Msg_Unit,*)'Error: SetBasis_First.bin does not exist'
   STOP 'DVR_to_ABM'
ENDIF


!
!-------------------------------------------------------------------------
!  check size of nang
!-------------------------------------------------------------------------

IF(nang > mxang )THEN
  WRITE(Out_Unit,*) 'nang =',nang,' > mxang =',mxang
  !STOP 'Stopping in DVR_to_ABM'
ENDIF

!-----------------------------------------------------------------------
! Loop over surface functions for each mega.
!-----------------------------------------------------------------------
DO 4 mega = megamin, megamax
!-----------------------------------------------------------------------
! Calculate the Surface functions at rho.
!----------------------------------------------------------------------- !XXXX
   CALL readbas(rho, nbasis(mega), mega, nang, Phi_ABM, ThetAPH, chiaph, weight, vect, neigmin(mega), eigen(1,mega), nbasiss(mega), IthStart_ABM)
   WRITE(Out_Unit,*)'mega, nbasis, nbasiss', mega, nbasis(mega), nbasiss(mega)
   CALL weight_unset(nang,nbasiss(mega),Phi_ABM,weight)
   rhodvr=Basis_Dist(IthEnd_DVR)
   rhoabm=Basis_Dist(IthStart_ABM)
   rho=rhoabm
   rholast= rhodvr

!-----------------------------------------------------------------------
! Calculate overlap integrals of the surface functions at rho with
! the surface functions at the previous rho (rholast).
!-----------------------------------------------------------------------
   CALL cputime(cpu1)
   INQUIRE(file=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile, exist = there)
   IF(there)THEN
      filename='.DVR_to_ABM'
   ELSE
      CALL getenv('APH3D_HOME', prefix)
      filename = prefix(:lnblnk(prefix))//'/default_data/.DVR_to_ABM'
   ENDIF

   OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
   READ(In_Unit,NML=integopt)
   CLOSE(Unit=In_Unit)
   WRITE(Out_Unit,NML=integopt)
   IF(int_opt=='ABM_Integrals')THEN
      WRITE(Msg_Unit,*)'Using ABM Quadrature'
      WRITE(Out_Unit,*)'Using ABM Quadrature'
   ELSEIF(int_opt=='DVR_Integrals')THEN
      WRITE(Msg_Unit,*)'Using DVR Quadrature'
      WRITE(Out_Unit,*)'Using DVR Quadrature'
   ELSE
      WRITE(Msg_Unit,*)'Must use DVR or ABM Integrals'
      STOP 'main'
   ENDIF

   WRITE(Out_Unit,*)'mega, nbasis, nbasiss', mega, nbasis(mega),nbasiss(mega)
   IF(int_opt=='ABM_Integrals')THEN
      rholast = rhodvr
      CALL ovrdvr (rholast, rho, nang, nbasis(mega),chanl(1,mega),nangch, Phi_DVR, Phi_ABM, vect, s,amatrx, nbasiss(mega),   &
                   neigmin(mega), jrot(1,mega), ThetAPH, chiaph, weight )
   ELSEIF(int_opt=='DVR_Integrals')THEN
      WRITE(Msg_Unit,*)'Warning!  DVR integral option not debugged!'
      WRITE(Out_Unit,*)'Warning!  DVR integral option not debugged!'
      rholast = rhodvr
      CALL ovrbaf (rholast, rho, nbasis(mega), narran, vect, s, amatrx, neigmin(mega), megamin, mega, nbasiss(mega))
   ELSE
      WRITE(Msg_Unit,*)'Must use DVR or ABM Integrals'
      STOP 'main'
   ENDIF

   rhodvr = rholast
   WRITE(Msg_Unit,*)'Rho_DVR=',rhodvr,'   Rho_ABM=',rhoabm
   WRITE(Out_Unit,*)'Rho_DVR=',rhodvr,'   Rho_ABM=',rhoabm

   WRITE(Out_Unit,*)'For mega =',mega,'  the overlap matrix is:'
   CALL MxOut(s,nbasiss(mega),nbasiss(mega))

   CALL cputime(cpu2)
   tovr = cpu2 - cpu1 -overhead
   WRITE(Out_Unit,*)'Time in ovrdvr = ', tovr

   4  CONTINUE


WRITE(Out_Unit,*)' rholast= ', rholast
WRITE(Out_Unit,*)' neigmin=',(neigmin(mega), mega=megamin, megamax)


CLOSE(dvr7)
CLOSE(dvr11)
CLOSE(dvr12)
CLOSE(dvr13)
CLOSE(dvr14)
CLOSE(dvr15)
CLOSE(dvr16)
CLOSE(dvr18)
CLOSE(dvr23)
CLOSE(dvr22)
CLOSE(Ovr_Unit)
CLOSE(SBFrst_Unit)
CLOSE(SFunFBR_Unit)
CLOSE(Vfrst_Unit)
DEALLOCATE(nbasis)
DEALLOCATE(neigmin)
DEALLOCATE(nbasiss)
DEALLOCATE(Phi_DVR)
DEALLOCATE(Phi_ABM)
DEALLOCATE(s)
DEALLOCATE(vect)
DEALLOCATE(amatrx)
DEALLOCATE(ThetAPH)
DEALLOCATE(chiaph)
DEALLOCATE(weight)
DEALLOCATE(eigen)
DEALLOCATE(chanl, nvib, jrot)
DEALLOCATE(NAngch)


CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/BasisForScatt.txt',Form='FORMATTED',ACCESS='APPEND')
WRITE(Out_Unit,*)'Completed DVR_To_ABM'
WRITE(Out_Unit,*)

RETURN

END
