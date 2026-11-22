SUBROUTINE basbaf (nbasis, ndvr, neigmin, vect, amatrx, rho, ithrho, test, nbasiss)
USE Numeric_Kinds_Module
USE FileUnits_Module
USE DVR_Module
USE DVR2_Module
USE Convrsns_Module
USE Opts_Module
USE RhoSur_Module
USE Crays_Module
USE Csbasis_Module
USE DoEnlvls_Module
USE Totj_Module
USE FBR_Module
USE Point_RSF_Module
USE Parms_Module
USE InputFile_Module
IMPLICIT NONE
CHARACTER(LEN=3) test
LOGICAL there
INTEGER nbasis, neigmin, ithrho, nbasiss, ndvr, nabm, Nang_DVR, naph, index, itheta, ichi
INTEGER IErr
REAL(Kind=WP_Kind) amatrx(nbasiss,nbasiss), vect(nbasiss,nbasiss), rho
REAL(Kind=WP_Kind), ALLOCATABLE::  weightd(:)
REAL(Kind=WP_Kind), ALLOCATABLE::  chi1(:), thaph(:)
REAL(Kind=WP_Kind), ALLOCATABLE::  Phi_DVR(:,:), Phi_ABM(:,:)

WRITE(Out_Unit,*)'basbaf called:'
WRITE(Out_Unit,*)'     rho=', rho
WRITE(Out_Unit,*)'  nbasis=', nbasis
WRITE(Out_Unit,*)' neigmin=', neigmin
WRITE(Out_Unit,*)' nbasiss=', nbasiss
WRITE(Out_Unit,*)'ndvptmax=', ndvptmax
WRITE(Out_Unit,*)' nsfnmax=', nsfnmax
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile)
READ(In_Unit,NML=options)
CLOSE(Unit=In_Unit)
WRITE(Out_Unit,NML=options)

ALLOCATE(weightd(ndvptmax))
ALLOCATE(chi1(ndvptmax))
ALLOCATE(thaph(ndvptmax))
ALLOCATE(Phi_DVR(ndvptmax,nsfnmax))
ALLOCATE(Phi_ABM(ndvptmax,nsfnmax))

rhosurf = rho

!----------------------------------------------------------------------
!  READ in parameters of fbr(dvr) surface functions to be used.
!  lmax+1 is the actual ntheta used in the dvr calculations.
!  mmax is the actual nchi.
!  Nang_DVR is the total number of dvr points used.
!  naph=nsfunc is the number of dvr surface functions used.
!----------------------------------------------------------------------
INQUIRE(file=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,exist=there)
IF(there)THEN
  OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
  READ(In_Unit,NML=fbr, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'FBR')
      NSfunc=105 !TmpmodGregParker1
      lmax=20
      mmax=40
      nsym=1
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist FBR'
      WRITE(Msg_Unit,*)'ERROR with Namelist FBR'
      STOP 'BasBaf: FBR'
   ENDIF
  CLOSE(Unit=In_Unit)
  WRITE(Out_Unit,NML=fbr)
  naph = nsfunc
  Nang_DVR = mmax*(lmax+1)
  WRITE(Out_Unit,*)'    naph=', naph
  WRITE(Out_Unit,*)'Nang_DVR=', Nang_DVR
ELSE
  WRITE(Msg_Unit,*)'Error InputFile does not exist', InputDIR(1:LEN(TRIM(InputDIR)))//InputFile
  STOP 'basbaf: No Input Data'
ENDIF
!----------------------------------------------------------------------
! check array sizes for fbr(dvr) surface functions
!----------------------------------------------------------------------
IF(lmax+1>nthetmax)THEN
  WRITE(Out_Unit,*)' lmax + 1 (.fbrdata) > nthetmax (dvr.parms) '
  WRITE(Out_Unit,*)' make nthetmax bigger '
  STOP 'basbaf: lmax+1>nthetmax'
ENDIF
IF(mmax>nchimax)THEN
  WRITE(Out_Unit,*)' mmax (.fbrdata) > nchimax (dvr.parms) '
  WRITE(Out_Unit,*)' make nchimax bigger '
  STOP 'basbaf: mmax>nchimax'
ENDIF
!
IF(Nang_DVR>ndvptmax)THEN
  WRITE(Out_Unit,*) 'Nang_DVR =',Nang_DVR,' > ','ndvptmax =',ndvptmax
  STOP 'basbaf: Nang_DVR>ndvptmax'
ENDIF

index=0
DO itheta=1,lmax+1
   DO ichi=1,mmax
      index=index+1
      thaph(index)=acos(pttheta(itheta))
      chi1(index)=acos(ptchi(ichi))
      weightd(index)=1.d0       !Temp FIX TEMP tempmod
   ENDDO
ENDDO
      

!---------------------------------------------------------------------
!  construct the fbr(dvr) functions at the dvr quadrature points.
!---------------------------------------------------------------------

CALL intdvr(Nang_DVR, chi1, thaph, Phi_DVR, naph, rho)
rhosurf = rho
WRITE(Out_Unit,*) 'intdvr done',nbasiss

!----------------------------------------------intdvr done-------------------------

IF(test=='ABM')THEN
   ndvr = nbasiss
   nabm = nbasiss
ELSEIF(test=='DVR')THEN
   ndvr = naph
   nabm = naph
ELSE
   ndvr = naph
   nabm = nbasiss
ENDIF
!---------------------------------------------------------------------
! calculate abm primitive basis functions at dvr nodal points.
! construct overlap matrix.
!---------------------------------------------------------------------

CALL intbaf (Nang_DVR, ndvr, nabm, Phi_DVR, Phi_ABM, amatrx, weightd, thaph, chi1, test, nbasis, neigmin, nbasiss)

DEALLOCATE(weightd)
DEALLOCATE(chi1)
DEALLOCATE(thaph)
DEALLOCATE(Phi_DVR)
DEALLOCATE(Phi_ABM)

WRITE(Out_Unit,*)' Finished calculating overlaps'

RETURN
 91   WRITE(Out_Unit,*)'END of surface function file, first'
STOP 'basbaf'
ENDSUBROUTINE basbaf
