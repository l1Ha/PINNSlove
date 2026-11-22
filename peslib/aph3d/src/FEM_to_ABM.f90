SUBROUTINE FEM_to_ABM
USE FileUnits_Module
USE opts_Module
USE rhos_Module
USE Narran_Module
USE Masses_Module
USE Numbers_Module
USE chiang_Module
USE Parms_Module
USE convrsns_Module
USE converge_Module
USE totj_Module
USE quantb_Module
USE gaussb_Module
USE Integopt_Module
USE AzBzCz_Module
USE Boundary_Module
!
!   $RCSfile: abf_main.f,v $ $Revision: 1.13 $
!   $Date: 89/10/09 11:19:06 $
!   $State: Stable $
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
!     phinow   present primitive basis at present quadrature points.
!     philast  previous rho prim basis at  present quadrature points.
!     vect     contains the coeffs of the surf fcns in the prim. basis.
!      unit 5  is used for NAMELIST input.
!      unit 6  is used for standard output.
IMPLICIT NONE
LOGICAL :: there

!            C H A R A C T E R    S T R I N G S
CHARACTER(LEN=10) today, hour, curzone
INTEGER dtvalues(8)
CHARACTER(LEN=200) :: filename, prefix
INTEGER :: nbasis, nangch, nang, maxn, maxl, neigmin, nbasiss, lnblnk
REAL(Kind=WP_Kind) :: rholast, vect, amatrx, eigen, philast, phinow, s, ThetAPH, chiaph, cpu1, cpu2
REAL(Kind=WP_Kind) :: overhead, tovr, weight
!            D I M E N S I O N S

DIMENSION nangch (narran + 1), philast (mxang, mxbasis), phinow ( mxang, mxbasis), ThetAPH (mxang), chiaph (mxang)
DIMENSION s (mxbasis, mxbasis), vect (mxbasis, mxbasis), amatrx (mxbasis, mxbasis), eigen (mxbasis, 0:mxmega)
DIMENSION nbasis (0:mxmega), neigmin (0:mxmega), nbasiss (0:mxmega), weight (mxang)
!            E X T E R N A L S

EXTERNAL constant, cputime
!-----------------------------------------------------------------------
! UNIT=Out_Unit sfunout is the standard output file.
!-----------------------------------------------------------------------
!OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/FEM_to_ABM.txt', form = 'formatted', status = 'unknown')
!-----------------------------------------------------------------------
! UNIT=Ovr_Unit Ovrlp is an unformatted file for storing the overlap
!        matrix elements.
!-----------------------------------------------------------------------

WRITE(Out_Unit,*)'Called    FEM_To_ABM'
CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/FEM_To_ABM.txt',Form='FORMATTED')
OPEN(Unit=Ovr_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Ovrlp_FEM_to_ABM.bin', form = 'unformatted',status='unknown')
!-------------------------------------------------------------
!-------------------------------------------------------------
! Determine the time and todays date.
!-------------------------------------------------------------
CALL Date_And_Time(today, hour, curzone, dtvalues)
!-------------------------------------------------------------
!   WRITE the time and date to unit 6.
!-------------------------------------------------------------
CALL Print_Date_And_Time(Today, Hour, CurZone, DtValues, Out_Unit)
!-----------------------------------------------------------------------
! Establish all conversion factors and physical constants
!-----------------------------------------------------------------------

CALL constant
!-----------------------------------------------------------------------
! Determine overhead time for CALL timing routine.
!-----------------------------------------------------------------------
CALL cputime (cpu1)
CALL cputime (cpu2)

overhead = cpu2 - cpu1
INQUIRE(file = OutDIR(1:LEN(TRIM(OutDIR)))//'setbasis', exist = there)
IF(there)THEN
   OPEN(Unit=Setbasis_Bin_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Setbasis.bin', form = 'unformatted', status = 'old')
   READ(Setbasis_Bin_Unit) nang, nbasis, maxn, maxl, nangch, nbasiss, neigmin
   READ(Setbasis_Bin_Unit) minvib, jmin, maxvib, jmax, chanl, nvib, jrot, jskip
   READ(Setbasis_Bin_Unit) weau, ralpha, re, rx, noscil, nhermt, nglegn, nlegndre, balpha, calpha, dalpha, anharm, wexeau, intwt
   READ(Setbasis_Bin_Unit) jtot, parity, mega, megamin, megamax, symmetry, jeven
   READ(Setbasis_Bin_Unit) rhoabm, scheme
   rho = rhoabm
   CLOSE(unit=Setbasis_Bin_Unit)
ELSE
   WRITE(Msg_Unit, * ) 'Error: setbasis does not exist'
   STOP 'FEM_to_ABM'

ENDIF
!-----------------------------------------------------------------------
! Loop over surface functions for each mega.
!-----------------------------------------------------------------------
DO 4 mega = megamin, megamax
!-----------------------------------------------------------------------
! Calculate the Surface functions at rho.
!-----------------------------------------------------------------------
   CALL readbas(rho, nbasis(mega), mega, nang, phinow, ThetAPH, chiaph, weight, vect, neigmin(mega), eigen(1,mega), nbasiss(mega), IthStart_ABM)
   WRITE(Out_unit, * ) 'mega, nbasis, nbasiss', mega, nbasis (mega) , &
         nbasiss (mega)
   CALL weight_unset (nang, nbasiss (mega), phinow, weight)

   rhoabm = rho
!-----------------------------------------------------------------------
! Calculate overlap integrals of the surface functions at rho with
! the surface functions at the previous rho (rholast).
!-----------------------------------------------------------------------
   CALL cputime (cpu1)
   INQUIRE(file = '.FEM_to_ABM', exist = there)
   IF(there)THEN
      filename = '.FEM_to_ABM'
   ELSE
      CALL getenv ('APH3D_HOME', prefix)
      filename = prefix (:lnblnk (prefix) ) //'/default_data/.FEM_to_ABM '
   ENDIF
   
   OPEN(Unit=FEM_TO_ABM_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//filename, status = 'old')
   READ(FEM_TO_ABM_Unit, integopt)
   IF(int_opt=='ABM_Integrals')THEN
      WRITE(Msg_Unit, * ) 'Using ABM Quadrature'
      WRITE(Out_unit, * ) 'Using ABM Quadrature'
   ELSEIF (int_opt=='FEM_Integrals')THEN
      WRITE(Msg_Unit, * ) 'Using FEM Quadrature'
      WRITE(Out_unit, * ) 'Using FEM Quadrature'
   ELSE
      WRITE(Msg_Unit, * ) 'Must use FEM or ABM Integrals'
      STOP 'main'
   ENDIF

   CLOSE(unit=FEM_TO_ABM_Unit)
   WRITE(Out_unit, * ) 'mega, nbasis, nbasiss', mega, nbasis (mega) ,&
                        nbasiss (mega)
   IF(int_opt=='ABM_Integrals')THEN
      rholast = rhofem
      CALL ovrfem (rholast, rho, nang, nbasis (mega), chanl (1,  mega), &
           nangch, philast, phinow, vect, s, amatrx, nbasiss (mega), &
           neigmin (mega), jrot (1, mega), ThetAPH, chiaph,  weight)
   ELSEIF (int_opt=='FEM_Integrals')THEN
      rholast = rhofem
      CALL ovrbaf (rholast, rho, nbasis(mega), narran, vect, s, amatrx, neigmin(mega), megamin, mega, nbasiss(mega) )
   ELSE
      WRITE(Msg_Unit, * ) 'Must use FEM or ABM Integrals'
      STOP 'main'

   ENDIF
   rhofem = rholast
WRITE(Msg_Unit,  * ) 'Rho_FEM=', rhofem, '   Rho_ABM=', rhoabm

WRITE(Out_unit,  * ) 'Rho_FEM=', rhofem, '   Rho_ABM=', rhoabm
WRITE(Out_unit,  * ) 'For mega =', mega, '  the overlap matrix is:'

   CALL MxOut(s, nbasiss (mega), nbasiss (mega) )
   CALL cputime (cpu2)
   tovr = cpu2 - cpu1 - overhead

   WRITE(Out_unit, * ) 'Time in ovrfem = ', tovr


    4 ENDDO
WRITE(Out_unit, * ) ' rholast= ', rholast
WRITE(Out_unit, * ) ' neigmin=', (neigmin (mega) , mega = megamin, megamax)


CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/BasisForScatt.txt',Form='FORMATTED',ACCESS='APPEND')
WRITE(Out_Unit,*)'Completed FEM_To_ABM'
WRITE(Out_Unit,*)

RETURN
ENDSUBROUTINE FEM_to_ABM
