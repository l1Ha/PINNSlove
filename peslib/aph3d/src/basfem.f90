SUBROUTINE basfem (theta, chi, npoints, phi, nmode, rho, ithrho)
USE FileUnits_Module
USE Narran_Module
USE STST_Module
USE Rhosur_Module
USE Masses_Module
USE Convrsns_Module
USE TotJ_Module
USE Opts_Module
USE Crays_Module
USE CSBasis_Module
USE DoEnlvls_Module
USE Point_rsf_Module
USE Point_Msher_Module
IMPLICIT NONE
LOGICAL :: there
REAL(Kind=WP_Kind) :: massa, massb, massc, theta, chi, phi, rho, fuzz
INTEGER :: npoints, nmode, ithrho
INTEGER :: numnp, ndisce, neq, nel, nlambd, nidmx
DIMENSION chi(npoints), theta(npoints), phi(1)
OPEN(Unit=Msher_Bin_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//'.options')
READ(Msher_Bin_Unit, options)

CLOSE(unit=Msher_Bin_Unit)
rhosurf = rho
WRITE(Msg_Unit, * ) 'Construct FEM surface functions at quadrature ', 'points'
!-------------------------------------------------------------------
! READ in surface function DATA.
!-------------------------------------------------------------------
INQUIRE(file = OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/elmdef', exist = there)
IF(there)THEN
   OPEN(Unit=ElmDef_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Elmdef_FEM.bin', form = 'unformatted', status = 'old')
ELSE
   WRITE(Msg_Unit, * ) 'Error: elmdef does not exist'
   STOP 'basfem'

ENDIF
INQUIRE (file = OutDIR(1:LEN(TRIM(OutDIR)))//'nmodes', exist = there)
IF(there)THEN
   OPEN(Unit=Nmodes_FEM_Bin_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BibOut/Nmodes_FEM.bin', form = 'unformatted', status  = 'old')
ELSE
   WRITE(Msg_Unit, * ) 'Error: nmodes does not exist'
   STOP 'basfem'

ENDIF

    1 CONTINUE
READ(Nmodes_FEM_Bin_Unit, END = 91)jtot,mega,parity,symmetry,numnp,massa,massb,massc,rho,nel,nmode,ithrho,nlambd,ndisce,nidmx,neq,betaw,seq,chncof
WRITE(Out_Unit,  * )'rhoread=', rho, '   rhosurf=', rhosurf
WRITE(Out_Unit, * )jtot,mega,parity,symmetry,numnp,massa,massb,massc,rho,nel,nmode,ithrho,nlambd,ndisce,nidm,neq,betaw,seq,chncof
mass(1) = massa / amutoau
mass(2) = massb / amutoau
mass(3) = massc / amutoau
CALL massfac (.true.)
!-----------------------------------------------------------------------
!     Determine the array-positions in the w-array
!-----------------------------------------------------------------------
IF(nel>=0)THEN
   ALLOCATE(th(9,nel))
   ALLOCATE(ch(9,nel))
   ALLOCATE(id(9,int((nel+1)/2)))
ELSE
   WRITE(Out_Unit, * ) 'Error nel must be greater than or equal to 0'
   STOP 'basfem'
ENDIF
IF(numnp>0.AND.nmode>0)THEN
   ALLOCATE(f(numnp ,nmode))
   ALLOCATE(freq(nmode))
   ALLOCATE(xn(numnp))
ELSE
WRITE(Out_Unit,  * ) 'Error numnp and nmode  must be greater than 0'
   STOP 'basfem'
ENDIF
IF(nidm>0.AND.ndisce>0)THEN
   ALLOCATE(beta(nidm,ndisce))
   ALLOCATE(idi(nidm,int((ndisce+1)/2)))
   ALLOCATE(nid(int((ndisce+1)/2)))
   ALLOCATE(idni(int((numnp+1)/2)))
ELSE
   WRITE(Out_Unit, * ) 'Error nidm and ndisce must be greater than 0 '
   STOP 'basfem'
ENDIF

CALL rsffem (th, ch, nel, id, f, beta, idi, nid, nidm, freq, xn, idni, nmode, numnp, ndisce, neq, symmetry)
fuzz = 1.e-7
!---------------------------------------------------------------------
!     check to make sure that the correct rho value has been READ.
!     IF not go back and READ in another rho value.
!---------------------------------------------------------------------
IF(rhosurf > rho .AND. ABS(rhosurf-rho) > fuzz)THEN
   DEALLOCATE(th)
   DEALLOCATE(ch)
   DEALLOCATE(id)
   DEALLOCATE(f)
   DEALLOCATE(freq)
   DEALLOCATE(xn)
   DEALLOCATE(beta)
   DEALLOCATE(idi)
   DEALLOCATE(nid)
   DEALLOCATE(idni)
   GOTO 1
ENDIF

IF(rhosurf<rho.AND.abs (rhosurf - rho) >fuzz)THEN
   WRITE(Out_Unit, 7) rhosurf, rho
    7 FORMAT   (/, '***error***: rhosurf < rho', 2e14.7)
   STOP 'basfem'
ENDIF

rhosurf = rho
!---------------------------------------------------------------------
! Calculate FEM surface functions at the desired chi and theta points.
!---------------------------------------------------------------------
CALL intfem (npoints, chi, theta, phi, th, ch, id, nmode, f, numnp, nel)
DEALLOCATE(th)
DEALLOCATE(ch)
DEALLOCATE(id)
DEALLOCATE(f)
DEALLOCATE(freq)
DEALLOCATE(xn)
DEALLOCATE(beta)
DEALLOCATE(idi)
DEALLOCATE(nid)
DEALLOCATE(idni)
CLOSE(ElmDef_Unit)
CLOSE(Nmodes_FEM_Bin_Unit)
WRITE(Msg_Unit, * ) 'FEM surface functions constructed'
RETURN
   91 WRITE(Out_Unit, * ) 'END of surface function file, first'
STOP 'basfem'
ENDSUBROUTINE basfem
