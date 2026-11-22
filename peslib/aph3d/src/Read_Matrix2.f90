SUBROUTINE Read_Matrix2(amat, energies, chanl, elect, nvib, jrot,   &
                        lorb, n, etot, ithmat, TypeOfMat, jtot,     &
                        proj, approx, symm, inunit, outunit,        &
                        nstate, more)
USE FileUnits_Module
USE Region_Module, ONLY: iregion
USE Oops_Module,   ONLY: ioops
IMPLICIT NONE
LOGICAL ifirst
INTEGER n, nstate, jtot, chanl(NState), elect(NState), nvib(NState), jrot(NState)
INTEGER  lorb(NState), ithmat, jthmat, nu_jtot, inunit, outunit
CHARACTER(LEN=20) TypeOfMat, nu_type, approx, nu_approx, proj, nu_proj
CHARACTER(LEN=20) symm, nu_symm
REAL(Kind=WP_Kind) etot, amat(NState,NState), energies(NState)
LOGICAL more
DATA ifirst/.true./

READ(inunit,END=6)jthmat, nu_type, n, etot, nu_jtot, nu_proj, nu_approx, nu_symm
!WRITE(msg_unit,*)'inunit=',inunit,' typeofmat= ',typeofmat
!WRITE(msg_Unit,*)'Read_Matrix2'
!WRITE(msg_Unit,*)'jthmat=',jthmat
!WRITE(msg_Unit,*)'nu_type=',nu_type
!WRITE(msg_Unit,*)'n=',n
!WRITE(msg_Unit,*)'etot=',etot
!WRITE(msg_Unit,*)'nu_jtot=',nu_jtot
!WRITE(msg_Unit,*)'nu_proj=',nu_proj
!WRITE(msg_Unit,*)'nu_approx=',nu_approx
!WRITE(msg_Unit,*)'nu_symm=',nu_symm

IF (ifirst)THEN
   ifirst=.false.
   TypeOfMat=nu_type
   jtot=nu_jtot
   proj=nu_proj
   approx=nu_approx
   symm=nu_symm
ENDIF

IF(jthmat/=ithmat+1)THEN
   WRITE(msg_unit,*)'Error: jthmat/=ithmat+1'
   WRITE(msg_unit,*)'ithmat=',ithmat
   WRITE(msg_unit,*)'jthmat=',jthmat
   STOP 'Read_Matrix2'
ENDIF

IF(TypeOfMat/=nu_type)THEN
   WRITE(msg_unit,*)'ERROR: TypeOfMat/=nu_type'
   WRITE(msg_unit,*)'TypeOfMat =',TypeOfMat
   WRITE(msg_unit,*)'nu_type=',nu_type
   STOP 'Read_Matrix2'
ENDIF

IF(n>nstate)THEN
   WRITE(msg_unit,*)'ERROR: n>nstate'
   WRITE(msg_unit,*)'n =',n
   WRITE(msg_unit,*)'nstate=',nstate
   STOP 'Read_Matrix2'
ENDIF

IF(jtot/=nu_jtot)THEN
   WRITE(msg_unit,*)'ERROR: jtot/=nu_jtot'
   WRITE(msg_unit,*)'jtot =',jtot
   WRITE(msg_unit,*)'nu_jtot=',nu_jtot
   STOP 'Read_Matrix2'
ENDIF

IF(proj/=nu_proj.AND.jthmat/=1)THEN
   WRITE(msg_unit,*)'ERROR: proj/=nu_proj'
   WRITE(msg_unit,*)'proj =',proj
   WRITE(msg_unit,*)'nu_proj=',nu_proj
   STOP 'Read_Matrix2'
ENDIF

IF(approx/=nu_approx.AND.jthmat/=1)THEN
   WRITE(msg_unit,*)'ERROR: approx/=nu_approx'
   WRITE(msg_unit,*)'approx =',approx
   WRITE(msg_unit,*)'nu_approx=',nu_approx
   STOP 'Read_Matrix2'
ENDIF

IF(symm/=nu_symm.AND.jthmat/=1)THEN
   WRITE(msg_unit,*)'ERROR: symm/=nu_symm'
   WRITE(msg_unit,*)'symm =',symm
   WRITE(msg_unit,*)'nu_symm=',nu_symm
   STOP 'Read_Matrix2'
ENDIF
ithmat = jthmat
proj = nu_proj
approx = nu_approx
symm = nu_symm

CALL ReadMat (amat, energies, chanl, elect, nvib, jrot, lorb, n, inunit)

RETURN
!
! This is reached at the end of the input file.
!
6  more = .false.
RETURN
ENDSUBROUTINE Read_Matrix2
