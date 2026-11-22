!-------------------------------------------------------------------
! SUBROUTINE: read_potmat
!
! package   : CID
!
! Language  : Fortran 90
!
! author    : F. Colavecchia (flavioc@lanl.gov)
!
! date      : Apr-02       version: 0.1
! revision  :              version:
!
! purpose   :  Read potmat matrices
!
! input     :    distance -> rho  value for which the potmat matrix is needed
!                n        -> size of the potmat matrix
!
! output    :    pm1  -> upper triangular potmat matrix at sector_left  for DBS propagator
!                pm2  -> upper triangular potmat matrix at sector_right for DBS propagator
!                or
!                pm1  -> eigenvectors at distance for SVD propagator
!                pm2  -> dummy
!                eigen-> eigenenergies of the basis.
!
! modules   : time_Module -> timing
!             FileNames_Module-> file names
!             boundary_Module -> last sector rho value
!             jtot_Module -> Jtotal propagation value
! common    :
!
!
! notes     :
!
!-------------------------------------------------------------------
SUBROUTINE read_potmat(distance,pm1,pm2,eigen, n)
USE Time_Conv_Module
USE FileNames_Module
USE jtot_Module
USE boundary_Module
USE BasisRegions_Module
USE method_Module
IMPLICIT NONE
REAL(Kind=WP_Kind) distance
REAL(Kind=WP_Kind) pm1(n,n),pm2(n,n)
REAL(Kind=WP_Kind) eigen(n)
INTEGER n

!
CHARACTER(LEN=11) :: Blank
CHARACTER(LEN=21), PARAMETER ::  procname='read_potmat'
LOGICAL, save            ::  debug
LOGICAL firstcall,isthere
data firstcall /.true./
SAVE firstcall
!
REAL(Kind=WP_Kind) rho2,rholast
REAL(Kind=WP_Kind) rhos(3)
INTEGER nsfunc2,nrho2
INTEGER Jtot2
CHARACTER(LEN=3) Parity2,symm2
LOGICAL fGP2
REAL(Kind=WP_Kind) phi02,theta02
INTEGER mvec2,i,j
!
! Initialize
!
pm1 = 0d0
pm2 = 0d0
!
! Open file in first CALL
!
IF(firstcall)THEN
    CALL get_debug(procname,debug)
    IF(bmethod=='svd')THEN
        INQUIRE(file=OutDIR(1:LEN(TRIM(OutDIR)))//potmatsvd_name,EXIST=isthere)
        IF(isthere)THEN
            OPEN(Unit= potmatsvd_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//potmatsvd_name, status='old',form='unformatted')
        ELSE
            CALL errormsg(msg_unit,procname,' svd file does not exist')
        ENDIF
!            potmat_unit=potmatsvd_unit
    ELSEIF(bmethod=='dbs')  then
        INQUIRE(file=OutDIR(1:LEN(TRIM(OutDIR)))//potmatdbs_name,EXIST=isthere)
        IF(isthere)THEN
            OPEN(Unit= potmatdbs_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//potmatdbs_name, status='old',form='unformatted')
        ELSE
            CALL errormsg(msg_unit,procname,' dbs file does not exist')
        ENDIF
!            potmat_unit=potmatdbs_unit
    ELSE
        CALL errormsg(msg_unit,procname,'method not implemented')
    ENDIF
    firstcall = .false.
ENDIF
!
! Read Header
!
READ(potmat_unit) rho2,rhos,nsfunc2
READ(potmat_unit) Jtot2,Parity2,symm2
READ(potmat_unit) fGP2,phi02,theta02,mvec2
IF(debug)THEN
    WRITE(Out_Unit,*) 'bmethod = ',bmethod,' coord =', coord
    WRITE(Out_Unit,*) 'pm1 matrix at',rho2
    WRITE(Out_Unit,*) 'rho2   = ',rho2,   '  rhos(1)= ',rhos(1)
    WRITE(Out_Unit,*) 'rhos(2)= ',rhos(2),'  rhos(3)= ',rhos(3)
    WRITE(Out_Unit,*) 'Jtot   = ',jtot2
    WRITE(Out_Unit,*) 'parity = ',parity2,' symmetry= ',symm2
    IF(fGP2)THEN
        WRITE(Out_Unit,*) 'phi02  =',phi02
        WRITE(Out_Unit,*) 'theta02=',theta02
        WRITE(Out_Unit,*) 'mvec   =',mvec2
    ENDIF
ENDIF
!
! Check header
!
IF(nsfunc2 > n)THEN
    WRITE(msg_unit,*) 'nsfunc2= ',nsfunc2
    WRITE(msg_unit,*) 'n      = ',n
    CALL errormsg(msg_unit,procname, 'nsfunc2 gt n')
ENDIF

IF(Jtot2 /= Jtot)THEN
    WRITE(msg_unit,*) 'Jtot2 = ',Jtot2
    WRITE(msg_unit,*) 'Jtot  = ',Jtot
    CALL errormsg(msg_unit,procname, 'Jtot2 ne Jtot')
ENDIF

IF(ABS(rhos(2)-distance) > 1.0E-9)THEN
    WRITE(msg_unit,*) 'rhos(2)  = ',rhos(2)
    WRITE(msg_unit,*) 'distance = ',distance
    CALL errormsg(msg_unit,procname, 'rhos(2) ne distance')
ENDIF

READ(potmat_unit) (eigen(j),j=1,nsfunc2)

IF(bmethod=='dbs')THEN
    DO i=1,nsfunc2
        READ(potmat_unit) (pm1(i,j),j=i,nsfunc2)
        READ(potmat_unit) (pm2(i,j),j=i,nsfunc2)
    ENDDO
ELSE
    !
    !   APH basis is adiabatic, so we do not need to
    !   read the eigenvectors.
    !
    IF(coord=='del')THEN
        DO i=1,nsfunc2
            READ(potmat_unit) (pm1(i,j),j=1,nsfunc2)
        ENDDO
    ELSE
        CALL unit(nsfunc2,pm1)
    ENDIF
ENDIF
!
!   Fill in the lower triangle of the matrices
!
IF(bmethod=='dbs')THEN
    DO i=1,nsfunc2
        DO j=1,i
            pm1(i,j) = pm1(j,i)
            pm2(i,j) = pm2(j,i)
        ENDDO
    ENDDO
ENDIF

IF(debug)THEN
    CALL MatrixOut(pm1,nsfunc2,nsfunc2,'Pm1_Matrix','Pm1 matrix as read',Blank, Blank,.False.,Blank,.False.)
    IF(bmethod=='dbs')THEN
        CALL MatrixOut(pm2,nsfunc2,nsfunc2,'Pm2_Matrix','Pm2 matrix as read',Blank, Blank,.False.,Blank,.False.)
    ENDIF
    WRITE(Out_Unit,*) 'eigenvalues =',eigen
ENDIF

!
! Close file for last CALL
!
IF(rhos(2)==basis_dist(NSectors))THEN
    IF(bmethod=='dbs') CLOSE(potmatdbs_unit)
    IF(bmethod=='svd') CLOSE(potmatsvd_unit)
ENDIF

RETURN
END
