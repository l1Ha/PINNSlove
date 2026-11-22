SUBROUTINE Read_Overlap(distance, overlap, n)
!-------------------------------------------------------------------
! SUBROUTINE: read_overlap
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
! purpose   :  Read overlap matrices
!
! input     :    distance -> rho  value for which the overlap matrix is needed
!                n        -> size of the overlap matrix
!
! output    :    overlap  -> overlap matrix
!
! modules   : time_Module -> timing
!             FileNames_Module-> file names
!             nrho_Module -> nrho value
!             jtot_Module -> Jtotal propagation value
! common    :
!
!
! notes     :
!
!-------------------------------------------------------------------
USE Time_Conv_Module
USE FileNames_Module
USE jtot_Module
USE nrho_Module
IMPLICIT NONE
REAL(Kind=WP_Kind) distance
REAL(Kind=WP_Kind) overlap(n,n)
INTEGER n

!
CHARACTER(LEN=11) :: Blank
CHARACTER(LEN=50), PARAMETER:: Overlap_File='Overlap.bin'
CHARACTER(LEN=21), PARAMETER::  procname='read_overlap'
LOGICAL, save            ::  debug
LOGICAL firstcall,isthere,opn
data firstcall /.true./
SAVE firstcall
!
REAL(Kind=WP_Kind) rho2,rholast
INTEGER nsfunc2,nrho2
INTEGER Jtot2
CHARACTER(LEN=3) Parity2,symm2
LOGICAL fGP2
REAL(Kind=WP_Kind) phi02,theta02
INTEGER mvec2,i,j
!
! Initialize
!
overlap = 0d0
!
! Open file in first CALL if not already opened
!
IF(firstcall)THEN
   CALL get_debug(procname,debug)
   INQUIRE(file=OutDIR(1:LEN(TRIM(OutDIR)))//overlap_file,EXIST=isthere,OPENED=opn)
   IF(isthere)THEN
      IF(.NOT.opn)THEN
            OPEN(Unit= overlap_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//overlap_file, status='old',form='unformatted')
      ENDIF
   ELSE
      WRITE(msg_unit,*) 'overlap_file = ',overlap_file
      CALL errormsg(msg_unit,procname,' Overlap file does not exist')
   ENDIF
   firstcall = .false.
ENDIF
!
! Read Header
!
READ(overlap_unit) rho2,rholast,nrho2,nsfunc2
READ(overlap_unit) Jtot2,Parity2,symm2
READ(overlap_unit) fGP2,phi02,theta02,mvec2
!
! Debug
!
IF(debug)THEN
   WRITE(Out_Unit,*) 'rho2 =',rho2,' rholast=', rholast
   WRITE(Out_Unit,*) 'nrho2=',nrho2,'nsfunc2=',nsfunc2
   WRITE(Out_Unit,*) 'Jtot2=',Jtot2
   WRITE(Out_Unit,*) 'Par. =',Parity2
   WRITE(Out_Unit,*) 'Symm.=',symm2
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

IF(ABS(rho2-distance) > 1.0E-9)THEN
   WRITE(msg_unit,*) 'rho2     = ',rho2
   WRITE(msg_unit,*) 'distance = ',distance
   CALL errormsg(msg_unit,procname, 'rho2 ne distance')
ENDIF

DO i=1,nsfunc2
   READ(overlap_unit) (overlap(i,j),j=1,nsfunc2)
ENDDO
!
! Debug
!
IF(debug)THEN
   CALL MatrixOut(overlap,nsfunc2,nsfunc2,'Overlap-read', 'overlap',Blank, Blank,.false., Blank, .false.)
ENDIF
!
! Complete the overlap
!
IF(nsfunc2<n)THEN
   DO i=nsfunc2+1,n
      DO j=nsfunc2+1,n
            overlap(i,j) = 0d0
            IF(i==j) overlap(i,j) = 1d0
      ENDDO
   ENDDO
ENDIF

!
! Close file for last CALL
!
IF(nrho2==nrho) CLOSE(overlap_unit)

RETURN
ENDSUBROUTINE Read_Overlap
