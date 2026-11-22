SUBROUTINE APHlogder
!      main program for log derivative propagation
!
USE FileUnits_MODULE
USE region_MODULE
USE InputFile_Module
USE pbasis_MODULE
USE Energy_Module
IMPLICIT NONE
LOGICAL there
INTEGER kenergy, jenergy, test, lnblnk, ncall, mcall, IErr
!       INTEGER nlam, lammin, lammax, megamin, megamax, coupmega
INTEGER jref, parity, mega, nrho
!
INTEGER jtot
LOGICAL symmetry,jeven
!
CHARACTER(LEN=5) ename
CHARACTER(LEN=5) blank
CHARACTER(LEN=110) R_file, RE_File
CHARACTER(LEN=125) todo, working, done, evalue, dworking, ddone
CHARACTER(LEN=200) pmatfile, ovrfile, sbfrst, jtotfile, corfile, asymfile
CHARACTER(LEN=256) filename
REAL(Kind=WP_Kind) cpu1, cpu2, smx_time,sge_time
!#include <Inc.region>
!#include <Inc.unitno>
NAMELIST/pmaovr/ pmatfile, ovrfile, sbfrst, jtotfile, corfile, asymfile
!       COMMON / basis / nlam(0:1000), lammin, lammax, megamin, megamax,
!     >                  coupmega
DATA blank /'     '/
DATA Jref/0/

WRITE(Out_Unit,*)'Called    APHLogder'
CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/APHLogder.txt',Form='FORMATTED')
iregion='aph    '
!      WRITE(Msg_Unit,*)'megamin1=',megamin,lammin,lammax,megamax,coupmega
!      megamin=0
!      STOP 'here'
! --------------------------------------------------------------------
! If we are in the aph region, we don't need the rho switch variable
! since we are propagating using only the logder method.
! --------------------------------------------------------------------

!      CALL FileUnits
!  -------------------------------------------------------------------
!  name and open logder.machinename.out file on unit Out_Unit
!  -------------------------------------------------------------------
!temp      test = get_hostnm(evalue)
test=0
evalue='Parker-APH_Dell_Precision_T7910'
filename='Output/logder.'//evalue(:lnblnk(evalue))//'.txt'
WRITE(Msg_Unit,*)'filename=',filename,':'
!OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//filename, status='unknown')
WRITE(Out_Unit,*) ' log-derivative propagation '
!      WRITE(Msg_Unit,*)'megamin2=',megamin
!---------------------------------------------------------------------
!  read data
!---------------------------------------------------------------------
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//TRIM(InputFile),status='old')
   READ(In_Unit,  NML=TotEnergy, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'TotEnergy')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist TotEnergy'
      WRITE(Msg_Unit,*)'ERROR with Namelist TotEnergy'
      STOP 'APHLogder: TotEnergy'
   ENDIF
   IF(eV_Input)THEN
      Efirst=Efirst_eV/autoeV
      DeltaEng=DeltaEng_eV/autoeV
   ENDIF
CLOSE(Unit=In_Unit)
WRITE(Out_Unit,NML=TotEnergy)
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
READ(In_Unit,NML=pmaovr, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'PmaOvr')
      pmatfile='BinOut/PotMatrx.bin'
      ovrfile='BinOut/Ovrlp.bin'
      sbfrst='BinOut/SetBasis_First.bin'
      jtotfile='BinOut/Diag_Jtot.bin'
      corfile='BinOut/Coriolis.bin'
      asymfile='BinOut/Asym_Top.bin'
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist PmaOvr'
      WRITE(Msg_Unit,*)'ERROR with Namelist PmaOvr'
      STOP 'APHLogder: PmaOvr'
   ENDIF
CLOSE(Unit=In_Unit)
WRITE(Out_Unit,NML=pmaovr)
WRITE(Out_Unit,*)'pmatfile=',pmatfile(:lnblnk(pmatfile))
WRITE(Out_Unit,*)'ovrfile=',ovrfile(:lnblnk(ovrfile))
WRITE(Out_Unit,*)'sbfrst=',sbfrst(:lnblnk(sbfrst))
WRITE(Out_Unit,*)'jtotfile=',jtotfile(:lnblnk(jtotfile))
WRITE(Out_Unit,*)'corfile=',corfile(:lnblnk(corfile))
WRITE(Out_Unit,*)'asymfile=',asymfile(:lnblnk(asymfile))
OPEN(Unit=pmat_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//pmatfile,form='unformatted')
OPEN(Unit=jtot_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//jtotfile,form='unformatted')
OPEN(Unit=cor_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//corfile,form='unformatted')
OPEN(Unit=asym_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//asymfile,form='unformatted')
OPEN(Unit=sbfrst_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//sbfrst,form='unformatted')
REWIND(sbfrst_unit)
! revised by X.L.
!      READ(sbfrst_unit)itemp1,itemp2,itemp3,itemp4
!      WRITE(Msg_Unit,*)'line1',itemp1,itemp2,itemp3,itemp4
!      READ(sbfrst_unit)itemp1,itemp2,itemp3,itemp4,itemp5,itemp6,itemp7
!      WRITE(Msg_Unit,*)'line2',itemp1,itemp2,itemp3,itemp4,itemp5,itemp6,itemp7
READ(sbfrst_unit)
READ(sbfrst_unit)
READ(sbfrst_unit)
READ(sbfrst_unit)
!      READ(sbfrst_unit)jtot,parity,mega,itemp1

!      jtot, parity, mega, megamin, megamax,  symmetry, jeven
!      READ(sbfrst_unit) jref, parity, mega, megamin, megamax, coupmega
!      READ(sbfrst_unit) jref, parity, mega, megamin, megamax,symmtemp,
!     > jeventemp
READ(SBFrst_Unit)jtot,parity,mega,megamin,maxmega,symmetry,jeven
coupmega=jref
!      READ(sbfrst_unit)
READ(sbfrst_unit) nrho
WRITE(Out_Unit,*) jref, parity, mega, megamin, maxmega, coupmega
WRITE(Out_Unit,*) nrho
CLOSE(Unit=sbfrst_unit)
!----------------------------------------------------------------------
!  The following is for distributed processing.
!  runs through the energies, finds which energy next needs to be done,
!  and puts an empty file, Rvwxyz, where vwxyz=ename=jenergy, into
!  the working subdirectory.
!----------------------------------------------------------------------
ename =blank

DO 2 kenergy = 1, nenergy
    IF(kenergy<10)THEN
        WRITE(ename,'(i1)') kenergy
    ELSEIF(kenergy<100)THEN
        WRITE(ename,'(i2)') kenergy
    ELSEIF(kenergy<1000)THEN
        WRITE(ename,'(i3)') kenergy
    ELSEIF(kenergy<10000)THEN
        WRITE(ename,'(i4)') kenergy
    ELSEIF(kenergy<100000)THEN
        WRITE(ename,'(i5)') kenergy
    ELSE
        WRITE(Msg_Unit,*)'kenergy to large=',kenergy
        STOP 'main'
    ENDIF

    R_file = 'R'//ename(1:LEN(TRIM(ename)))
    RE_file = 'RE'//ename(1:LEN(TRIM(ename)))
    ename=blank
    
    todo=   'todo/'//R_file(1:LEN(TRIM(R_File)))
    working='working/'//R_file(1:LEN(TRIM(R_File)))
    done=   'done/'//R_file(1:LEN(TRIM(R_File)))
    
    dworking = 'working/'//RE_file(1:LEN(TRIM(RE_File)))
    ddone =    'done/'//R_file(1:LEN(TRIM(RE_File)))

    INQUIRE(file=OutDIR(1:LEN(TRIM(OutDIR)))//working,exist=there)
!          IF(there) GOTO 2
    !IF(there) WRITE(Msg_Unit,*)' Warning, you are not checking', ' the working dir'

    INQUIRE(file=OutDIR(1:LEN(TRIM(OutDIR)))//done,exist=there)
    IF(there) GOTO 2

!temp          INQUIRE(file=OutDIR(1:LEN(TRIM(OutDIR)))//todo,exist=there)
!temp          IF(.NOT.there)THEN
!temp             OPEN(Unit=todo_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//todo,status='new',
!temp     >            form='unformatted')
!temp             CLOSE(Unit=todo_unit)
!temp          ENDIF

!temp          INQUIRE(file=OutDIR(1:LEN(TRIM(OutDIR)))//todo,exist=there)
!temp          IF(there)THEN
        jenergy = kenergy
!temp            test = rename(todo, working)
!temp             IF(test==0) GOTO 3
!temp          ENDIF

!temp    2  CONTINUE

!temp       WRITE(Msg_Unit,*)'All propagations are now complete'
!temp       STOP 'complete'

3 WRITE(Msg_Unit,*)
WRITE(Msg_Unit,'("Starting jenergy=",I5)')jenergy
CALL cputime(cpu1)
REWIND pmat_unit
REWIND jtot_unit
REWIND cor_unit
REWIND asym_unit
!---------------------------------------------------------------------
!  CALL propagator
!---------------------------------------------------------------------
!
ncall = 0
smx_time = 0.d0
mcall = 0
sge_time = 0.d0
CALL logprop(jenergy, working, dworking, ncall, smx_time, mcall, sge_time, ovrfile, jref, parity, nrho)
!
!---------------------------------------------------------------------
!  propagation done
!---------------------------------------------------------------------
CALL cputime(cpu2)
WRITE(Msg_Unit,34)jenergy,cpu2-cpu1
WRITE(Out_Unit,34)jenergy,cpu2-cpu1
WRITE(Out_Unit,'(A,I5)')'No. of calls to smxinv=',ncall
WRITE(Out_Unit,'(A,f12.5)')'Cpu time in smxinv=',smx_time
WRITE(Out_Unit,'(A,I5)')'No. of calls to dgemm=',mcall
WRITE(Out_Unit,'(A,f12.5)')'Cpu time in dgemm=',sge_time
!--------------------------------------------------------------------
!  move the completed R matrix from working to done subdirectory.
!--------------------------------------------------------------------
!temp       test = rename(working, done)
IF(test==0)THEN
    WRITE(Out_Unit,*)'Renamed R_file'
ELSE
    WRITE(Msg_Unit,*)'Error: Cannot rename file=',working
ENDIF

INQUIRE(file=OutDIR(1:LEN(TRIM(OutDIR)))//dworking,exist=there)
IF(there)THEN
!temp          test = rename(dworking, ddone)
    IF(test==0)THEN
        WRITE(Out_Unit,*)'Renamed deriv R_file'
    ELSE
        WRITE(Msg_Unit,*)'Error: Cannot rename file=',dworking
    ENDIF
ENDIF
2 CONTINUE
WRITE(Msg_Unit,*)
WRITE(Msg_Unit,*)'All propagations are now complete'
WRITE(Out_Unit,*)'All propagations are now complete'

CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/BasisForScatt.txt',Form='FORMATTED',ACCESS='APPEND')
WRITE(Out_Unit,*)'Completed APHLogder'
WRITE(Out_Unit,*)

RETURN
!---------------------------------------------------------------------
!  go back to start another energy.
!---------------------------------------------------------------------
34  FORMAT('Completed jenergy=',i5,5x,'CPU time=',f12.5)
ENDSUBROUTINE APHLogder
