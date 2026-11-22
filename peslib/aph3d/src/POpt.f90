SUBROUTINE POpt(subr, little, medium, full, ithcall, ithsub)
!
!  P U R P O S E   O F    S U B R O U T I N E
!----------------------------------------------------------------------
!  this routine sets the following parameters for printing options:
!  little......minimal output option
!  medium......medium output option
!  full........debug output option
!
!  on entering the following should be defined:
!  subr........routine name (8 Hollerith CHARACTERs at most)
!  ithcall.....ith-CALL from routine subr (initialize to zero)


!  I N P U T    A R G U M E N T S
!  subr
!  little
!  medium
!  full
!  ithcall
!  ithsub
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!
USE Numeric_Kinds_Module
USE FileUnits_Module
USE Prntng_Module, ONLY: Option, IPrtAll, Iter, Iprt, nsubs
USE SubList_Module
IMPLICIT NONE
LOGICAL Little, medium, Full
INTEGER i, ithcall, prtopt, Ithsub
CHARACTER(LEN=*) subr
CHARACTER(LEN=21) subq
CHARACTER(LEN=21) name1, name2, blank

blank = '          '
!----------------------------------------------------------------------
!  check to see if parameters have already been set
!  and will not be altered (ithcall<0).
!----------------------------------------------------------------------
IF(option == 'all     ')THEN
   IF(iprtall == 0)THEN
      little = .false.
      medium = .false.
      full   = .false.
   ELSEIF (iprtall == 1)THEN
      little = .true.
      medium = .false.
      full   = .false.
   ELSEIF (iprtall == 2)THEN
      little = .true.
      medium = .true.
      full   = .false.
   ELSEIF (iprtall == 3)THEN
      little = .true.
      medium = .true.
      full   = .true.
   ELSE
      WRITE(Out_Unit,*)'Invalid iprtall'
      STOP 'popt'
   ENDIF
   RETURN
ENDIF
IF(ithcall < 0) RETURN
IF(ithcall > 0)THEN
   ithcall=ithcall+1
   IF(Ithsub>0)THEN
      IF(iter(3,ithsub)+iter(2,ithsub)+iter(1,ithsub) < ithcall)little=.false.
      IF(iter(3,ithsub)+iter(2,ithsub) < ithcall) medium=.false.
      IF(iter(3,ithsub) < ithcall) full=.false.
      IF( .NOT. little) ithcall=-1
   ENDIF
   WRITE(Out_Unit,13) subr, ithcall, little, medium, full
   RETURN
ENDIF
!----------------------------------------------------------------------
!  set default printing options.
!----------------------------------------------------------------------
prtopt=0
little=.false.
medium=.false.
full=.false.
ithcall=-1
!----------------------------------------------------------------------
!  find routine `subr' in the list of routines and determine print
!  options.
!----------------------------------------------------------------------
return
DO i=1, nsubs
   name1 = blank
   name2 = blank
   subq = subs(i)
   name1 = TRIM(subr)   !subr(:LEN(subr)) // char(0) !GregParker TEMPMOD
   name2 = TRIM(subq)   !subq(:LEN(subq)) // char(0) !GregParker TEMPMOD
   IF(name1 == name2)THEN
      prtopt=iprt(i)
      IF(prtopt/=0)THEN
         WRITE(Out_Unit,'("Turning on print for routine: ",A, "Printing level = ", i2)')name1,prtopt
      ENDIF
      ithsub=i
      EXIT
   ENDIF
ENDDO
!----------------------------------------------------------------------
!  reset print options if needed.
!----------------------------------------------------------------------
IF(prtopt == 0)THEN       !  return if no printing is desired.
   RETURN
ELSEIF (prtopt == 1)THEN  !  minimal printing option.
   little=.true.
ELSEIF (prtopt == 2)THEN  !  medium printing option.
   little=.true.
   medium=.true.
ELSEIF (prtopt == 3)THEN  !  full debug printing option.
   little=.true.
   medium=.true.
   full=.true.
ELSEIF (prtopt == 4)THEN  !  printing level based on the number of calls made.
   IF(iter(3,ithsub)+iter(2,ithsub)+iter(1,ithsub) > 0)little=.true.
   IF(iter(3,ithsub)+iter(2,ithsub) > 0) medium=.true.
   IF(iter(3,ithsub) > 0) full=.true.
   ithcall=1
   WRITE(Out_Unit,13) subr, ithcall, little, medium, full
ELSE                      !  prtopt does not have correct value STOP!!
   WRITE(Out_Unit,*) 'error in routine popt ','(prtopt incorrect) ', prtopt
   STOP 'popt'
ENDIF

RETURN
 13   FORMAT(1x, "No. of times ", a8, " has been called=", i5, "  little=", l1, "  medium=", l1, "  full=", l1)
ENDSUBROUTINE POpt
