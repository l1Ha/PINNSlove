FUNCTION IGlobal (Command, Name, Value)
!
!-----------------------------------------------------------------------
!
! >>> IGlobal   -- INTEGER global variable saver routine
!
! *** Calling Sequence:         Result  = IGlobal(Command,Name,Value)
!
! --- on entry the calling routine supplies --
!
!     Command   -- CHARACTER string indicating whether the indicated glo
!                       variable is to be saved, retrieved, or deleted
!                       Possible Commands are:
!                               Command = 'save'        to save a global
!                               Command = 'get'         to retrieve one
!                               Command = 'delete'      to remove one
!                               Command = 'list'        to WRITE EVERYTH
!               NOTE:  only the first CHARACTER of Command is significan
!
!     Name      -- CHARACTER string (8 or fewer CHARACTERs) indicating t
!                       unique name of the global to be saved/retrieved/
!
!     Value     -- INTEGER which is an input variable if Command='save'
!                       If Command /= 'save', Value is not read on inp
!
! --- on exit the routine returns --
!
!     Value     -- is an output INTEGER variable only if Command='get',
!                       which case the remembered value of the Name-d
!                       global variable is returned.
!
!     IGlobal   -- is a CHARACTER string indicating the success/failure
!                       the operation.  Possible returns are:
!                               IGlobal = 'OK'    if all was OK
!                               IGlobal = 'CANT'  if no internal space l
!                                               a 'save' operation, or i
!                                               Named variable is not he
!                                               a 'delete' or 'get'.
!                               IGlobal = 'WHAT'  if Command unrecognize
!
!-----------------------------------------------------------------------
!
!
CHARACTER(LEN=4) :: IGlobal
CHARACTER(LEN=*) :: Command
CHARACTER(LEN=*) :: Name
INTEGER :: Value
!
INTEGER :: MAXSave
PARAMETER (MAXSave = 100)
!
INTEGER :: NSaved
! * * * * *                             !       loop index
INTEGER :: i
! * * * * *                             !       new value of NSaved
INTEGER :: new
! * * * * *                             !       value of output unit if
INTEGER :: iout
!
INTEGER :: Save (MAXSave)
CHARACTER(LEN=8) :: VNames (MAXSave)
!
SAVE
!
DATA NSaved / 5 /
!
DATA VNames (1) / 'StdIn' /
DATA Save (1) / 9 /
!
DATA VNames (2) / 'StdOut' /
DATA Save (2) / 9 /
!
DATA VNames (3) / 'StdErr' /
DATA Save (3) / 9 /
!
DATA VNames (4) / 'StdDbg' /
DATA Save (4) / 9 /
!
DATA VNames (5) / 'StdWrt' /
DATA Save (5) / 9 /
!
!
!---------------------------- execution begins here --------------------
!
!
!                        --- blank Names NOT allowed ---
!
IF(Name==' ')THEN
   IGlobal = 'CANT'
   RETURN
ENDIF
!
!                    --- special case if Command = 'list' ---
!
IF(Command (1:1) =='l')THEN
   DO 100 i = 1, NSaved
      IF(Name==VNames (i) )THEN
         iout = Save (i)
         goto 101
      ENDIF
  100    ENDDO
   iout = 0
   IGlobal = 'CANT'
  101    CONTINUE
   IF(iout==0) RETURN
!
WRITE(iout,  * ) '   === INTEGER Globals ==='
   DO 200 i = 1, NSaved
      IF(VNames (i) /=' ')THEN
         WRITE(iout, '(5x,a8,'' = '',i10)') VNames (i) , Save (i)
      ENDIF
  200    ENDDO
   IGlobal = 'OK'
   RETURN
ENDIF
!
!                     --- has Name already been defined? ---
!
DO 1 i = 1, NSaved
   IF(Name==VNames (i) ) goto 2
    1 ENDDO
i = 0
    2 CONTINUE
!
!---------------------------- saving a value ---------------------------
!
IF(Command (1:1) =='s')THEN
!
   IF(i==0)THEN
      new = NSaved+1
      IF(new>MAXSave)THEN
!
!               --- arrays full, look for a blank entry ---
!
         DO 10 i = 1, NSaved
            IF(VNames (i) ==' ') goto 11
   10          ENDDO
         i = 0
   11          CONTINUE
         IF(i==0)THEN
            IGlobal = 'CANT'
            RETURN
         ELSE
            new = NSaved
         ENDIF
      ENDIF
!
!                --- save value for the newly defined global ---
!
      NSaved = new
      VNames (NSaved) = Name
      Save (NSaved) = Value
      IGlobal = 'OK'
!
!           --- save a new value for a previously defined global ---
!
   ELSE
      Save (i) = Value
      IGlobal = 'OK'
   ENDIF
!
!----------------------------- retrieving a value ----------------------
!
ELSEIF (Command (1:1) =='g')THEN
!
   IF(i==0)THEN
      IGlobal = 'CANT'
   ELSE
      IGlobal = 'OK'
      Value = Save (i)
   ENDIF
!
!----------------------------- deleting a global -----------------------
!
ELSEIF (Command (1:1) =='d')THEN
!
   IF(i==0)THEN
      IGlobal = 'CANT'
   ELSE
      IGlobal = 'OK'
      VNames (i) = ' '
   ENDIF
!
!---------------------------- unrecognized command ---------------------
!
ELSE
   IGlobal = 'WHAT'
ENDIF
!
RETURN
END FUNCTION IGlobal
