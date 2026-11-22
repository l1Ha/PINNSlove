SUBROUTINE mem_picture (nelemx, kase)
USE FileUnits_Module
USE point_picture_Module
CHARACTER(LEN=*) kase
LOGICAL little, medium, full
INTEGER ithcall, ithsub, nelemx

!-----------------------------------------------------------------------
!  Determine printing options.
!-----------------------------------------------------------------------
DATA ithcall/0/, ithsub/0/
DATA little/.false./, medium/.false./, full/.false./
CALL popt ('mem_picture', little, medium, full, ithcall, ithsub)

IF(little)THEN
WRITE(Out_unit,*)'Starting to ',kase,' memory in mem_picture'
WRITE(Out_unit,*)'      name    length     start    finish'
ENDIF

IF(kase=='ALLOCATE'.or.kase=='release')THEN
   IF(nelemx>=1)THEN
      ALLOCATE(igloss(nelemx))
      ALLOCATE(kolr(nelemx))
      ALLOCATE(center(3*nelemx))
      ALLOCATE(rad(nelemx))

      IF(little)THEN
      WRITE(Out_unit,*)'Memory ',kase,'d in mem_picture'
      ENDIF

   ELSE
      WRITE(Out_unit,*)'Error nelemx must be greater than 1'
      STOP 'mem_picture'
   ENDIF
ELSE
   WRITE(Out_unit,*)'Error kase /= ALLOCATE or release'
   STOP 'mem_picture'
ENDIF

RETURN
END
