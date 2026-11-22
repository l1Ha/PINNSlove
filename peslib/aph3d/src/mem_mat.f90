SUBROUTINE mem_mat( numnp, nq, nquad, kase)
USE FileUnits_Module
USE point_mat_Module
CHARACTER(LEN=*) kase
INTEGER Nq, Numnp, Nquad
LOGICAL little, medium, full
INTEGER ithcall, ithsub

!-----------------------------------------------------------------------
!  Determine printing options.
!-----------------------------------------------------------------------
DATA ithcall/0/, ithsub/0/
DATA little/.false./, medium/.false./, full/.false./
CALL popt ('mem_mat ', little, medium, full, ithcall, ithsub)

IF(little)THEN
WRITE(Out_unit,*)'Starting to ',kase,' memory in mem_mat'
WRITE(Out_unit,*)'      name    length     start    finish'
ENDIF

IF(kase=='ALLOCATE'.or.kase=='release')THEN
   IF(nq>=1.AND.nquad>=3.AND.numnp>=9)THEN
      IF(kase=='ALLOCATE')THEN
         ALLOCATE(fbar( nq,nquad,nquad))
         ALLOCATE(dfbar(nq,nquad,nquad))
         ALLOCATE(hamil(nq,nq))
         ALLOCATE(ff(numnp,nq))
      ELSEIF(kase=='release')THEN
         DEALLOCATE(fbar)
         DEALLOCATE(dfbar)
         DEALLOCATE(hamil)
         DEALLOCATE(ff)
      ENDIF
   ELSE
      IF(nq<1) WRITE(Out_unit,*)'Error nq must be >= 1'
      IF(nquad<3) WRITE(Out_unit,*)'Error nquad must be >= 3'
      IF(numnp<9) WRITE(Out_unit,*)'Error numnp must be >= 9'
      STOP 'mem_mat'
   ENDIF

   IF(little)THEN
   WRITE(Out_unit,*)'Memory ',kase,'d in mem_mat'
   ENDIF
ELSE
   WRITE(Out_unit,*)'Error kase /= ALLOCATE or release'
   STOP 'mem_mat'
ENDIF

RETURN
END
