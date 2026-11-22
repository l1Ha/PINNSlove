SUBROUTINE mem_ovr(numnp, nq, nquad, nel1, numnp1, ntheta,nchi, kase)
USE FileUnits_Module
USE Point_Ovr_Module
CHARACTER(LEN=*) kase
LOGICAL little, medium, full
INTEGER Numnp, Nq, NQuad, Nel1, Numnp1, Ntheta, Nchi
INTEGER ithcall, ithsub

!-----------------------------------------------------------------------
!  Determine printing options.
!-----------------------------------------------------------------------
DATA ithcall/0/, ithsub/0/
DATA little/.false./, medium/.false./, full/.false./
CALL popt ('mem_ovr ', little, medium, full, ithcall, ithsub)

IF(little)THEN
WRITE(Out_unit,*)'Starting to ',kase,' memory in mem_ovr'
WRITE(Out_unit,*)'      name    length     start    finish'
ENDIF

IF(kase=='ALLOCATE'.or.kase=='release')THEN
   IF(kase=='ALLOCATE')THEN
   IF(nq>=1.AND.nquad>=3.AND.numnp>=9)THEN
      ALLOCATE(fb( nq,nquad,nquad))
      ALLOCATE(fbp(nq,nquad,nquad))
      ALLOCATE(sovrlap(nq,nq))
   ELSE
      IF(nq<1) WRITE(Out_unit,*)'Error nq must be >= 1'
      IF(nquad<3) WRITE(Out_unit,*)'Error nquad must be >= 3'
      IF(numnp<9) WRITE(Out_unit,*)'Error numnp must be >= 9'
      STOP 'mem_ovr'
   ENDIF

   IF(nel1>=1.AND.numnp1>=9)THEN
      ALLOCATE(idlast( 9,nel1))
      ALLOCATE(thlast( numnp1))
      ALLOCATE(chlast( numnp1))
      ALLOCATE(philast(max(numnp,numnp1)*nq))
   ELSE
      IF(nel1<1)WRITE(Out_unit,*)'Error nel1 must be >= 1'
      IF(numnp1<9)WRITE(Out_unit,*)'Error numnp1 must be >= 9'
      STOP 'mem_ovr'
   ENDIF

   IF(ntheta>=1.AND.nchi>=1)THEN
      ALLOCATE(f3(numnp))
      ALLOCATE(chi1(numnp))
      ALLOCATE(ThetAPH(numnp))
      ALLOCATE(igrid(2*ntheta*nchi))
   ELSE
      WRITE(Out_unit,*)'Error ntheta and nchi must be >= 1'
      STOP 'mem_ovr'
   ENDIF

   IF(little)THEN
   WRITE(Out_unit,*)'Memory ',kase,'d in mem_ovr'
   ENDIF

   ELSE
   IF(nq>=1.AND.nquad>=3.AND.numnp>=9)THEN
      DEALLOCATE(fb)
      DEALLOCATE(fbp)
      DEALLOCATE(sovrlap)
   ELSE
      IF(nq<1) WRITE(Out_unit,*)'Error nq must be >= 1'
      IF(nquad<3) WRITE(Out_unit,*)'Error nquad must be >= 3'
      IF(numnp<9) WRITE(Out_unit,*)'Error numnp must be >= 9'
      STOP 'mem_ovr'
   ENDIF

   IF(nel1>=1.AND.numnp1>=9)THEN
      DEALLOCATE(idlast)
      DEALLOCATE(thlast)
      DEALLOCATE(chlast)
      DEALLOCATE(philast)
   ELSE
      IF(nel1<1)WRITE(Out_unit,*)'Error nel1 must be >= 1'
      IF(numnp1<9)WRITE(Out_unit,*)'Error numnp1 must be >= 9'
      STOP 'mem_ovr'
   ENDIF

   IF(ntheta>=1.AND.nchi>=1)THEN
      DEALLOCATE(f3)
      DEALLOCATE(chi1)
      DEALLOCATE(ThetAPH)
      DEALLOCATE(igrid)
   ELSE
      WRITE(Out_unit,*)'Error ntheta and nchi must be >= 1'
      STOP 'mem_ovr'
   ENDIF
   ENDIF

   IF(little)THEN
   WRITE(Out_unit,*)'Memory ',kase,'d in mem_ovr'
   ENDIF
ELSE
   WRITE(Out_unit,*)'Error kase /= ALLOCATE or release'
   STOP 'mem_ovr'
ENDIF

RETURN
END
