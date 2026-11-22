SUBROUTINE mem_msher(ntheta, nchi, nel, numnp, kase)
USE FileUnits_Module
USE point_msher_Module
USE Nidm_Module
IMPLICIT NONE
LOGICAL little, medium, full
INTEGER ntheta, nchi, nel, numnp, malloc, ithcall, ithsub
CHARACTER(LEN=*) kase
EXTERNAL malloc

!-----------------------------------------------------------------------
!  Determine printing options.
!-----------------------------------------------------------------------
DATA ithcall/0/, ithsub/0/
DATA little/.false./, medium/.false./, full/.false./
CALL popt ('mem_msher', little, medium, full, ithcall, ithsub)

WRITE(Out_Unit,*)'Nidm=',Nidm
IF(Nidm/=5)THEN
   !STOP 'NIDM/=5'
ENDIF
IF(kase=='rel2')THEN
   IF(little)THEN
   WRITE(Out_unit,*)'Starting to release memory in mem_msher'
   ENDIF
   DEALLOCATE(nid)
   DEALLOCATE(beta)
   DEALLOCATE(nodj)
   IF(little)THEN
   WRITE(Out_unit,*)'Memory released in mem_msher'
   ENDIF
   RETURN
ENDIF

IF(kase=='ALLOCATE'.or.kase=='release')THEN
   IF(little)THEN
   WRITE(Out_unit,*)'Starting to ',kase,' memory in mem_msher'
   WRITE(Out_unit,*)'      name    length     start    finish'
   ENDIF
   IF(nel>=1.AND.numnp>=9)THEN
      IF(kase=='ALLOCATE')THEN
         ALLOCATE(lnod(2*nel))
         ALLOCATE(scr1(2* numnp))
         ALLOCATE(scr2(2*numnp))
         IF(ntheta>=1.AND.nchi>=1)THEN
            IF(.NOT.allocated(iscr1))ALLOCATE(iscr1(ntheta*nchi))  !ALLOCATE(iscr1(ntheta*(nchi/2+1))) !GregParker TEMPMOD
         ELSE
            WRITE(Out_unit,*)'ntheta and nchi must be >= 1'
            STOP 'mem_msher'
         ENDIF
         IF(.NOT.allocated(iscr2))ALLOCATE(iscr2( nel))
         IF(.NOT.allocated(Virtualx))ALLOCATE(Virtualx(max(ntheta,9),max(nchi,nel) ))
      ELSE
         DEALLOCATE(lnod)
         DEALLOCATE(scr1)
         DEALLOCATE(scr2)
         IF(ntheta>=1.AND.nchi>=1)THEN
         IF(ALLOCATED(iscr1))WRITE(Out_Unit,*)'iscr1 is allocated'
!            DEALLOCATE(iscr1)
         ELSE
            WRITE(Out_unit,*)'ntheta and nchi must be >= 1'
            STOP 'mem_msher'
         ENDIF
         IF(ALLOCATED(iscr2))WRITE(Out_Unit,*)'iscr2 is allocated'
!         DEALLOCATE(iscr2)
         IF(ALLOCATED(Virtualx))WRITE(Out_Unit,*)'Virtualx is allocated'
!         DEALLOCATE(Virtualx)
      ENDIF
      IF(kase=='ALLOCATE')THEN
         ALLOCATE(nid(2*numnp/3))
         ALLOCATE(beta(nidm,(2*numnp/3+1)))
         ALLOCATE(nodj(max((nidm+1)*(2*numnp/3+1),9*nel)))
      ENDIF
      IF(kase=='ALLOCATE')THEN
         ALLOCATE(idteq(numnp))
         ALLOCATE(idiv(nel))
         ALLOCATE(lcor(2*nel))
         ALLOCATE(mapr(nel))
         ALLOCATE(idivo(2*nel))
         ALLOCATE(lcoro(2*nel))
         ALLOCATE(mapnod(numnp))
         ALLOCATE(sinlist(nel))
         ALLOCATE(philist(nel))
         ALLOCATE(iwork33(2*nel))
      ELSE
         DEALLOCATE(idteq)
         DEALLOCATE(idiv)
         DEALLOCATE(lcor)
         DEALLOCATE(mapr)
         DEALLOCATE(idivo)
         DEALLOCATE(lcoro)
         DEALLOCATE(mapnod)
         DEALLOCATE(sinlist)
         DEALLOCATE(philist)
         DEALLOCATE(iwork33)
      ENDIF
   ELSE
      WRITE(Out_unit,*)'nel must be >=1 and numnp must be >=9'
      STOP 'mem_msher'
   ENDIF
   IF(little)THEN
      WRITE(Out_unit,*)'Memory ',kase,'d in mem_msher'
   ENDIF
ELSE
   WRITE(Out_unit,*)'Error kase /= ALLOCATE or release'
   STOP 'mem_msher'
ENDIF

RETURN
END
