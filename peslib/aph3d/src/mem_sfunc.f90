SUBROUTINE mem_sfunc( nel, numnp, ndisce, nq, nidmjunk, neq, nwk, nnc, nwork, nlancz, kase, pmax)
USE point_sfunc_Module
USE FileUnits_Module
USE Nidm_Module
USE Point_RSF_Module
CHARACTER(LEN=*) kase
LOGICAL little, medium, full
INTEGER ithcall, ithsub, pmax, nwk, nel, ndisce, neq, nlancz, nnc
INTEGER numnp, nwork, maxw, maxiw, nidmjunk, nq

!-----------------------------------------------------------------------
!  Determine printing options.
!-----------------------------------------------------------------------
DATA ithcall/0/, ithsub/0/
DATA little/.false./, medium/.false./, full/.false./
CALL popt ('mem_sfunc', little, medium, full, ithcall, ithsub)

IF(kase=='ALLOCATE'.or.kase=='release')THEN
   IF(nwk==0.or.kase=='release')THEN
      IF(little)THEN
         WRITE(Out_unit,*)'Starting to ',kase, ' memory in mem_sfunc for nwk=0'
         WRITE(Out_unit,*)'      name    length     start    finish'
      ENDIF
      IF(kase=='ALLOCATE')THEN
         ALLOCATE(lm(nel*18))
         WRITE(Out_Unit,*)'Is idi allocated?:',ALLOCATED(idi)
         ALLOCATE(idi(ndisce,5))
         ALLOCATE(maxa(neq+1))
         ALLOCATE(mht(neq))
         ALLOCATE(nloc(neq))
         ALLOCATE(nsit(nq))
         ALLOCATE(yz(nel*18))
      ELSE
         DEALLOCATE(lm)
         DEALLOCATE(idi)
         DEALLOCATE(maxa)
         DEALLOCATE(mht)
         DEALLOCATE(nloc)
         DEALLOCATE(nsit)
         DEALLOCATE(yz)
      ENDIF
   ENDIF
   IF(nwk==0)ALLOCATE(bghaml(1000000))
   IF(nwk/=0)THEN
         IF(little)THEN
            WRITE(Out_unit,*)'Starting to ',kase, ' memory in mem_sfunc for nwk>0'
         ENDIF
      IF(nlancz==0)THEN
         IF(kase=='ALLOCATE')THEN
            ALLOCATE(ar(nnc))
            ALLOCATE(br(nnc))
            ALLOCATE(vec(nq,nq))
            ALLOCATE(eigv(nq))
            ALLOCATE(d(nq))
            ALLOCATE(tt(neq))
            ALLOCATE(w(neq))
            ALLOCATE(nloc(neq))
            ALLOCATE(rtolv(nq))
            ALLOCATE(evc1(nq))
            ALLOCATE(evc2(numnp))
            ALLOCATE(ww(neq))
            ALLOCATE(nsit(nq))
            IF(ALLOCATED(bghaml))DEALLOCATE(bghaml)
            ALLOCATE(bghaml(2*nwk))
         ELSE
            DEALLOCATE(ar)
            DEALLOCATE(br)
            DEALLOCATE(vec)
            DEALLOCATE(eigv)
            DEALLOCATE(d)
            DEALLOCATE(tt)
            DEALLOCATE(w)
            DEALLOCATE(nloc)
            DEALLOCATE(rtolv)
            DEALLOCATE(evc1)
            DEALLOCATE(evc2)
            DEALLOCATE(ww)
            DEALLOCATE(nsit)
         ENDIF
      ELSEIF(nlancz==1)THEN
         IF(kase=='ALLOCATE')THEN
            ALLOCATE(ar(nq*4))
            ALLOCATE(d(nwork))
            ALLOCATE(tt(neq))
            !ALLOCATE(nloc(nq))
            ALLOCATE(evc2(numnp))
            IF(ALLOCATED(bghaml))DEALLOCATE(bghaml)
            ALLOCATE(bghaml(2*nwk))
         ELSE
            DEALLOCATE(ar)
            DEALLOCATE(d)
            DEALLOCATE(tt)
            DEALLOCATE(nloc)
            DEALLOCATE(evc2)
            DEALLOCATE(bghaml)
         ENDIF
      ELSEIF(nlancz==2)THEN
!  note.  the following must match what is in sfunc.
         IF(kase=='ALLOCATE')THEN
            ALLOCATE(ar(max(neq,(2*neq+6*pmax+4)/2)))
            IF(ALLOCATED(bghaml))DEALLOCATE(bghaml)
            ALLOCATE(bghaml(max(3*nwk,2*nwk+(2*pmax+5)*neq +pmax*(pmax+5))))
         ELSE
            DEALLOCATE(ar)
            DEALLOCATE(bghaml)
         ENDIF
      ELSE
         WRITE(Out_unit,*)'nlancz must be >=0 and <=2'
         STOP 'mem_sfunc'
      ENDIF
   ENDIF
   IF(little)THEN
      WRITE(Out_unit,*)'Memory ',kase,'d in mem_sfunc'
   ENDIF

ELSE
   WRITE(Out_unit,*)'Error kase /= ALLOCATE or release'
   STOP 'mem_sfunc'
ENDIF
maxw=max(3*nwk,2*nwk+(2*pmax+5)*neq+pmax*(pmax+5))
maxiw=max(2*neq,2*neq+6*pmax +4)
WRITE(Out_Unit,*)'maxw,maxiw,nwk,pmax,neq=',maxw,maxiw,nwk,pmax,neq

RETURN
END
