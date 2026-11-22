SUBROUTINE intbaf(Nang_DVR, ndvr, nabm, Phi_DVR, Phi_ABM, overlap, weight, th, ch, test, nbasis, neigmin, nbasiss)
USE Numeric_Kinds_Module
USE Numbers_Module
USE FileUnits_Module
IMPLICIT NONE
LOGICAL, PARAMETER:: Medium=.True.
CHARACTER(LEN=3) test
INTEGER Nang_DVR, ndvr, nabm, nbasis, neigmin, nbasiss, idvr, iabm, iang
REAL(Kind=WP_Kind) Phi_DVR(Nang_DVR,ndvr)
REAL(Kind=WP_Kind) overlap(ndvr,nabm), ovr 
REAL(Kind=WP_Kind) Phi_ABM(Nang_DVR,nbasiss)
REAL(Kind=WP_Kind) weight(Nang_DVR), th(Nang_DVR), ch(Nang_DVR)

WRITE(Out_Unit,*)'intbaf called:'
WRITE(Out_Unit,*)'Nang_DVR,ndvr,nabm=',Nang_DVR,ndvr,nabm
WRITE(Out_Unit,*)'nbasis,neigmin,nbasiss=',nbasis, neigmin, nbasiss

IF(test=='DVR'.or.test=='ABM')THEN
   IF(ndvr/=nabm)THEN
      WRITE(Out_Unit,*)'Error: ndvr should be equal to nabm in this case'
     STOP 'intbaf'
   ENDIF
ENDIF

CALL vsets (ndvr*nabm, overlap, 1, Zero)

!--------------------------------------------------------------------
! calculate ABM primitive basis functions at DVR points.
!--------------------------------------------------------------------
CALL bascal(th, ch, nbasis, neigmin, Phi_ABM, nbasiss, Nang_DVR)

IF(test=='ABM')THEN
   Phi_DVR=0.D0
   CALL vcopy (Nang_DVR*nabm, Phi_DVR, Phi_ABM)  ! Phi_DVR is now equal to Phi_ABM
   WRITE(Out_Unit,*)'Only using ABM Surface Functions'
ELSEIF(test=='DVR')THEN
   Phi_ABM=0.D0
   CALL vcopy (Nang_DVR*ndvr, Phi_ABM, Phi_DVR)  ! Phi_ABM is now equal to Phi_DVR
   WRITE(Out_Unit,*)'Only using DVR Surface functions'
ELSE
   WRITE(Out_Unit,*)'Using both DVR and ABM Surface functions'
ENDIF

IF(Medium)THEN
   WRITE(Out_Unit,'(/2x,"In IntBaf: Phi_DVR surface function"/)')
   CALL MxOutf(Phi_DVR, Nang_DVR, ndvr, Nang_DVR, ndvr)
   WRITE(Out_Unit,'(/2x,"In IntBaf: Phi_ABM surface function"/)')
   CALL MxOutf(Phi_ABM, Nang_DVR, nbasiss, Nang_DVR, nbasiss)
ENDIF

WRITE(Out_Unit,*)'Calculating Overlaps'
WRITE(Msg_Unit,*)'Calculating Overlaps'
!--------------------------------------------------------------------
! construct overlap matrix.
!--------------------------------------------------------------------
DO idvr = 1, ndvr
   DO iabm = 1, nabm
      ovr=0.d0
      DO iang = 1, Nang_DVR
         ovr = ovr + Phi_DVR(iang,idvr)*Phi_ABM(iang,iabm)*weight(iang)
      ENDDO
      overlap(idvr,iabm)=ovr
   ENDDO
ENDDO

IF(Medium)THEN
   WRITE(Out_Unit,'(/2x,"IntBaf: surface function overlap matrix elements"/)')
   CALL MxOutf(overlap, ndvr, nabm, ndvr, nabm)
ENDIF

WRITE(Out_Unit,*)'Finished Calculating Overlaps'

RETURN
ENDSUBROUTINE intbaf
