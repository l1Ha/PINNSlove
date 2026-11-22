SUBROUTINE ReadDelBas_Unit(nt, xsq, xk, eint, nchanacc)
USE FileUnits_Module
USE Narran_Module
USE Arrch_Module
USE Masses_Module
USE NumNuj_Module
USE Oops_Module
USE Qall_Module
USE Storage_Module
USE VFunc_Module
USE TotalEng_Module
USE Region_Module
USE GaussB_Module, ONLY : nhermt
IMPLICIT NONE
INTEGER I, Index, Nchanacc, NT
REAL(Kind=WP_Kind) eshift
REAL(Kind=WP_Kind) xsq(nt+nchanacc), xk(nt+nchanacc), eint(nt+nchanacc)

IF(nt+nchanacc>nvibrot.or.nt<=0.or.nvibrot<=0)THEN
   WRITE(Msg_Unit,*)"Error in IO_Basis nt,nchanacc,nvibrot:",nt, nchanacc, nvibrot
   WRITE(Msg_Unit,*)"Error in IO_Basis nt,nchanacc,nvibrot:",nt, nchanacc, nvibrot
   STOP "ERROR: Stopping in IO_Basis"
ENDIF


xktot=usys2*Etot
READ(bas_unit)nt,nnuj(karran)
WRITE(*,*)nt,nnuj(karran)
READ(bas_unit)xkold,ioops
WRITE(*,*)xkold,ioops
READ(bas_unit)(mvib3(i+nchanacc),jrot3(i+nchanacc),lorb3(i+nchanacc),   &
               mega3(i+nchanacc),kchan(i+nchanacc),xksq3(i+nchanacc),   &
               nuj3(i+nchanacc), i=1,nt)
WRITE(*,*)(mvib3(i+nchanacc),jrot3(i+nchanacc),lorb3(i+nchanacc),   &
               mega3(i+nchanacc),kchan(i+nchanacc),xksq3(i+nchanacc),   &
               nuj3(i+nchanacc), i=1,nt)
                  
IF(karran==1)THEN
  Index=0
ELSEIF(karran==2)THEN
  Index=nnuj(1)
ELSEIF(karran==3)THEN
  Index=nnuj(1)+nnuj(2)
ENDIF
   
READ(bas_unit)(chinuj(Index+i), i=1,nnuj(karran))   
WRITE(*,*)(chinuj(Index+i), i=1,nnuj(karran))
READ(bas_unit)(xsq(i),xk(i),eint(i), i=1,nt)
WRITE(*,*)(xsq(i),xk(i),eint(i), i=1,nt)
eshift=xktot-xkold
DO i=1,nt
   xsq(i)=xsq(i)+eshift
   xk(i)=sqrt(ABS(xsq(i)))
   xksq3(i+nchanacc)=xksq3(i+nchanacc)+eshift
ENDDO
WRITE(*,*)(xksq3(i+nchanacc),xsq(i),xk(i),eint(i), i=1,nt)

RETURN
ENDSUBROUTINE ReadDelBas_Unit