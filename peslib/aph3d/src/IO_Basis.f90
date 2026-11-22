SUBROUTINE IO_Basis(nt, xsq, xk, eint, nchanacc)
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
USE Rhovalue_Module, ONLY: Rhoval
USE GaussB_Module, ONLY : nhermt
IMPLICIT NONE
INTEGER I, Index, Nchanacc, NT, nb
REAL(Kind=WP_Kind) eshift
REAL(Kind=WP_Kind) xsq(nt+nchanacc), xk(nt+nchanacc), eint(nt+nchanacc)

IF(nt+nchanacc>nvibrot.or.nt<=0.or.nvibrot<=0)THEN
   WRITE(Msg_Unit,*)"Error in IO_Basis nt,nchanacc,nvibrot:",nt, nchanacc, nvibrot
   WRITE(Msg_Unit,*)"Error in IO_Basis nt,nchanacc,nvibrot:",nt, nchanacc, nvibrot
   STOP "ERROR: Stopping in IO_Basis"
ENDIF
IF(irdindep)THEN
   xktot=usys2*Etot
   READ(bas_unit)nt, nnuj(karran), nchanacc, nb, nhermt(karran), ioops, iregion
   READ(bas_unit)xkold,ioops
   READ(bas_unit)(mvib3(i+nchanacc),jrot3(i+nchanacc),lorb3(i+nchanacc),   &
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
   READ(bas_unit)(xsq(i+nchanacc),xk(i+nchanacc),eint(i+nchanacc), i=1,nt)
   eshift=xktot-xkold
   
!IS something wrong here????? tmpmodgregparker FIX THIS ERROR
   DO i=nchanacc+1,nt+nchanacc
      !xsq(i)=xsq(i)   !+eshift
      !xk(i)=sqrt(ABS(xsq(i)))
      !xksq3(i)=xsq(i)   !xksq3(i+nchanacc)+eshift
   ENDDO 
   IF(karran==3)THEN
      !WRITE(*,'(A,6e23.15)')"eint5=",eint(1),eint(36),eint(71),Energy3_Jacobi(1),Energy3_Jacobi(36),Energy3_Jacobi(71)
      do i=1,nt+nchanacc  !tmpmod new loop added here?
         xsq(i)=usys2*(etot-eint(i))
         xk(i)=sqrt(ABS(xsq(i)))
         xksq3(i)=xsq(i)
      enddo
      !WRITE(*,'(A,6e23.15)')"xsq5 =",xsq(1) ,xsq(36) ,xsq(71) ,xksq3(1)         ,xksq3(36)         ,xksq3(71)
      !if(eint(1)/=eint(36).or.xksq3(1)/=xksq3(36).or.eint(1)/=eint(71).or.xksq3(1)/=xksq3(71))THEN
      !   stop "error in routine input5"
      !endif
      !if(Energy3_Jacobi(1)/=Energy3_Jacobi(36).or.xsq(1)/=xsq(36))THEN
      !   stop "error in routine input6"
      !endif
   ENDIF
   !WRITE(*,*)"IO_Basis", "  ioops=", ioops, "  iregion=", iregion, "  karran=", karran
   !WRITE(*,*)" eint(1+nchanacc)=", eint(1+nchanacc), "irdindep=", irdindep
ELSEIF(iwrindep)THEN
   !WRITE(*,*)"IO_Basis", "  ioops=", ioops, "  iregion=", iregion, "  karran=", karran
   !WRITE(*,*)" eint(1+nchanacc)=", eint(1+nchanacc), "irdindep=", irdindep
   !WRITE(*,*)" ioops=",ioops,"iregion=",iregion, " Loc6"
   xktot=usys2*Etot
   IF(ioops)THEN
      WRITE(bas_unit)nt, nnuj(karran), nchanacc, nnuj(karran)/nhermt(karran), nhermt(karran), ioops, iregion, Rhoval
   ELSE   
      WRITE(bas_unit)nt, nnuj(karran), nchanacc, nnuj(karran)/nhermt(karran), nhermt(karran), ioops, iregion
   ENDIF
   WRITE(bas_unit)xktot,ioops
   WRITE(bas_unit)(mvib3(i+nchanacc),jrot3(i+nchanacc),lorb3(i+nchanacc),   &
                   mega3(i+nchanacc),kchan(i+nchanacc),xksq3(i+nchanacc),   &
                   nuj3(i+nchanacc), i=1,nt)                
   IF(karran==1)THEN
     Index=0
   ELSEIF(karran==2)THEN
     Index=nnuj(1)
   ELSEIF(karran==3)THEN
     Index=nnuj(1)+nnuj(2)
   ENDIF
   WRITE(bas_unit)(chinuj(Index+i), i=1,nnuj(karran))
   WRITE(bas_unit)(xsq(i+nchanacc),xk(i+nchanacc),eint(i+nchanacc), i=1,nt)
ENDIF

RETURN
ENDSUBROUTINE IO_Basis
