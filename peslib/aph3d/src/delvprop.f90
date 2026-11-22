SUBROUTINE delvprop(method,idname,start,finish,nrchanl)
USE Numeric_Kinds_Module
USE Time_Conv_Module
USE dbs_Module
USE method_Module
USE FileNames_Module
!USE G1Matrix_Module
USE readWRITE_Module
USE TotalEng_Module
!
!  working is the file Rvwxyz in the working subdirectory into which the
!  R matrix will be put.
!
IMPLICIT NONE

!Input/Output variables
CHARACTER(LEN=11) :: Blank
CHARACTER(LEN=*) idname,method
REAL(Kind=WP_Kind) start, finish
INTEGER nrchanl

CHARACTER(LEN=21), PARAMETER:: ProcName='delvprop'
LOGICAL EDeriv, firstcall
LOGICAL debug
INTEGER i
data firstcall /.true./
SAVE firstcall

REAL(Kind=WP_Kind), ALLOCATABLE  ::  rmat(:,:), drmat(:,:)

IF(firstcall)THEN
   CALL get_debug(procname, debug)
   firstcall = .false.
ENDIF

!------------------------------------------------------------------
! The maximum energy will determine the wavelength
! of the most oscillatory state which will be used to set the
! step size in the log-derivative propagation.
!------------------------------------------------------------------

WRITE(Out_Unit,*)                  '#####################################'
WRITE(Out_Unit,*)                  '#      DelvesProp: Start            #'
WRITE(Out_Unit,*)                  '#          method: ',method,'       #'
WRITE(Out_Unit,'(1x,a19,e10.3,a6)')'#          Energy: ',etot  ,' Ha.   #'
WRITE(Out_Unit,*)                  '#            File: ',rmat_name(1:5),'          #'
WRITE(Out_Unit,*)                  '#####################################'
WRITE(Out_Unit,*)

WRITE(msg_unit,*)'Starting propagation for ',idname
ALLOCATE(rmat(nrchanl,nrchanl), drmat(nrchanl,nrchanl))
!
!  Set all elements of the LOG_derivative matrix initially to zero.
!  Then, set the diagonal elements to a large positive number.
!
rmat=0.d0
DO i=1,nrchanl
   rmat(i,i)=1.d+30
ENDDO
!
! Propagate the coupled equations from start  to finish.
!
IF(method(1:3)/='svd')THEN
   bmethod = 'dbs'
   CALL dbsprop(nrchanl,rmat,drmat,start,finish,nsteps,ksqmax,idname,method,frst_WRITE)
ELSE
   bmethod = 'svd'
   CALL svdprop(Etot, nrchanl, nrchanl, start, finish, rmat)
ENDIF
!
!   Write the R_Matrix to a file for subsequent analysis.
!
!   CALL WRITE_rmat(rmat, drmat, nrchanl, idname, finish)
!
WRITE(Out_Unit,*)'R_Matrix at finish ', finish,' Bohr'
CALL MatrixOut(rmat, nrchanl, nrchanl,'R_Matrix','R_Matrix at Finish',Blank, Blank,.False.,Blank,.False.)
!
!    IF(method=='dlogder ')THEN
!        EDeriv=.true.
!    ELSE
!        EDeriv=.false.
!    ENDIF
!      WRITE(Out_Unit,*)'here at loc 3: nopen=',nopen
!      CALL main_asymptotic(nrchanl, nopen, rmat, drmat, finish, EDeriv)
!      WRITE(Out_Unit,*)'here at loc 4: nopen=',nopen


DEALLOCATE(rmat, drmat)
IF(debug)THEN
   WRITE(Out_Unit,*)                '###################################'
   WRITE(Out_Unit,*)                '#      DelvesProp: End            #'
ENDIF

RETURN
end
