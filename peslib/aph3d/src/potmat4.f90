SUBROUTINE PotMat4(nchanl, rho, wmat, hstep, ideriv)
USE Numeric_Kinds_Module
USE TotalEng_Module
USE jtot_Module
USE fileunits_Module
USE boundary_Module
USE quantum_numbers_Module
USE Masses_Module
USE psi_Module
USE numbers_Module
USE int1d_Module
USE popt_Module
USE Time_Module
IMPLICIT NONE
CHARACTER(LEN=11) :: Blank
CHARACTER(LEN=21), PARAMETER:: ProcName='potmat4'
!---------------------------------------------------------------------
! written by g. a. parker
! calculates the interaction matrix
!---------------------------------------------------------------------
INTEGER nchanl, i, ideriv
REAL(Kind=WP_Kind) dt1, dt2, hstep, pmat, ksquared
REAL(Kind=WP_Kind) work(nltheta+1), bfpotmat(nstate,nstate,0:mxjrot)
REAL(Kind=WP_Kind) pj(nbtheta,0:mxjrot,0:mxjrot)
REAL(Kind=WP_Kind) vprime(nbtheta,nltheta+1), sfpotmat(nstate,nstate)
REAL(Kind=WP_Kind) xpt(nbtheta), weight(nbtheta)
REAL(Kind=WP_Kind)  rho, wmat(nchanl,nchanl)
REAL(Kind=WP_Kind) timeinfo
!  HP's won't allow a matrix as big as pmat unless in a common block.
!  This array is called pmat here.  It is vleg in potmatn, vleg_mat,
!  and bfpmatrix.
COMMON /bigblock/ pmat(nltheta+1,0:mxjrot,0:mxjrot,0:mxjrot)
SAVE /bigblock/

IF(ideriv/=0)THEN
   WRITE(msg_unit,*)'Error ideriv/=0'
   STOP 'potmat'
ENDIF

!      dt1=(switch-lthetzro)/(n1-1)
!      dt2=(pi/two-switch)/(n2-1)
!  the above statements were corrected by RTP on 99/2/12
!  to make them consistent with q1d.  Up to this point, the
!  ltheta index runs on 0,nltheta.  It shifts in potmatn.
dt1=(switch-lthetzro)/n1
dt2=(pi/two-switch)/n2
!
DO i=1,nchanl
   IF(jrot(i)>mxjrot)THEN
      WRITE(Out_Unit,*)'Error jrot(i)>mxjrot',i,jrot(i)
      WRITE(msg_unit,*)'Error jrot(i)>mxjrot',i,jrot(i)
      STOP
   ENDIF
ENDDO

t1 = TimeInfo(ProcName,'PotMatN4')
WRITE(Out_Unit,*)' rho=', rho, ' rho_basis=', rho_basis
CALL potmatn4(jrot,lorb,nvib,dt1,dt2,psi_new,                     &
             pmat,nchanl,work,bfpotmat,sfpotmat,mxjrot,          &
             xpt, weight, rho, vprime, pj, jtot, j_min, symmetry)
t2 = TimeInfo(procname,'PotMatN4')
CALL WRITEtime(procname,t2-t1)
!      WRITE(msg_unit,*)' rho=', rho, ' rho_basis=', rho_basis
!
!      WRITE(Out_Unit,*)'ksquared, eigen_new(i),eigen_new(i)*usys2,',
!     >          'ksquared-eigen_new(i)*usys2'
!      DO i=1,nchanl
!         WRITE(Out_Unit,10)ksquared, eigen_new(i),
!     >                     eigen_new(i)*usys2,
!     >                     ksquared-eigen_new(i)*usys2
!      ENDDO


ksquared=usys2*etot
!      WRITE(Out_Unit,*)'energies   xksq   wvec'
DO i=1,nchanl
   energies(i)=eigen_new(i)
   xksq(i)=usys2*(etot-energies(i))
   wvec(i)=sqrt(abs(xksq(i)))
!         WRITE(Out_Unit,*)energies(i),xksq(i),wvec(i)
ENDDO
t1 = TimeInfo(ProcName,'Coupling_Matrix')
CALL coupling_matrix(nchanl, wmat, sfpotmat, eigen_new, ksquared, usys2, rho, rho_basis)
t2 = TimeInfo(procname,'Coupling_Matrix')
CALL WRITEtime(procname,t2-t1)
10   FORMAT(1x,5e15.7)

medium=.false.
IF(medium)THEN
   WRITE(Out_Unit,170) rho
   CALL MatrixOut(wmat, nchanl, nchanl, 'wmat',Blank, Blank,Blank,.False.,Blank,.False.)
ENDIF
!
t1 = TimeInfo(ProcName,'EnLvls')
CALL enlvls(rho, nchanl, wmat)
t2 = TimeInfo(procname,'EnLvls')
CALL WRITEtime(procname,t2-t1)
    
RETURN
!----------------***end-potmat***--------------------------------------
!
170 FORMAT(1x,"the potential at r=",f10.5)
ENDSUBROUTINE PotMat4
