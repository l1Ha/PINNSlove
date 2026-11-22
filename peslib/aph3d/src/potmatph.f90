SUBROUTINE potmatph (n, rmid, w, hstep, ider)
USE FileUnits_MODULE
USE nstate_MODULE
USE Narran_Module
USE Qall_Module
USE Storage_Module
USE melm_MODULE
USE region_MODULE
USE Oops_Module
USE etamat_MODULE
USE dip_MODULE
IMPLICIT NONE
!---------------------------------------------------------------------
! written by g. a. parker
! calculates the interaction matrix and its first and second
! derivatives.
!---------------------------------------------------------------------
LOGICAL little, medium, full, numder
INTEGER ithcll, ithsub, ideriv, ider, n, kk, nnp1, ii, iii
INTEGER i, index, nterm, ms, in, j, ij
REAL(Kind=WP_Kind) rold, r, rmid, hstep, rsqinv, xksq, rho, ddr2, dr2, r2
REAL(Kind=WP_Kind) w(*)   ! PORT: assumed-size
! PORT (2026-10-02): optional dump of the coupled-channel interaction matrix
! W(r) in the surface-function basis, for building/validating the 2x2 inner
! propagator.  Enabled by APH3D_WDUMP=1, window APH3D_WDUMP_LO/HI.
LOGICAL wdump_on, wdump_first
REAL(Kind=WP_Kind) wdump_lo, wdump_hi
CHARACTER(LEN=64) wdump_env
INTEGER wdump_len, wdump_stat
DATA wdump_first/.true./
REAL(Kind=WP_Kind) ww0(nstate,nstate),ww1(nstate,nstate),ww2(nstate,nstate)
REAL(Kind=WP_Kind) w0(nstate*nstate), w1(nstate*nstate), w2(nstate*nstate)
REAL(Kind=WP_Kind)g0(nstate,nstate),g1(nstate,nstate), g2(nstate,nstate)
DATA rold /-995.d+21/
DATA ithcll/0/,ithsub/0/
DATA little/.false./,medium/.false./,full/.false./
CALL popt('potmat  ',little,medium,full,ithcll,ithsub)

IF(wdump_first)THEN
   wdump_first = .false.
   wdump_on = .false.; wdump_lo = -1.d10; wdump_hi = 1.d10
   CALL get_environment_variable('APH3D_WDUMP', wdump_env, wdump_len, wdump_stat)
   IF(wdump_stat==0 .AND. wdump_len>0)THEN
      IF(wdump_env(1:1)=='1' .OR. wdump_env(1:1)=='T' .OR. wdump_env(1:1)=='t') wdump_on=.true.
   ENDIF
   CALL get_environment_variable('APH3D_WDUMP_LO', wdump_env, wdump_len, wdump_stat)
   IF(wdump_stat==0 .AND. wdump_len>0) READ(wdump_env(1:wdump_len),*) wdump_lo
   CALL get_environment_variable('APH3D_WDUMP_HI', wdump_env, wdump_len, wdump_stat)
   IF(wdump_stat==0 .AND. wdump_len>0) READ(wdump_env(1:wdump_len),*) wdump_hi
ENDIF

ideriv=ider
IF(ider==-1) ideriv=0
r=rmid
nnp1=n*(n+1)/2

IF(irdindep)THEN
READ(frst_unit) (w(kk),kk=1,nnp1)
ii=0
DO iii=1,n
   ii=ii+iii
   w(ii)=w(ii)-xktot
ENDDO
CALL Wdump2(n, w, r, ider)
RETURN
ENDIF

!      WRITE(Out_Unit,*)'potmat: iregion=',iregion
IF(iregion/='aph    ')THEN
IF(.NOT.numder) GOTO 20
IF(ideriv==0) r=rmid
IF(ideriv==1) r=rmid-0.5d0*hstep
IF(ideriv==2) r=rmid+0.5d0*hstep
GOTO 30
!20 IF(ideriv-1) 30,100,130
20 IF(ideriv==1)goto 100
IF(Ideriv==2)goto 130
30 IF(r==rold) GOTO 50
!---------------------------------------------------------------------
! determine the interaction matrix
!---------------------------------------------------------------------
rold=r
CALL potvib (n,r,w0,w1,w2)
IF(ideriv/=0) GOTO 70
IF(r/=0.d0) rsqinv=1.d0/(r*r)
IF(r==0.d0) rsqinv=1.d+30
!----------------------------------------------------------------------
! determine the diagonal contributions.
!----------------------------------------------------------------------
index=0
DO i=1,n
index=index+i
xksq=xksq3(i)
IF(iregion=='jacobi ')THEN
IF(melmnt=='probf') w0(index)=(w0(index)-xksq)
IF(melmnt=='space') w0(index)=(w0(index)-xksq)+lorb3(i)*(lorb3(i)+1)*rsqinv
IF(melmnt=='molbf') w0(index)=(w0(index))+lorb3(i)*(lorb3(i)+1)*rsqinv
IF(full)WRITE(Out_Unit,*)i,xksq,index,w0(index),lorb3(i),rsqinv
ELSE
rho=r
xksq=xktot+(xksq-xktot)*(rhoj/rho)**2+0.25d0/rho**2
w0(index)=(w0(index)-xksq)
ENDIF
ENDDO
!----------------------------------------------------------------------
! now include the off-diagonal contributions.
!----------------------------------------------------------------------
50 nterm=n*(n+1)/2
DO i=1,nterm
IF(melmnt=='probf') w(i)=w0(i)+etamx(i)*rsqinv
IF(melmnt=='space') w(i)=w0(i)
IF(melmnt=='molbf') w(i)=w0(i)+etamx(i)
ENDDO
ms=1
GOTO 160
!----------------------------------------------------------------------
! determine the interaction potential and its derivatives.
! first calculate the diagonal contributions.
!----------------------------------------------------------------------
70 index=0
DO i=1,n
index=index+i
xksq=xksq3(i)
IF(iregion=='jacobi ')THEN
IF(melmnt=='probf') w0(index)=(w0(index)-xksq)
IF(melmnt=='space') w0(index)=(w0(index)-xksq)+lorb3(i)*(lorb3(i)+1)*rsqinv
IF(melmnt=='molbf') w0(index)=(w0(index))+lorb3(i)*(lorb3(i)+1)*rsqinv
ELSE
rho=r
xksq=xktot+(xksq-xktot)*(rhoj/rho)**2+0.25d0/rho**2
w0(index)=(w0(index)-xksq)
ENDIF
ENDDO
!---------------------------------------------------------------------
! now include the off-diagonal contributions.
!---------------------------------------------------------------------
rsqinv=1.d0/(r*r)
in=0
DO j=1,n
DO i=1,j
in=in+1
ij=i+(j-1)*n
w(ij)=0.0d0
IF(kchan(i)/=kchan(j)) GOTO 90
ji=j+(i-1)*n
IF(melmnt=='probf') w(ij)=w0(in)+etamx(in)*rsqinv
IF(melmnt=='space') w(ij)=w0(in)
IF(melmnt=='molbf') w(ij)=w0(in)+etamx(in)
90 w(ji)=w(ij)
ENDDO
ENDDO
GOTO 160
100 IF(r==rold) GOTO 110
!---------------------------------------------------------------------
! determine the first derivative of the interaction matrix
!---------------------------------------------------------------------
rold=r
CALL potvib (n,r,w0,w1,w2)
110 dr2=-2.d0/(r*r*r)
in=0
DO j=1,n
DO i=1,j
in=in+1
ij=i+(j-1)*n
w(ij)=0.0d0
IF(kchan(i)/=kchan(j)) GOTO 120
ji=j+(i-1)*n
IF(iregion=='jacobi ')THEN
IF(melmnt=='probf') w(ij)=w1(in)+etamx(in)*dr2
IF(i==j)THEN
IF(melmnt=='space') w(ij)=w1(in)+lorb3(j)*(lorb3(j)+1)*dr2
IF(melmnt=='molbf') w(ij)=w1(in)+lorb3(j)*(lorb3(j)+1)*dr2
ELSE
IF(melmnt=='space') w(ij)=w1(in)
IF(melmnt=='molbf') w(ij)=w1(in)
ENDIF
ELSE
rho=r
xksq=(xksq3(i)-xktot)*(rhoj/rho)**2*(-2/rho)
xksq=xksq-0.5d0/rho**3
IF(i==j)w(ij)=w1(in)
IF(i/=j)w(ij)=w1(in)
ENDIF
120 w(ji)=w(ij)
ENDDO
ENDDO
ms=0
GOTO 160
!---------------------------------------------------------------------
! determine the second derivative of the interaction matrix
!---------------------------------------------------------------------
130 IF(r==rold) GOTO 140
rold=r
CALL potvib (n,r,w0,w1,w2)
140 r2=r*r
ddr2=6.d0/(r2*r2)
in=0
DO j=1,n
DO i=1,j
in=in+1
ij=i+(j-1)*n
w(ij)=0.0d0
IF(kchan(i)/=kchan(j)) GOTO 150
ji=j+(i-1)*n
IF(iregion=='jacobi ')THEN
IF(melmnt=='probf') w(ij)=w2(in)+etamx(in)*ddr2
IF(i==j)THEN
IF(melmnt=='space') w(ij)=w2(in)+lorb3(j)*(lorb3(j)+1)*ddr2
IF(melmnt=='molbf') w(ij)=w2(in)+lorb3(j)*(lorb3(j)+1)*ddr2
ELSE
IF(melmnt=='space') w(ij)=w2(in)
IF(melmnt=='molbf') w(ij)=w2(in)
ENDIF
ELSE
rho=r
xksq=(xksq3(i)-xktot)*(rhoj/rho)**2*(6/rho**2)
xksq=xksq-1.5d0/rho**4
IF(i==j)w(ij)=w1(in)
IF(i/=j)w(ij)=w1(in)
ENDIF
150 w(ji)=w(ij)
ENDDO
ENDDO
!---------------------------------------------------------------------
! read in aph potential matrix elements, dipole matrix elements
!---------------------------------------------------------------------
ELSE
CALL aphmat(rmid,n,w,ww0,ww1,ww2,ider)
!
!        WRITE(Out_Unit,*)'lphoto,rmid,enddip=',lphoto,rmid,enddip
IF(lphoto.AND.rmid<=enddip)THEN
    CALL aphmatph(rmid,n,g0,g1,g2)
ENDIF
!
ENDIF
ms=0
160 IF(iwrindep.AND.(ider==-1))THEN
    WRITE(frst_unit) (w(kk),kk=1,nnp1)
ENDIF

IF(full)THEN
IF(ideriv==0) ms=1
IF(ideriv==0) WRITE(Out_Unit,170) r
IF(ideriv==1) WRITE(Out_Unit,180) r
IF(ideriv==2) WRITE(Out_Unit,190) r
CALL MxOutD(w,n,n,1)
ENDIF

IF(wdump_on .AND. ideriv==0 .AND. r>=wdump_lo .AND. r<=wdump_hi)THEN
   OPEN(Unit=94, File='wdump.txt', Status='UNKNOWN', Position='APPEND')
   WRITE(94,'(A,ES22.14,A,I6,A,I4)') '# W(r) at r=', r, '  n=', n, '  ider=', ider
   DO i = 1, n
      WRITE(94,'(200ES22.14)') (w(i+(j-1)*n), j = 1, n)
   ENDDO
   CLOSE(94)
ENDIF
IF(ideriv==0) CALL Wdump2(n, w, r, ider)
IF(ideriv==0)CALL enlvls(r,n,w)
RETURN
!----------------***end-potmat***---------------------------------------
!
170 FORMAT(1x,"the potential at r=",f10.5)
180 FORMAT(1x,"the first derivative of the potential at r=",f10.5)
190 FORMAT(1x,"the second derivative of the potential at r=",f10.5)

CONTAINS
SUBROUTINE Wdump2(nn, ww, rr, idd)
! PORT (2026-10-02): gated dump of the coupled-channel matrix W(r) in the
! surface-function basis (APH3D_WDUMP / APH3D_WDUMP_LO / APH3D_WDUMP_HI).
IMPLICIT NONE
INTEGER, INTENT(IN) :: nn, idd
REAL(Kind=WP_Kind), INTENT(IN) :: ww(*), rr
INTEGER :: iw, jw
IF(.NOT.wdump_on) RETURN
IF(rr < wdump_lo .OR. rr > wdump_hi) RETURN
OPEN(Unit=94, File='wdump.txt', Status='UNKNOWN', Position='APPEND')
WRITE(94,'(A,ES22.14,A,I6,A,I4)') '# W(r) at r=', rr, '  n=', nn, '  ider=', idd
DO iw = 1, nn
   WRITE(94,'(500ES22.14)') (ww(iw+(jw-1)*nn), jw = 1, nn)
ENDDO
CLOSE(94)
END SUBROUTINE Wdump2
END
