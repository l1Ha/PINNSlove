SUBROUTINE potth (nt, rho, u, up, upp)
!-----------------------------------------------------------------------
! calculates the matrix elements of (1/rho*cos(theta))**2.
! written by g. a. parker
!-----------------------------------------------------------------------
USE FileUnits_MODULE
USE Narran_Module
USE nstate_MODULE
USE parms_MODULE
USE approx_MODULE
USE tlogd_MODULE
USE VFunc_Module
USE Qall_Module
USE GaussQuady_Module
USE pow_MODULE
USE melm_MODULE
USE region_MODULE
USE Masses_Module
USE Oops_Module
USE thetas_MODULE
IMPLICIT NONE
LOGICAL little, medium, full
INTEGER ithcll, ithsub, nt, nqp, in, jp, megp, lbp, ichan, knujp, nq
INTEGER j, meg, lb, knuj, k, ms, n
REAL(Kind=WP_Kind) rho, rho2, rho3, rho4, cf, bf, sum, theta
REAL(Kind=WP_Kind) u(1), up(1), upp(1)
DATA ithcll/0/,ithsub/0/
DATA little/.false./,medium/.false./,full/.false./
CALL popt('potth   ',little,medium,full,ithcll,ithsub)
IF(iregion/='delves ')RETURN
rho2=rho*rho
rho3=rho*rho2
rho4=rho*rho3
in=0
DO nqp=1,nt
jp=jrot3(nqp)
megp=mega3(nqp)
lbp=lorb3(nqp)
ichan=kchan(nqp)
knujp=(nuj3(nqp)-1)*nhermt(ichan)
cf=c(ichan)/rhoj
bf=b(ichan)/rhoj
DO nq=1,nqp
j=jrot3(nq)
meg=mega3(nq)
lb=lorb3(nq)
knuj=(nuj3(nq)-1)*nhermt(ichan)
in=in+1
IF(kchan(nqp)/=kchan(nq))GOTO 300
IF(j/=jp)GOTO 300
IF(lb/=lbp)GOTO 300
sum=0.0
DO 290 k=1,nhermt(ichan)
theta=cf*xpth(k,ichan)+bf
sum=sum+wpth(k,ichan)*chinuj(knuj+k)*chinuj(knujp+k)/cos(theta)**2
290 CONTINUE
sum=sum*lb*(lb+1)
u(in)=u(in)+sum/usys2/rho2
up(in)=up(in)-2*sum/usys2/rho3
upp(in)=upp(in)+6*sum/usys2/rho4
300 CONTINUE
ENDDO
ENDDO
!-----------------------------------------------------------------------
! if desired print the matrix elements of the potential
! and its derivatives.
!-----------------------------------------------------------------------
IF(full)THEN
ms=1
WRITE(Out_Unit,410)
CALL MxOutL(  u, n, n, ms, melmnt, melmnt)
IF(numder)GOTO 310
WRITE(Out_Unit,420)
CALL MxOutL( up, n, n, ms, melmnt, melmnt)
WRITE(Out_Unit,430)
CALL MxOutL(upp, n, n, ms, melmnt, melmnt)
ENDIF
310 RETURN
!-----------------------------------------------------------------------
!----------------***end-potth***---------------------------------------
!-----------------------------------------------------------------------
410 FORMAT(//,"potential energy matrix elements")
420 FORMAT(//,"first derivative  of the potential energy matrix elements")
430 FORMAT(//,"second derivative of the potential energy matrix elements")
END
