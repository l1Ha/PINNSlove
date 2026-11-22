SUBROUTINE Equilth (k,rho,we,wexe,usys,zet,seq, th3,delt,g2,g3)
!
!  This routine calculates the delves' theta where the oscillators
!  are centered; that is, where z=a tan(theta)-b cot(theta) +c is
!  zero.  It puts this at the maximum of the wavefunction for
!  a Morse oscillator determined from the asymptotic diatomic
!  except that (1) things may be shifted by the empirical
!  scaling factors
!  used, and (2) it is modified by terms which assure that the proper
!  boundary conditions are satisfied at theta = 0 and pi/2.
!  In May 1991, the additional factor was (sin 2 theta)**zet.
!  In May 1992, the additional factor was (cos theta)**zet.
!  atomic units.
!  R. T Pack.
!
USE Numeric_Kinds_Module
USE FileUnits_Module
USE Numbers_Module
IMPLICIT NONE
INTEGER k, i, ithcall, isav, j
REAL(Kind=WP_Kind) rho, we, wexe, usys, zet, seq, cosm
REAL(Kind=WP_Kind) sin, cos, asin, th3, fudge
REAL(Kind=WP_Kind) betm, d, d2, exp
REAL(Kind=WP_Kind) f1, f2, f3, test, si, co, co2, co3, delt, thetam
REAL(Kind=WP_Kind) ex, ga2, ga3, g2, g3, tan, gamm, exu
REAL(Kind=WP_Kind) th1(3), th2(3)
DATA th1, th2 /6*zero/
DATA ithcall /0/
DATA test /1.d-9/
!
WRITE(Out_Unit,*)'begining of routine equilth'
!  begin executable statements
IF(ithcall==0)THEN
   ithcall=1
ENDIF
!
!  determine parameters of asymptotic Morse wavefunction
betm = sqrt(two*usys*wexe)
d = we/(two*wexe)
!
!  parameters of Morse wavefunction in theta, at finite rho.
!  thetam is theta of minimum of Morse potential.
thetam=halfpi
cosm=zero
IF(rho>seq)THEN
   thetam=asin(seq/rho)
   cosm=cos(thetam)
   fudge = delt +(one-delt)*cosm
   gamm=betm*rho*fudge
   d2=d/fudge
ELSE
   exu=exp(-betm*(rho-seq))
   fudge=exu*(exu-one)*seq/betm + (rho*delt)**2
   fudge=sqrt(fudge)
   gamm=betm*fudge
   d2=d*rho/fudge
ENDIF
!  initial estimates
IF(th1(k)==zero)THEN
   th1(k)=0.9*halfpi
   IF(thetam<th1(k)) th1(k)=thetam
ENDIF
f1=d2*(one-exp(-gamm*(th1(k)-thetam))) - half + zet*tan(th1(k))/gamm
!
IF(th2(k)==zero) th2(k)=half*(fourthpi+th1(k))
f2=d2*(one-exp(-gamm*(th2(k)-thetam))) - half + zet*tan(th2(k))/gamm
!  put the better estimate into th2
IF(ABS(f1)<ABS(f2))THEN
   f3=f2
   th3=th2(k)
   f2=f1
   th2(k)=th1(k)
   f1=f3
   th1(k)=th3
ENDIF
!  use secant method to find the zero of the slope function.
DO i=1,50
   isav=i
   th3=(th1(k)*f2-th2(k)*f1)/(f2-f1)
   IF(th3>halfpi) th3=halfpi-(1.d-7)*i
   IF(th3<zero) th3=(1.d-7)*i
   f3=d2*(one-exp(-gamm*(th3-thetam))) - half + zet*tan(th3)/gamm
!  try to prevent overshooting
   DO j=1,50
      IF(ABS(f3)<ABS(f2))GOTO 5
      th3=half*(th2(k)+th3)
      f3=d2*(one-exp(-gamm*(th3-thetam))) - half + zet*tan(th3)/gamm
   ENDDO
   STOP 'overshot in Equilith'
 5    CONTINUE
!
!        WRITE(Out_Unit,61) k, isav, th3, f3
   IF(ABS(f3)<test) GOTO 3
!  this leaves next to last estimates in as first estimates for the
!  next rho to keep them spaced a bit.  To get the last two, move
!  the above statement to just before 2
   IF(ABS(f1)<ABS(f2)) GOTO 1
   th1(k)=th2(k)
   f1=f2
    1    CONTINUE
   th2(k)=th3
   f2=f3
ENDDO

 3 WRITE(Out_Unit,61) k, isav, th3, f3
61    FORMAT(1x,' for channel ',i3,' after ',i3,' iters, theq=',d14.7,' and err=',d14.7)
IF(ABS(f3)>test) STOP 'equilth'
!
!  now calculate the second and third derivatives at thetaeq for use
!  in determining the parameters of z.
!
si=sin(th3)
co=cos(th3)
co2=co*co
co3=co2*co
ex=exp(-gamm*(th3-thetam))
ga2=gamm*gamm
ga3=ga2*gamm
g2=d2*ex*ga2+zet/co2
g3=-d2*ex*ga3+two*zet*si/co3
!
WRITE(Out_Unit,*)'END of routine equilth'
RETURN
ENDSUBROUTINE Equilth
