SUBROUTINE Critical(rvals, hstep, mu, v, vp, vpp, jval, elec)
USE Numeric_Kinds_Module
USE Numbers_Module
USE CriticalDist_Module
USE FileUnits_Module
USE FileUnits_OneDim_Module
USE OneDim_Parms_Module
IMPLICIT NONE
LOGICAL root
INTEGER i, jval, elec, isigma0, ire, irb, isigma1, irc, ird
REAL(KIND=WP_Kind) r, mu, hstep, v0, v1, v2
REAL(KIND=WP_Kind) rvals(0:nmax), v(0:nmax), vp(0:nmax), vpp(0:nmax)

IF(jval.eq.jmin)CALL TimeInfoSub('Critical','Begin Critical')
sigma0=Zero
v0sigma0=Zero
v1sigma0=Zero
v2sigma0=Zero
sigma1=rvals(nmax)
v0sigma1=Zero
v1sigma1=Zero
v2sigma1=Zero
re=Zero
v0re=Zero
v1re=Zero
v2re=Zero
rb=Zero
v0rb=Zero
v1rb=Zero
v2rb=Zero
rc=rvals(nmax)
v0rc=Zero
v1rc=Zero
v2rc=Zero
rd=rvals(nmax)
v0rd=Zero
v1rd=Zero
v2rd=Zero
k=Zero
omega=Zero
alpha=Zero
beta=Zero
de=Zero
vasy=v(nmax)
v0rc=vasy

isigma0=-nmax
ire=-nmax
irb=-nmax
isigma1=-nmax
irc=-nmax
ird=-nmax

! For single well potentials there will be at most 2 distances which cross the asymptotic value of the potential
DO i=1,nmax
   CALL findroots(i, nmax, rvals, hstep, v, vasy, jval, mu, r, v0, v1, v2, root)
   IF(sigma0.eq.0.and.root)THEN
      sigma0=r ! The first distance where the potential crosses the asymptotic value of the potential is sigma0.
      v0sigma0=v0
      v1sigma0=v1
      v2sigma0=v2
      isigma0=i
   ELSE IF(root.and.irb>0)THEN
      sigma1=r ! The second distance where the potential crosses the asymptotic value of the potential is sigma1.
      v0sigma1=v0
      v1sigma1=v1
      v2sigma1=v2
      isigma1=i
   ENDIF

   ! For single well potentials there will be at most 1 extememe for the centrifical barrier.
   ! and
   ! For single well potentials there will be at most 1 minimum of the effective potential.
   IF(isigma0>0)CALL findextreme(i, nmax, rvals, hstep, vp, jval, mu, r, v0, v1, v2, root)
   IF(re.eq.0.and.root.and.isigma0>0)THEN
      re=r     ! Re is the equilibrium (minimum) of the potential.
      v0re=v0
      v1re=v1
      v2re=v2
      ire=i
   ELSE IF(root.and.irb>0)THEN
      rc=r    ! Rc is local maximum of the potential where the first derivative of the potential is zero
      v0rc=v0
      v1rc=v1
      v2rc=v2
      irc=i
   ENDIF

   ! For single well potentials there will be at most 2 Inflection Points (rb,rd)
   IF(ire>0)CALL findinflection(i, nmax, rvals, hstep, vpp, jval, mu, r, v0, v1, v2, root)
   IF(rb.eq.0.and.root.and.ire>0)THEN
      rb=r   ! Rb is the first inflection point.
      v0rb=v0
      v1rb=v1
      v2rb=v2
      irb=i
   ELSE IF(root.and.irb>0)THEN
      rd=r   ! Rd is the second inflection point.
      v0rd=v0
      v1rd=v1
      v2rd=v2
      ird=i
   ENDIF
ENDDO

IF((re.eq.Zero).and.(rb.eq.Zero).and.(rc.eq.Zero))THEN
   CALL TimeInfoSub('Critical','  End Critical')
   RETURN
ENDIF

WRITE(Out_Unit,*)
WRITE(Out_Unit,*)
WRITE(Out_Unit,'(a,i5,a,i5)')"Critical Distance in Bohr Information for jval=",jval," elec=", elec
WRITE(Out_Unit,'(4(a,es23.15))')'sigma0=',sigma0,'  V(sigma)=',v0sigma0,'  V1(sigma)=',v1sigma0, '  V2(sigma)=',v2sigma0
WRITE(Out_Unit,'(4(a,es23.15))')'    re=',    re,'     V(re)=',    v0re,'     V1(re)=',v1re,     '     V2(re)=',v2re
WRITE(Out_Unit,'(4(a,es23.15))')'    rb=',    rb,'     V(rb)=',    v0rb,'     V1(rb)=',    v1rb, '     V2(rb)=',v2rb
WRITE(Out_Unit,'(4(a,es23.15))')'sigma1=',sigma1,'  V(sigma)=',v0sigma1,'  V1(sigma)=',v1sigma1, '  V2(sigma)=',v2sigma1
WRITE(Out_Unit,'(4(a,es23.15))')'    rc=',    rc,'     V(rc)=',    v0rc,'     V1(rc)=',    v1rc, '     V2(rc)=',v2rc
WRITE(Out_Unit,'(4(a,es23.15))')'    rd=',    rd,'     V(rd)=',    v0rd,'     V1(rd)=',    v1rd, '     V2(rd)=',v2rd

IF(re.gt.Zero)THEN
   WRITE(Out_Unit,*)
   WRITE(Out_Unit,'(A)')"Morse Parameters"
   k=v2re
   omega=SQRT(k/mu)
   alpha=SQRT(mu*omega)
   de=MAX(v0rc,vasy)-v0re
   be=One/(Two*mu*re**2)
   beta=SQRT(v2re/(Two*ABS(de)))
   WRITE(Out_Unit,'(5(a,es23.15))')'   k=',   k, ' alpha=', alpha
   WRITE(Out_Unit,'(5(a,es23.15))')'   Be=', be, ' Omega=', omega
   WRITE(Out_Unit,'(5(a,es23.15))')'   De=', de, '  Beta=', beta, '  Re=',re, ' v2(rd)=', v0rd
ENDIF


IF(jval.eq.jmin.and.elec.eq.MinElec)THEN
   WRITE(Crit_Unit,'(2(a5,","),7(a6,","))')"elec","jval","sigma0","re","rb","sigma1","rc","rd","ri"
ENDIF
WRITE(Crit_Unit,'(2(I5,","),7(E23.15,","))')elec,jval,sigma0,re,rb,sigma1,rc,rd,ri

vmax(elec)=v0rc
vinf(elec)=vasy
vmin(elec)=v0re
rmin(elec)=re

RETURN
END SUBROUTINE Critical
