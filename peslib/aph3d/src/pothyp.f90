SUBROUTINE pothyp (arg1, arg2, arg3, pot)
USE FileUnits_Module
USE Arrch_Module
USE Masses_Module
USE Numeric_Kinds_Module
REAL(Kind=WP_Kind) r(3)
REAL(Kind=WP_Kind) arg1, arg2, arg3, cgamaa, cgamab, cgamac, pot
REAL(Kind=WP_Kind) r1, r2, r3, RA, RB, RC, RA0, RB0, RC0
REAL(Kind=WP_Kind) RAB, RBC, RCA, U

IF(karran==1)THEN
   ra=arg1
   rbc=arg2
   cgamaa=arg3
   rb0=rbc*mass(3)/(mass(2)+mass(3))
   rc0=rbc*mass(2)/(mass(2)+mass(3))
   r1=rbc
   r2=sqrt(ra**2+rc0**2-2*ra*rc0*cgamaa)
   r3=sqrt(ra**2+rb0**2+2*ra*rb0*cgamaa)
ELSEIF(karran==2)THEN
   rb=arg1
   rca=arg2
   cgamab=arg3
   rc0=rca*mass(1)/(mass(3)+mass(1))
   ra0=rca*mass(3)/(mass(3)+mass(1))
   r1=sqrt(rb**2+rc0**2+2*rb*rc0*cgamab)
   r2=rca
   r3=sqrt(rb**2+ra0**2-2*rb*ra0*cgamab)
ELSEIF(karran==3)THEN
   rc=arg1
   rab=arg2
   cgamac=arg3
   ra0=rab*mass(2)/(mass(1)+mass(2))
   rb0=rab*mass(1)/(mass(1)+mass(2))
   r1=sqrt(rc**2+rb0**2-2*rc*rb0*cgamac)
   r2=sqrt(rc**2+ra0**2+2*rc*ra0*cgamac)
   r3=rab
ENDIF
r(1)=r1
r(2)=r2
r(3)=r3
CALL surface (u,r)
pot=u
RETURN
END
