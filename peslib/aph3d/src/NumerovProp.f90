SUBROUTINE NumerovProp(nmax, mu, hstep, v, Emin, Emax, NEnergy, rvals, Psi)
USE Numeric_Kinds_Module
USE Numbers_Module
USE FileUnits_Common_Module
USE FileUnits_OneDim_Module
USE QState_Module
USE DiatomicPot_Module, ONLY: DiatomicPot
USE CriticalDist_Module, ONLY: vasy
IMPLICIT NONE
INTEGER n, nmax, LeftTurn, KEnergy, NEnergy
REAL(KIND=WP_Kind) E, mu, hstep, tn, qn, Scale, Emin, Emax, DeltaE, PlotPsi, SCLength
REAL(KIND=WP_Kind) z, zjn, zjnp, zyn, zynp, a, b, kmat, k, phase, c, d, rmat, f
REAL(KIND=WP_Kind) Psi(0:nmax), v(0:nmax), rvals(0:nmax)
REAL(KIND=WP_Kind), ALLOCATABLE:: u(:)
COMPLEX(KIND=WP_Kind) smat
COMPLEX(KIND=WP_Kind), PARAMETER:: i=(zero,one), zone=(one,zero)

IF(ALLOCATED(u))THEN
  IF(SIZE(Psi).ne.nmax+1)STOP "Num"
ELSE
   ALLOCATE(u(0:nmax))
ENDIF

OPEN(Unit=Numerov_Prop_Unit,File=TRIM(OutDIR)//"DiatomicPlots/"//TRIM(DiatomicPot)//"_Numerov_Prop.csv")
OPEN(Unit=Continuum_Unit,File=TRIM(OutDIR)//"DiatomicPlots/"//TRIM(DiatomicPot)//"_Continuum.csv")
OPEN(Unit=ContinuumBin_Unit,File=TRIM(OutDIR)//"BinOut/Continuum.bin",FORM="unformatted")

IF(jval==0)THEN
   WRITE(Continuum_Unit,'(A5,7A17)')"jval,  ", "E(au),","-Log10(E-Vasy),", "K,  ", "Real(S),  ", "Imag(S),  ", "Phase,  ", "SCLength"
ELSE
   WRITE(Continuum_Unit,*)"&"
ENDIF
DeltaE=(EMax-EMin)/NEnergy
DO KEnergy=NEnergy-1,0,-1
   E=EMin+DeltaE*KEnergy+vasy
   u=Zero
   Psi=Zero
   u(0)=Zero
   u(1)=One
   LeftTurn=0
   DO n=1,nmax-1
      tn=(hstep**2/Twelve)*(Two*mu)*(v(n)-E)
      IF(tn.ge.Zero)THEN
         IF(SQRT(Twelve*tn)<700.d0)THEN
            qn=Two*COSH(SQRT(Twelve*tn))
         ELSE
            qn=One
            u(n-1)=Zero
            u(n)=one
         ENDIF
      ELSE
         qn=Two*COS(SQRT(-Twelve*tn))
         IF(LeftTurn.eq.0)LeftTurn=n-1 ! Find left turning point
      ENDIF
      u(n+1)=qn*u(n)-u(n-1)
      psi(n)=u(n)*qn/(Two+Ten*tn)
      IF(u(n).gt.1.d+20)THEN    ! If wavefunction gets too large rescale the wavefunction
         Scale=u(n)
         u=u/Scale
         Psi=Psi/Scale
      ENDIF
   ENDDO
   k=SQRT(two*mu*E)
   z=k*rvals(nmax-1)
   CALL rbes (jval,z,zjn,zjnp,zyn,zynp)
   a=z*zjn/SQRT(k)
   b=z*zyn/SQRT(k)
   
   z=k*rvals(nmax-2)
   CALL rbes (jval,z,zjn,zjnp,zyn,zynp)
   c=z*zjn/SQRT(k)
   d=z*zyn/SQRT(k)
   
   rmat=psi(nmax-1)/psi(nmax-2)
   
   kmat=(a-c*rmat)/(b-d*rmat)
   smat=(zone+i*kmat)/(zone-i*kmat)
   phase=ATAN2(REAL(smat),IMAG(smat))/two
   SCLength=-Tan(Phase)/k
   WRITE(Continuum_Unit,'(i5,",",7(es17.9,","))')jval, E, -LOG10(Emin+DeltaE*KEnergy), kmat, smat, phase, SCLength
   f=(a*d-c*b)/(psi(nmax-1)*d-psi(nmax-2)*b)
   psi=psi*f
   IF(kenergy==0.and.jval==0)THEN
      WRITE(Numerov_Prop_Unit,'(1x,A5,2A16,A10,I5,",",2(A10,ES15.7,","))')" n,"," rvals(n),"," PlotPsi,", "jval=",jval, "E=",E, "E-vasy=",DeltaE*KEnergy+vasy
   ELSE
      ! PORT FIX (2026-10-01): the original format '(...3(",")...)' consumed
      ! list items of the wrong type (integer as A); rebuilt to match the list.
      WRITE(Numerov_Prop_Unit,'(1x,A1,A5,I5,1x,A2,ES15.7,1x,A7,ES15.7)')"&", "jval=",jval, "E=",E, "E-vasy=", Emin+DeltaE*KEnergy
   ENDIF
   DO n=0,nmax-1
      PlotPsi=Psi(n)
      IF(ABS(PlotPsi)<1.d-30)PlotPsi=Zero
      WRITE(Numerov_Prop_Unit,'(I5,",",ES15.7,",",ES15.7,",",ES15.7)')n, rvals(n), PlotPsi
   ENDDO
   IF(Kenergy==0)THEN
      WRITE(ContinuumBin_Unit)rvals
      WRITE(ContinuumBin_Unit)v
   ENDIF
   WRITE(ContinuumBin_Unit)jval, E, kmat, smat, phase
   WRITE(ContinuumBin_Unit)psi
ENDDO


RETURN
ENDSUBROUTINE NumerovProp
