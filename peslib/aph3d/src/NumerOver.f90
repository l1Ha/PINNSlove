SUBROUTINE Numerover(nmax, mu, hstep, v, E, Psi, Diff, nu)
USE Numeric_Kinds_Module
USE Numbers_Module
USE FileUnits_OneDim_Module
IMPLICIT NONE
SAVE
INTEGER n, nmax, LeftTurn, RightTurn, nu, i, nnodes
REAL(KIND=WP_Kind) E, mu, hstep, tn, qn, Scale, Norm, Diff, diffsave
REAL(KIND=WP_Kind) Psi(0:nmax), v(0:nmax)
REAL(KIND=WP_Kind), ALLOCATABLE:: u(:), PsiRight(:), PsiLeft(:)

IF(ALLOCATED(u))THEN
  IF(SIZE(Psi).ne.nmax+1)STOP "Num"
ELSE
   ALLOCATE(u(0:nmax+1),PsiRight(0:nmax+1),PsiLeft(0:nmax+1))
ENDIF

LeftTurn=0
RightTurn=nmax
DO n=1,nmax-1
   IF(LeftTurn==0.and.((v(n)-E)*(v(n+1)-E))<Zero)LeftTurn=n  
   IF(n>LeftTurn.and.RightTurn==nmax.and.((v(n)-E)*(v(n+1)-E))<Zero)THEN
      RightTurn=n+1
      EXIT
   ENDIF
ENDDO
LeftTurn = MAX(1,LeftTurn)
RightTurn= MIN(nmax-1,RightTurn)

! Outward propagation
u=Zero
Psi=Zero
PsiRight=Zero
u(0)=Zero
u(1)=One
DO n=1,rightturn
   tn=(hstep**2/Twelve)*(Two*mu)*(v(n)-E)
   IF(tn.ge.Zero)THEN
      IF(SQRT(Twelve*tn)<700.d0)THEN
         qn=Two*COSH(SQRT(Twelve*tn))
      ELSE
         qn=(Two+Ten*tn)/(One-tn)
      ENDIF
   ELSE
      qn=Two*COS(SQRT(-Twelve*tn))
   ENDIF
    psi(n)=u(n)*qn/(Two+Ten*tn)
    psi(n)=u(n)/(One-tn)
    u(n+1)=qn*u(n)-u(n-1)
   IF(ABS(u(n+1)).gt.1.d+20)THEN    ! If wavefunction gets too large rescale the wavefunction
      Scale=u(n+1)
      u=u/Scale
      Psi=Psi/Scale
   ENDIF
ENDDO

! Determine the norm of the wavefunction in the classically allowed region
Norm=hstep*Dot_Product(Psi(LeftTurn:RightTurn),Psi(LeftTurn:RightTurn))
Psi=Psi/Sqrt(Norm)
PsiLeft=Psi

! Inward propagation
u=Zero
PsiRight=Zero
u(nmax)=Zero
u(nmax-1)=One
DO n=nmax-1,LeftTurn,-1
   tn=(hstep**2/Twelve)*(Two*mu)*(v(n)-E)
   IF(tn.ge.Zero)THEN
      IF(SQRT(Twelve*tn)<700.d0)THEN
         qn=Two*COSH(SQRT(Twelve*tn))
      ELSE
         qn=(Two+Ten*tn)/(One-tn)
      ENDIF
   ELSE
      qn=Two*COS(SQRT(-Twelve*tn))
   ENDIF
   psiright(n)=u(n)*qn/(Two+Ten*tn)
   psiright(n)=u(n)/(One-tn)
   u(n-1)=qn*u(n)-u(n+1)
   IF(ABS(u(n-1)).gt.1.d+20)THEN
      Scale=u(n-1)
      u=u/Scale
      PsiRight=PsiRight/Scale
   ENDIF
ENDDO

! Determine the norm of the wavefunction in the classically allowed region
Norm=hstep*Dot_Product(PsiRight(LeftTurn:RightTurn),PsiRight(LeftTurn:RightTurn))
PsiRight=PsiRight/Sqrt(Norm)

! Determine the relative phase difference of the outward and inward wavefunction
Scale=Sign(One,Dot_Product(Psi(LeftTurn:RightTurn),PsiRight(LeftTurn:RightTurn)))
PsiRight=PsiRight*Scale

! Determine the difference between the outward and inward propagated wavefunction
u=Zero
DO n=LeftTurn+1,RightTurn
  u(n)=(psi(n)*psiright(n-1)-psiright(n)*psi(n-1))/(psi(n)*psiright(n-1)+psiright(n)*psi(n-1))
  u(n)=ABS(u(n))
ENDDO
Diff=Sqrt(Dot_Product(u(LeftTurn:RightTurn),u(LeftTurn:RightTurn)))

! Combine the outward and inward wavefunction into a single wavefunction
DO n=0,nmax
   IF(n.lt.LeftTurn)THEN
      Psi(n)     =Psi(n)
   ELSE IF(n.gt.RightTurn)THEN
      Psi(n)     =PsiRight(n)
   ELSE
      Psi(n)     =(Psi(n)+PsiRight(n))/Two
   ENDIF
ENDDO

! This norm is calculated over the full range
Norm=hstep*Dot_Product(Psi,Psi)
Psi=Psi/Sqrt(Norm)
PsiRight=PsiRight/Sqrt(Norm)
PsiLeft=PsiLeft/Sqrt(Norm)
Diff=Diff

DO n=0,nmax
   IF(ABS(Psi(n))<1.d-30)THEN
      PSi(n)=Zero
      PSiLeft(n)=Zero
      PSiRight(n)=Zero
   ENDIF
ENDDO

! Count the number of nodes in the classically allowed region
nnodes=0
DO n=1,nmax
   IF(n>LeftTurn.and.n<RightTurn.and.Psi(n)*Psi(n-1)<Zero)THEN
      nnodes=nnodes+1
   ENDIF
ENDDO

!IF(nu/=nnodes)THEN
   !WRITE(*,*)"nu=",nu,"  Nnodes=",nnodes
   !STOP "In Numerover: nu/=nnodes"
!ENDIF


DIFFSave=Diff

RETURN
END
