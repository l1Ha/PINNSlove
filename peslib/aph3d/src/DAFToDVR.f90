SUBROUTINE DAFToDVR(NStates)
USE Numeric_Kinds_Module
USE Numbers_Module
USE FileUnits_OneDim_Module
USE OneDim_Parms_Module
USE QstateDAF_Module
USE ElectronicState_Module
IMPLICIT NONE
CHARACTER(LEN=1), PARAMETER:: Jobz='V', UpLo='U', Range='A'
INTEGER NBound, NQuasi, NStates, NumNu, nYTerms, i, istate, nup, lwork, info, LocPsiMax
INTEGER jvalue, electronic, kstate
REAL(KIND=WP_Kind) :: Sum1, Fac, Sum0, Sum2, Sum3, Sum4 !, hstep !, rvals(0:nMax)
REAL(KIND=WP_Kind), ALLOCATABLE:: Yfit(:,:), Eig(:), Psi(:,:), PsiMax(:), DVR(:,:), work(:)
REAL(KIND=WP_Kind), ALLOCATABLE:: Abscissas(:), Weights(:), PsiDVR(:,:), ScalePsi(:)

Nyterms=MaxYterms
ALLOCATE(YFit(jmin:jmax,0:NYterms-1))
ALLOCATE(Eig(0:NStates-1))
ALLOCATE(Psi(0:NMax,0:NStates-1))
REWIND(DAFBin_Unit)
Yfit=Zero
IState=0
kstate=0
DO WHILE(IState.lt.NStates)
   READ(DAFBin_Unit)electronic,jvalue,numnu,nbound,nquasi,nYTerms
   DO Nup=0,NumNu-1
      READ(DAFBin_Unit)Nu,Eig(kState)
      READ(DAFBin_Unit)Psi(:,kState)
      IState=IState+1
      IF(jval==jvalue.and.elec==electronic)kstate=kstate+1
   ENDDO
   READ(DAFBin_Unit)(YFit(jval,i),i=0,NYterms-1)
ENDDO

IF(Numnu.ne.0)THEN
   ALLOCATE(DVR(0:numnu-1,0:numnu-1))
   ALLOCATE(Abscissas(0:numnu-1))
   ALLOCATE(Weights(0:numnu-1))
   ALLOCATE(PsiMax(0:numnu-1))
   ALLOCATE(ScalePsi(0:numnu-1))
   ALLOCATE(PsiDVR(0:numnu-1,0:numnu-1))
   lwork=100000
   ALLOCATE(work(100000))
   
   WRITE(DAF_DVR_Unit,*)
   WRITE(DAF_DVR_Unit,*)"Diagonal matrix elements of r"
   DO nu=0,numnu-1
      DO nup=0,numnu-1
         SUM0=Zero
         Sum1=Zero
         SUM2=Zero
         SUM3=Zero
         SUM4=Zero
         DO I=0,nmax
            SUM0=SUM0+psi(i,nu)*psi(i,nup)
            Sum1=Sum1+psi(i,nu)*psi(i,nup)*rvals(i)
            SUM2=SUM2+psi(i,nu)*psi(i,nup)*rvals(i)**2
            SUM3=SUM3+psi(i,nu)*psi(i,nup)*rvals(i)**3
            SUM4=SUM4+psi(i,nu)*psi(i,nup)*rvals(i)**4
         ENDDO
         sum0=sum0*hstep
         sum1=Sum1*hstep
         sum2=sum2*hstep
         sum3=sum3*hstep
         sum4=sum4*hstep
         IF(nu==nup)WRITE(DAF_DVR_Unit,'(2I5,5es14.7)')nu, nup, sum0, Sum1, sum2, sum3, sum4
         DVR(nu,nup)=Sum1
      ENDDO
   ENDDO
   DO nu=0,numnu-2
      LocPsiMax=MAXLOC(ABS(Psi(:,nu)),1)
      Fac=SIGN(One, Psi(LocPsiMax,nu) )
      !Psi(:,nu)=Fac*Psi(:,nu)
      PsiMax(nu)=MAXVAL(ABS(Psi(:,nu)))
      ScalePsi(nu)=(Eig(nu+1)-Eig(nu))/(Two*PsiMax(nu))
   ENDDO
   IF(numnu>1)ScalePsi(numnu-1)=ScalePsi(numnu-2)
   WRITE(DAF_DVR_Unit,*)"FBR Matrix"
   CALL MxOut_Unit(DVR,numnu,numnu,DAF_DVR_Unit)
   WRITE(DAF_DVR_Unit,*)"Calling dsyev"
   CALL dsyev(Jobz, UpLo, numnu, DVR, Numnu, Abscissas, work, lwork, info)
   WRITE(DAF_DVR_Unit,*)
   WRITE(DAF_DVR_Unit,*)"DVR Abscissas"
   WRITE(DAF_DVR_Unit,*)Abscissas
   
   WRITE(DAF_DVR_Unit,*)"FBR to DVR transformation matrix"
   !OPEN(Unit=9321,file='da.csv')
   DO i=0,numnu-1
      WRITE(DAF_DVR_Unit,'(300(ES13.5,","))')rvals(i),(DVR(i,nu),nu=0,numnu-1)
      !WRITE(9321,'(300(ES13.5,","))')rvals(i),(DVR(i,nu),nu=0,numnu-1)
   ENDDO
   
   DO nu=0,numnu-1
      CALL AKIMA (3,nmax+1, rvals, psi(:,nu), Numnu, Abscissas, PsiDVR(:,nu))
   ENDDO
   
   IF(jval==jmin)WRITE(PsiDAFG_Unit,'(2A13,300(A13,I5,","))')"i,","rvals(i),",("Psi",nu,nu=0,numnu-1)
   DO i=0,nmax
      WRITE(PsiDAFG_Unit,'(I13,",",300(ES13.5,","))')i,rvals(i),(Psi(i,nu),nu=0,numnu-1)
   ENDDO
   WRITE(PsiDAFG_Unit,*)"&"
   
    IF(jval==jmin)WRITE(PsiDAFGS_Unit,'(3A13,300(A13,I5,","))')"i,","rvals(i),","Potential,",("Psi",nu,nu=0,numnu-1)
   DO i=0,nmax
      WRITE(PsiDAFGS_Unit,'(I13,",",300(ES13.5,","))')i,rvals(i),DiatomicPotential(i,electronic),((Eig(nu)+ScalePsi(nu)*Psi(i,nu)),nu=0,numnu-1)
   ENDDO
   WRITE(PsiDAFGS_Unit,*)"&"
   
    IF(jval==jmin)THEN
      WRITE(EigPltDAF_Unit,'(2A12,300(A12,I5,","))')"jval,","rvals,",("Eig",nu,nu=0,numnu-1)
   ELSE
      WRITE(EigPltDAF_Unit,*)"&"
   ENDIF
   WRITE(EigPltDAF_Unit,'(I12," ,",300(ES12.4,","))')jval,rvals(0),(Eig(nu),nu=0,numnu-1)
   WRITE(EigPltDAF_Unit,'(I12," ,",300(ES12.4,","))')jval,rvals(nmax),(Eig(nu),nu=0,numnu-1)
      
   IF(jval==jmin)WRITE(PsiDAFD_Unit,'(2A13,300(A15,I5,","))')"i,","Abscissas(i),",("PsiDVR",nu,nu=0,numnu-1)
   DO i=0,numnu-1
      WRITE(PsiDAFD_Unit,'(I12," ,",300(ES15.5,","))')i,Abscissas(i),(PsiDVR(i,nu),nu=0,numnu-1)
   ENDDO
   WRITE(PsiDAFD_Unit,*)"&"
   
   IF(jval==jmin)WRITE(PsiDAFDS_Unit,'(2A13,300(A15,I5,","))')"i,","Abscissas(i),",("PsiDVR",nu,nu=0,numnu-1)
   DO i=0,numnu-1
      WRITE(PsiDAFDS_Unit,'(I12," ,",300(ES15.5,","))')i,Abscissas(i),((Eig(nu)+ScalePsi(nu)*PsiDVR(i,nu)),nu=0,numnu-1)
   ENDDO
   WRITE(PsiDAFDS_Unit,*)"&"
  
   Sum1=0.d0
   DO i=0,numnu-1
      Weights(i)=Zero
      DO nu=0,numnu-1
         Weights(i)=Weights(i)+PsiDVR(i,nu)**2
      ENDDO
      Weights(i)=one/Weights(i)
      Sum1=Sum1+Weights(i)
   ENDDO
   WRITE(DAF_DVR_Unit,*)
   WRITE(DAF_DVR_Unit,*)"DVR Weights"
   WRITE(DAF_DVR_Unit,*)Weights
   WRITE(DAF_DVR_Unit,*)"Sum1 of Weights=",Sum1
   
   WRITE(DAF_DVR_Unit,*)
   WRITE(DAF_DVR_Unit,*)"Weight and Abscissa Check"   
   DO nu=0,numnu-1
      DO nup=0,numnu-1
         SUM0=Zero
         Sum1=Zero
         SUM2=Zero
         SUM3=Zero
         SUM4=Zero
         DO I=0,numnu-1    
            SUM0=SUM0+Weights(i)*PsiDVR(i,nu)*PsiDVR(i,nup)
            Sum1=Sum1+Weights(i)*PsiDVR(i,nu)*PsiDVR(i,nup)*rvals(i)
            SUM2=SUM2+Weights(i)*PsiDVR(i,nu)*PsiDVR(i,nup)*rvals(i)**2
            SUM3=SUM3+Weights(i)*PsiDVR(i,nu)*PsiDVR(i,nup)*rvals(i)**3
            SUM4=SUM4+Weights(i)*PsiDVR(i,nu)*PsiDVR(i,nup)*rvals(i)**4
         ENDDO
         IF(nu==nup)WRITE(DAF_DVR_Unit,'(2I5,5es14.7)')nu, nup, sum0, Sum1, sum2, sum3, sum4
      ENDDO
   ENDDO
   WRITE(DAF_DVR_Unit,*)
   
   DEALLOCATE(DVR)
   DEALLOCATE(Abscissas)
   DEALLOCATE(Weights)
   DEALLOCATE(PsiMax)
   DEALLOCATE(ScalePsi)
   DEALLOCATE(PsiDVR)
   DEALLOCATE(work)
ENDIF

DEALLOCATE(YFit)
DEALLOCATE(Eig)
DEALLOCATE(Psi)
   
RETURN
END SUBROUTINE DAFToDVR