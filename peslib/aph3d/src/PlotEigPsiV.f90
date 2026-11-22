SUBROUTINE PlotEigPsiV(Dist,NGrid,Psi,Eig,Neig,PltPsi_Unit,PsiFileName,PltEig_Unit,EigFileName)
USE Numeric_Kinds_Module
USE Numbers_Module
USE FileUnits_Common_Module
USE FileUnits_OneDim_Module
USE OneDim_Parms_Module
USE QstateDAF_Module
IMPLICIT NONE
CHARACTER(LEN=*)PsiFileName, EigFileName
INTEGER NGrid, Neig, i, LocPsiMax, PltPsi_Unit, PltEig_Unit
REAL(KIND=WP_Kind) :: Dist(0:NGrid), Psi(0:NGrid,0:Neig), Eig(0:Neig)
REAL(KIND=WP_Kind), ALLOCATABLE:: PsiMax(:), ScalePsi(:)

ALLOCATE(PsiMax(0:NEig))
ALLOCATE(ScalePsi(0:Neig))
DO nu=0,Neig-1
   LocPsiMax=MAXLOC(ABS(Psi(:,nu)),1)
   PsiMax(nu)=MAXVAL(ABS(Psi(:,nu)))
   ScalePsi(nu)=(Eig(nu+1)-Eig(nu))/(Two*PsiMax(nu))
ENDDO

OPEN(Unit=PltPsi_Unit,File=TRIM(OutDIR)//"DiatomicPlots/"//PsiFileName)
DO i=0,nmax
   WRITE(PltPsi_Unit,'(300ES12.4)')Dist(i),((Eig(nu)+ScalePsi(nu)*Psi(i,nu)),nu=0,Neig-1)
ENDDO
CLOSE(Unit=PltPsi_Unit)

OPEN(Unit=PltEig_Unit,FILE=TRIM(OutDIR)//"DiatomicPlots/"//EigFileName)
WRITE(PltEig_Unit,'(300ES12.4)')Dist(0),    (Eig(nu),nu=0,Neig-1)
WRITE(PltEig_Unit,'(300ES12.4)')Dist(NGrid),(Eig(nu),nu=0,Neig-1)
CLOSE(Unit=PltEig_Unit)
   
DEALLOCATE(PsiMax)
DEALLOCATE(ScalePsi)
   
RETURN
END SUBROUTINE PlotEigPsiV