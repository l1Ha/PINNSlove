SUBROUTINE CollinearFixedRho(Basis_Distxx,Nsectorsxx)
USE Numeric_Kinds_Module
USE Numbers_Module
USE Parms_Module
USE FileUnits_Module
USE Masses_Module
USE MassFactor2_Module
USE Narran_Module
USE BasisRegions_Module

IMPLICIT NONE
INTEGER, PARAMETER :: nc=181
INTEGER k, ic, isector, nsectorsxx
REAL(KIND=WP_Kind) pot, rho, theta, xchi, sin, cos, dc, Basis_Distxx(21)
REAL(KIND=WP_Kind), ALLOCATABLE :: PotChi(:)
INTRINSIC sqrt, sin, cos
EXTERNAL surface

OPEN(Unit=639,File=TRIM(OutDIR)//'PlotVee/Collinear/CollinearFixedRho.csv',Form='FORMATTED')
WRITE(639,'(A,",",3000(A4,1PE10.2,","))')"xchi", ("Rho=",Basis_Dist(ISector),ISector=1,NSectors)
theta=pi/two
dc=Pi/(nc-1)
ALLOCATE(PotChi(1:NSectors))
!-----------------------------------------------------------------------
!  Calculate the interparticle distances and store them in the array r.
!-----------------------------------------------------------------------
DO ic=0,nc-1
  DO ISector=1,NSectors
    rho=Basis_Dist(ISector)
    xchi=ic*dc
    chi3(1)=xchi
    chi3(2)=chi3(1)-chij(2,1)
    chi3(3)=chi3(1)-chij(3,1)
    DO k=1,narran
       s(k)=rho*sqrt(Half-Half*sin(theta)*cos(Two*chi3(k)))
       r(k)=d(k)*s(k)
    ENDDO

!-----------------------------------------------------------------------
! This CALL returns the interaction potential in Hartree atomic units.
!-----------------------------------------------------------------------
    CALL surface (pot, r)
    potchi(ISector)=pot
  ENDDO
  WRITE(639,'(3000(1PE14.6,","))')xchi*(180.d0/Pi), (potchi(ISector),ISector=1,NSectors)
ENDDO
DEALLOCATE(PotChi)
CLOSE(Unit=639)

RETURN
ENDSUBROUTINE CollinearFixedRho