SUBROUTINE Boundaries(RhoMin, RhoMax, DeltaRho1, DeltaRho2, NRho)
USE Numeric_Kinds_Module
USE FileUnits_Module
Use Boundary_Module

IMPLICIT NONE
SAVE
REAL(Kind=WP_Kind) RhoMin, RhoMax, DeltaRho1, DeltaRho2, Rho, RhoLast, RhoNext
INTEGER IthRho, NRho, IFlag, calcparm
LOGICAL, Parameter:: DEBUG=.FALSE.

WRITE(Out_Unit,*)
WRITE(Out_Unit,'("Boundaries: Start")')
WRITE(Out_Unit,*)'LastRho=',LastRho

!
! See IF NRho > Max_Sector
!
IF(Max_Sector<NRho)THEN
   WRITE(Out_Unit,'("Error NRho is larger than Max_Sector",2I10)')NRho, Max_Sector
   STOP 'Boundaries'
ENDIF

!
! Determine DeltaRho2
!
IFlag=CalcParm(NRho, RhoMin, DeltaRho1, DeltaRho2, RhoMax, Kase_CalcParm)
IF(IFlag/=0)THEN
   WRITE(Out_Unit,'("Error Boundaries: IFlag=",L2)')IFlag
   WRITE(Msg_Unit,'("Error Boundaries: IFlag=",L2)')IFlag
   STOP 'Boundaries: IFlag/=0'
ENDIF

!
! Determine where Surface Functions will be calculated
!
DO IthRho=1, NRho
   CALL RhoCalc(RhoMin, DeltaRho1, DeltaRho2, IthRho, Rho, RhoLast, RhoNext)
   Basis_Dist(IthRho+LastRho)=Rho
ENDDO

!
! Determine Sector Boundaries
!
Sector_Boundaries(1)=Max(0.9d0*Basis_Dist(1),2*Basis_Dist(1)-Basis_Dist(2))
IF(NRho>1)THEN
   DO IthRho=2+LastRho, NRho+LastRho
      Sector_Boundaries(IthRho)=(Basis_Dist(IthRho-1)+Basis_Dist(IthRho))/2.d0
   ENDDO
ENDIF
IF(NRho<LastRho)THEN
   Sector_Boundaries(NRho+1+LastRho)=Basis_Dist(NRho+LastRho) 
ELSE
   Sector_Boundaries(NRho+1+LastRho)=2*Basis_Dist(NRho+LastRho)-Basis_Dist(NRho+LastRho-1)
ENDIF

LastRho=LastRho+NRho

IF(DEBUG)THEN
   WRITE(Out_Unit,'("Start: DEBUG=.TRUE. in Boundaries")')
   DO IthRho=1, NRho
      WRITE(Out_Unit,10)IthRho,Sector_Boundaries(IthRho),Basis_Dist(IthRho),Sector_Boundaries(IthRho+1)
   ENDDO
   WRITE(Out_Unit,*)'END:   DEBUG=.TRUE. in Boundaries'
ENDIF


WRITE(Out_Unit,'("Boundaries: END")')
WRITE(Out_Unit,*)

RETURN
  10  FORMAT(5x,'IthRho=',i5,' StartOfSector=',  1pe21.14,' BasisEvaluation=',1pe21.14,' EndOfSector=',    1pe21.14)
ENDSUBROUTINE Boundaries
