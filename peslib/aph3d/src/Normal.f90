!*==NORMAL.f90  processed by SPAG 6.55Dc at 15:05 on 21 Aug 2007
SUBROUTINE Normal(Pot,NTheta,NChi,ITheta,IChi,Norm,ThetaVal,ChiVals)    !      0
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
!*--NORMAL6
INTEGER Ntheta, Nchi, ITheta, IChi
REAL(KIND=WP_Kind) Theta, Chi, t1, t2, t3, c1, c2, c3, VTheta, VChi
REAL(KIND=WP_Kind) ThetaVal(NTheta), ChiVals(NChi)
REAL(KIND=WP_Kind) Pot(NTheta,NChi), r
REAL(KIND=WP_Kind) PotVal(3,3)
REAL(KIND=WP_Kind) Norm(3)
 
Theta = ThetaVal(ITheta)                                                !      0
Chi = ChiVals(IChi)
r = Pot(ITheta,IChi)
 
IF(ITheta==1)THEN
   Norm(1) = 0.0                                                        !      0
   Norm(2) = 0.0
   Norm(3) = 1.0
 
ELSEIF(ITheta==NTheta)THEN
   T1 = ThetaVal(ITheta-1)                                              !      0
   T2 = ThetaVal(Itheta)
   T3 = 2*ThetaVal(Itheta) - ThetaVal(Itheta-1)
 
   IF(IChi/=1)THEN
       C1 = ChiVals(IChi-1)                                             !      0
       PotVal(1,1) = Pot(Itheta-1,IChi-1)
       PotVal(2,1) = Pot(Itheta,IChi-1)
       PotVal(3,1) = Pot(Itheta-1,IChi-1)
   ELSE
       C1 = -ChiVals(2)                                                 !      0
       PotVal(1,1) = Pot(Itheta-1,NChi-1)
       PotVal(2,1) = Pot(Itheta,NChi-1)
       PotVal(3,1) = Pot(Itheta-1,NChi-1)
   ENDIF
 
   C2 = ChiVals(IChi)                                                   !      0
   PotVal(1,2) = Pot(Itheta-1,IChi)
   PotVal(2,2) = Pot(Itheta,IChi)
   PotVal(3,2) = Pot(Itheta-1,IChi)
 
   IF(IChi/=NChi)THEN
       C3 = ChiVals(IChi+1)                                             !      0
       PotVal(1,3) = Pot(Itheta-1,IChi+1)
       PotVal(2,3) = Pot(Itheta,IChi+1)
       PotVal(3,3) = Pot(Itheta-1,IChi+1)
   ELSE
       C3 = 2*ChiVals(IChi) - ChiVals(IChi-1)                           !      0
       PotVal(1,3) = Pot(Itheta-1,2)
       PotVal(2,3) = Pot(Itheta,2)
       PotVal(3,3) = Pot(Itheta-1,2)
   ENDIF
 
   CALL deriv(T1,T2,T3,C1,C2,C3,PotVal,VTheta,VChi)                     !      0
   VChi = VChi/r
   Norm(1) = cos(Chi) + VChi*sin(Chi)
   Norm(2) = sin(Chi) - VChi*cos(Chi)
   Norm(3) = 0.0
 
ELSE
   T1 = ThetaVal(ITheta-1)                                              !      0
   T2 = ThetaVal(Itheta)
   T3 = ThetaVal(Itheta+1)
   IF(IChi/=1)THEN
       C1 = ChiVals(IChi-1)                                             !      0
       PotVal(1,1) = Pot(Itheta-1,IChi-1)
       PotVal(2,1) = Pot(Itheta,IChi-1)
       PotVal(3,1) = Pot(Itheta+1,IChi-1)
   ELSE
       C1 = -ChiVals(2)                                                 !      0
       PotVal(1,1) = Pot(Itheta-1,NChi-1)
       PotVal(2,1) = Pot(Itheta,NChi-1)
       PotVal(3,1) = Pot(Itheta+1,NChi-1)
   ENDIF
 
   C2 = ChiVals(IChi)                                                   !      0
   PotVal(1,2) = Pot(Itheta-1,IChi)
   PotVal(2,2) = Pot(Itheta,IChi)
   PotVal(3,2) = Pot(Itheta+1,IChi)
 
   IF(IChi/=Nchi)THEN
       C3 = ChiVals(IChi+1)                                             !      0
       PotVal(1,3) = Pot(Itheta-1,IChi+1)
       PotVal(2,3) = Pot(Itheta,IChi+1)
       PotVal(3,3) = Pot(Itheta+1,IChi+1)
   ELSE
       C3 = 2*ChiVals(IChi) - ChiVals(IChi-1)                           !      0
       PotVal(1,3) = Pot(Itheta-1,2)
       PotVal(2,3) = Pot(Itheta,2)
       PotVal(3,3) = Pot(Itheta+1,2)
   ENDIF
 
   CALL deriv(T1,T2,T3,C1,C2,C3,PotVal,VTheta,VChi)                     !      0
   VTheta = VTheta/r
   VChi = VChi/(r*sin(Theta))
   Norm(1) = sin(Theta)*cos(Chi) - VTheta*cos(Theta)*cos(Chi) + VChi*sin(Chi)
   Norm(2) = sin(Theta)*sin(Chi) - VTheta*cos(Theta)*sin(Chi) - VChi*cos(Chi)
   Norm(3) = cos(Theta) + VTheta*sin(Theta)
 
ENDIF
!
! Normalize the Normal vector.
!
r = (Norm(1)**2+Norm(2)**2+Norm(3)**2)**0.5                             !      0
Norm(1) = Norm(1)/r
Norm(2) = Norm(2)/r
Norm(3) = Norm(3)/r
 
ENDSUBROUTINE NORMAL
