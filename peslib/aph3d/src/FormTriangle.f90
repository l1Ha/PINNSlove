SUBROUTINE FormTriangle(ThetaVal,NTheta,ChiVals,NChi,Pot)
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
INTEGER NTheta, NChi, ITheta, IChi, NTriangle
REAL(KIND=WP_Kind) Pot(NTheta,NChi), ThetaVal(NTheta), ChiVals(NChi)
REAL(KIND=WP_Kind) Norm(3), r, Theta, Chi, x, y, z

OPEN(Unit=Norm_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'PlotVee/PlotNormal.txt')
!cccc      CALL open_spinsurface(2*(NChi-1)*(NTheta-1)-(NChi-1),'plot.binary')
NTriangle = 0
WRITE(Norm_Unit,*)'Number_Triangles = ', 2*(NChi-1)*(NTheta-1) - (NChi-1)
DO ITheta = 1, NTheta - 1
   DO IChi = 1, NChi - 1
       NTriangle = NTriangle + 1
       WRITE(Norm_Unit,*)'Triangle=', NTriangle
! First Triangle First Point
       r = Pot(ITheta,IChi)
       Theta = ThetaVal(ITheta)
       Chi = ChiVals(IChi)
       CALL Coord(x,y,z,Theta,Chi,r)
       CALL Normal(Pot,NTheta,NChi,ITheta,IChi,Norm,ThetaVal,ChiVals)
       WRITE(Norm_Unit,99001)x, y, z, Norm(1), Norm(2), Norm(3)
!cccc            CALL WRITE_vertex(x,y,z,Norm(1),Norm(2),Norm(3))
 
! First Triangle Second Point
       r = Pot(ITheta+1,IChi)
       Theta = ThetaVal(ITheta+1)
       Chi = ChiVals(IChi)
       CALL Coord(x,y,z,Theta,Chi,r)
       CALL Normal(Pot,NTheta,NChi,ITheta+1,IChi,Norm,ThetaVal,ChiVals)
       WRITE(Norm_Unit,99001)x, y, z, Norm(1), Norm(2), Norm(3)
!cccc            CALL WRITE_vertex(x,y,z,Norm(1),Norm(2),Norm(3))
 
! First Triangle Third Point
       r = Pot(ITheta+1,IChi+1)
       Theta = ThetaVal(ITheta+1)
       Chi = ChiVals(IChi+1)
       CALL Coord(x,y,z,Theta,Chi,r)
       CALL Normal(Pot,NTheta,NChi,ITheta+1,IChi+1,Norm,ThetaVal,ChiVals)
       WRITE(Norm_Unit,99001)x, y, z, Norm(1), Norm(2), Norm(3)
!cccc            CALL WRITE_vertex(x,y,z,Norm(1),Norm(2),Norm(3))
 
! Second Triangle First Point
       IF(ITheta/=1)THEN
           NTriangle = NTriangle + 1
           WRITE(Norm_Unit,*)'Triangle=', NTriangle
 
           r = Pot(ITheta,IChi)
           Theta = ThetaVal(ITheta)
           Chi = ChiVals(IChi)
           CALL Coord(x,y,z,Theta,Chi,r)
           CALL Normal(Pot,NTheta,NChi,ITheta,IChi,Norm,ThetaVal,ChiVals)
           WRITE(Norm_Unit,99001)x, y, z, Norm(1), Norm(2), Norm(3)
!cccc               CALL WRITE_vertex(x,y,z,Norm(1),Norm(2),Norm(3))
 
! Second Triangle Second Point
           r = Pot(ITheta+1,IChi+1)
           Theta = ThetaVal(ITheta+1)
           Chi = ChiVals(IChi+1)
           CALL Coord(x,y,z,Theta,Chi,r)
           CALL Normal(Pot,NTheta,NChi,ITheta+1,IChi+1,Norm,ThetaVal,ChiVals)
           WRITE(Norm_Unit,99001)x, y, z, Norm(1), Norm(2), Norm(3)
!cccc               CALL WRITE_vertex(x,y,z,Norm(1),Norm(2),Norm(3))
 
! Second Triangle Third Point
           r = Pot(ITheta,IChi+1)
           Theta = ThetaVal(ITheta)
           Chi = ChiVals(IChi+1)
           CALL Coord(x,y,z,Theta,Chi,r)
           CALL Normal(Pot,NTheta,NChi,ITheta,IChi+1,Norm,ThetaVal,ChiVals)
           WRITE(Norm_Unit,99001)x, y, z, Norm(1), Norm(2), Norm(3)
!cccc               CALL WRITE_vertex(x,y,z,Norm(1),Norm(2),Norm(3))
       ENDIF
 
   ENDDO
ENDDO
 
!cccc      CALL close_spinsurface()
CLOSE(UNIT=Norm_Unit)

99001 FORMAT(1x,6F12.5)
ENDSUBROUTINE FORMTRIANGLE
