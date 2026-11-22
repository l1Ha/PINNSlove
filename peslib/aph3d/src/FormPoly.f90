SUBROUTINE FormPoly (BeginVec, EndVec, NTheta, NChi)
USE Numeric_Kinds_Module
USE FileUnits_Module
INTEGER Ntheta, Nchi, Ichi, Itheta
REAL(KIND=WP_Kind) BeginVec(3,NTheta-1,NChi-1), EndVec(3,NTheta-1,NChi-1)

OPEN(Unit=p_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'PlotVee/FormPoly_p.txt')
OPEN(Unit=q_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'PlotVee/FormPoly_q.txt')
WRITE(q_unit,19)

DO IChi = 1, NChi-2
   WRITE(p_unit,*)'    3'
   WRITE(p_unit,9)BeginVec(1,1,IChi),BeginVec(2,1,IChi),BeginVec(3,1,IChi),   &
   EndVec(1,1,IChi),EndVec(2,1,IChi),EndVec(3,1,IChi),EndVec(1,1,IChi+1),EndVec(2,1,IChi+1),EndVec(3,1,IChi+1)
   WRITE(q_unit,20)BeginVec(1,1,IChi),BeginVec(2,1,IChi),BeginVec(3,1,IChi),   &
   EndVec(1,1,IChi),EndVec(2,1,IChi),EndVec(3,1,IChi),EndVec(1,1,IChi+1),EndVec(2,1,IChi+1),EndVec(3,1,IChi+1)

ENDDO
WRITE(p_unit,*)'    3'
WRITE(p_unit,9)BeginVec(1,1,NChi-1),BeginVec(2,1,NChi-1),BeginVec(3,1,NChi-1), &
EndVec(1,1,NChi-1),EndVec(2,1,NChi-1),EndVec(3,1,NChi-1),EndVec(1,1,1),EndVec(2,1,1),EndVec(3,1,1)
WRITE(q_unit,20)BeginVec(1,1,NChi-1),BeginVec(2,1,NChi-1),BeginVec(3,1,NChi-1), &
EndVec(1,1,NChi-1),EndVec(2,1,NChi-1),EndVec(3,1,NChi-1),EndVec(1,1,1),EndVec(2,1,1),EndVec(3,1,1)

DO IChi = 1, NChi-2
   DO ITheta = 2, NTheta-1
   WRITE(p_unit,*)'    4'
   WRITE(p_unit,8)BeginVec(1,ITheta,IChi),BeginVec(2,ITheta,IChi),BeginVec(3,ITheta,IChi),EndVec(1,ITheta,IChi),EndVec(2,ITheta,IChi),&
                  EndVec(3,ITheta,IChi),EndVec(1,ITheta,IChi+1),EndVec(2,ITheta,IChi+1),EndVec(3,ITheta,IChi+1),BeginVec(1,ITheta,IChi+1),&
                  BeginVec(2,ITheta,IChi+1),BeginVec(3,ITheta,IChi+1)
   WRITE(q_unit,22)BeginVec(1,ITheta,IChi),BeginVec(2,ITheta,IChi),BeginVec(3,ITheta,IChi),EndVec(1,ITheta,IChi),EndVec(2,ITheta,IChi),&
                  EndVec(3,ITheta,IChi),EndVec(1,ITheta,IChi+1),EndVec(2,ITheta,IChi+1),EndVec(3,ITheta,IChi+1),BeginVec(1,ITheta,IChi+1),&
                  BeginVec(2,ITheta,IChi+1),BeginVec(3,ITheta,IChi+1)
   ENDDO
ENDDO

DO ITheta = 2, NTheta-2
  WRITE(p_unit,*)'    4'
  WRITE(p_unit,8)BeginVec(1,ITheta,NChi-1),BeginVec(2,ITheta,NChi-1),BeginVec(3,ITheta,NChi-1),EndVec(1,ITheta,NChi-1),&
                  EndVec(2,ITheta,NChi-1),EndVec(3,ITheta,NChi-1),EndVec(1,ITheta,1),EndVec(2,ITheta,1),EndVec(3,ITheta,1),&
                  BeginVec(1,ITheta,1),BeginVec(2,ITheta,1),BeginVec(3,ITheta,1)
  WRITE(q_unit,22)BeginVec(1,ITheta,NChi-1),BeginVec(2,ITheta,NChi-1),BeginVec(3,ITheta,NChi-1),EndVec(1,ITheta,NChi-1),&
                  EndVec(2,ITheta,NChi-1),EndVec(3,ITheta,NChi-1),EndVec(1,ITheta,1),EndVec(2,ITheta,1),EndVec(3,ITheta,1),&
                  BeginVec(1,ITheta,1),BeginVec(2,ITheta,1),BeginVec(3,ITheta,1)
ENDDO

IF (NTheta>2)THEN
 DO ITheta = NTheta-1, NTheta-1
  WRITE(p_unit,*)'    4'
  WRITE(p_unit,8)BeginVec(1,ITheta,NChi-1),BeginVec(2,ITheta,NChi-1),BeginVec(3,ITheta,NChi-1),EndVec(1,ITheta,NChi-1),&
                  EndVec(2,ITheta,NChi-1),EndVec(3,ITheta,NChi-1),EndVec(1,ITheta,1),EndVec(2,ITheta,1),EndVec(3,ITheta,1),&
                  BeginVec(1,ITheta,1),BeginVec(2,ITheta,1),BeginVec(3,ITheta,1)
  WRITE(q_unit,23)BeginVec(1,ITheta,NChi-1),BeginVec(2,ITheta,NChi-1),BeginVec(3,ITheta,NChi-1),EndVec(1,ITheta,NChi-1),&
                  EndVec(2,ITheta,NChi-1),EndVec(3,ITheta,NChi-1),EndVec(1,ITheta,1),EndVec(2,ITheta,1),EndVec(3,ITheta,1),&
                  BeginVec(1,ITheta,1),BeginVec(2,ITheta,1),BeginVec(3,ITheta,1)
 ENDDO
ENDIF

WRITE(q_unit,21)

   9  FORMAT(1x,3f12.5,/1x,3f12.5,/1x,3f12.5)
   8  FORMAT(1x,3f12.5,/1x,3f12.5,/1x,3f12.5,/1x,3f12.5)
  19  FORMAT('Graphics3D[{')
  20  FORMAT(5x,'Polygon[{','{',f10.5,',',f10.5,',',f10.5,'},  ',/10x,&
  '{',f10.5,',',f10.5,  ',',f10.5,'},', /10x,'{',f10.5,',',f10.5,',',f10.5,'}','}],' )
  21  FORMAT(10x,'}]')
  22  FORMAT(5x,'Polygon[{','{',f10.5,',',f10.5,',',f10.5,'},  ',/10x,&
  '{',f10.5,',',f10.5,',',f10.5,'},',/10x,'{',f10.5,',',f10.5,',',f10.5,'},',/10x,'{',f10.5,',',f10.5,',',f10.5,'}','}],' )
  23  FORMAT(5x,'Polygon[{','{',f10.5,',',f10.5,',',f10.5,'},  ',/10x,&
  '{',f10.5,',',f10.5,',',f10.5,'},',/10x,'{',f10.5,',',f10.5,',',f10.5,'},',/10x,'{',f10.5,',',f10.5,',',f10.5,'}','}]' )

RETURN
END
