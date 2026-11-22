SUBROUTINE Graf_2d(NPoints, NCurves, XVals, Curves, Xmin, Xmax, Title, SubTitle, XLabel, YLabel, MTVFile)
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
INTEGER IPoint, JCurve, ltype, lcolr
INTEGER, PARAMETER::  LineWidth=4
INTEGER, INTENT(IN) :: NPoints, NCurves
REAL(Kind=WP_Kind), INTENT(IN):: Xmin, Xmax, XVals(NPoints), Curves(NCurves,NPoints)
CHARACTER(LEN=*), INTENT(IN):: Title, SubTitle, XLabel, YLabel, MTVFile
OPEN(Unit=MTV_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//MTVFile,form='formatted',status='unknown')
WRITE(MTV_Unit,*) '$ DATA=CURVE2D'
WRITE(MTV_Unit,'(a12,A,a2)') '%toplabel= "',TRIM(title),' "'
WRITE(MTV_Unit,'(a12,A,a2)') '%subtitle= "',TRIM(subtitle),' "'
WRITE(MTV_Unit,'(a12,A,a2)') '%xlabel  = "',TRIM(XLabel),' "'
WRITE(MTV_Unit,'(a12,A,a2)') '%ylabel  = "',TRIM(YLabel),' "'
WRITE(MTV_Unit,'("% xmin =",ES15.7)')xmin
WRITE(MTV_Unit,'("% xmax =",ES15.7)')xmax
WRITE(MTV_Unit,'("# First Curve")')
WRITE(MTV_Unit,'("% linelabel = C1")')
WRITE(MTV_Unit,'("% linetype  = 1")')
WRITE(MTV_Unit,'("% linecolor  =",I1,A2)')1,' "'
!WRITE(MTV_Unit,'("% markertype =",I2,A2)')1,' "'
WRITE(MTV_Unit,'("% linewidth  =",I2,A2)')LineWidth,' "'
DO IPoint=1,NPoints
   WRITE(MTV_Unit,'(1f10.6,e15.5)') XVals(IPoint),Curves(1,IPoint)
ENDDO
DO JCurve=2,NCurves
   WRITE(MTV_Unit,'("# Curve ",I5,A2)')JCurve,'"'
   IF(JCurve<10)THEN
      WRITE(MTV_Unit,'("% linelabel = C",I1,A2 )')JCurve,'"'
   ELSEIF(JCurve<=15)THEN
      WRITE(MTV_Unit,'("% linelabel = C",I2,A2 )')JCurve,'"'
   ENDIF
   ltype=MOD(JCurve,9 )+1
   lcolr=MOD(JCurve,7 )+1
   !mtype=MOD(JCurve,13)+1
   WRITE(MTV_Unit,'("% linetype   =",I1,A2)')ltype,' "'
   WRITE(MTV_Unit,'("% linecolor  =",I1,A2)')lcolr,' "'
   !WRITE(MTV_Unit,'("% markertype =",I2,A2)')mtype,' "'
   WRITE(MTV_Unit,'("% linewidth  =",I2,A2)')LineWidth,' "'
   DO IPoint=1,NPoints
      WRITE(MTV_Unit,'(1f10.6,e15.5)') XVals(IPoint),Curves(JCurve,IPoint)
   ENDDO
ENDDO

WRITE(MTV_Unit,*) '$END'
CLOSE(MTV_Unit)

RETURN
END
