SUBROUTINE VecDef(ThetaVal,NTheta,ChiVals,NChi,Pot,valmax)
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
INTEGER pic_unit, IChi, Itheta, IthVec, NChi, NTheta, NVecs, NWords
REAL(KIND=WP_Kind) BeginVec(3,512000), EndVec(3,512000)
REAL(KIND=WP_Kind) ThetaVal(NTheta), ChiVals(NChi), Pot(NTheta,NChi),           &
                 & VecList(800000), valmx
REAL(KIND=WP_Kind) valmax
 
 
pic_unit=92
OPEN(Unit=pic_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Pictre.bin',form='unformatted')
valmx = valmax
WRITE(Msg_Unit,*)
WRITE(Msg_Unit,*)"In VecDef"
WRITE(Msg_Unit,*)"Ntheta=", Ntheta," Nchi=", Nchi, " Valmx=", Valmx, " Thetaval(1)=", Thetaval(1), "Chivals(1)=", chivals(1)
WRITE(pic_unit)ntheta, nchi, valmx
DO ichi = 1, nchi
   DO itheta = 1, ntheta
       WRITE(pic_unit)itheta, thetaval(itheta), ichi, chivals(ichi),            &
                    & pot(itheta,ichi)
   ENDDO
ENDDO
CLOSE(UNIT=pic_unit)
 
IthVec = 0
DO Ichi = 1, NChi - 1
   DO ITheta = 1, NTheta - 1
       IthVec = IthVec + 1
       BeginVec(1,IthVec) = sin(ThetaVal(ITheta))*cos(ChiVals(IChi))            &
                          & *Pot(Itheta,IChi)
       BeginVec(2,IthVec) = sin(ThetaVal(ITheta))*sin(ChiVals(IChi))            &
                          & *Pot(Itheta,IChi)
       BeginVec(3,IthVec) = Pot(ITheta,IChi)*cos(ThetaVal(ITheta))
       EndVec(1,IthVec) = sin(ThetaVal(ITheta+1))*cos(ChiVals(IChi))            &
                        & *Pot(Itheta+1,IChi)
       EndVec(2,IthVec) = sin(ThetaVal(ITheta+1))*sin(ChiVals(IChi))            &
                        & *Pot(Itheta+1,IChi)
       EndVec(3,IthVec) = Pot(ITheta+1,IChi)*cos(ThetaVal(ITheta+1))
   ENDDO
ENDDO
 
DO ITheta = 2, NTheta
   DO IChi = 1, NChi - 1
       IthVec = IthVec + 1
       BeginVec(1,IthVec) = sin(ThetaVal(ITheta))*cos(ChiVals(IChi))            &
                          & *Pot(Itheta,IChi)
       BeginVec(2,IthVec) = sin(ThetaVal(ITheta))*sin(ChiVals(IChi))            &
                          & *Pot(Itheta,IChi)
       BeginVec(3,IthVec) = Pot(ITheta,IChi)*cos(ThetaVal(ITheta))
       EndVec(1,IthVec) = sin(ThetaVal(ITheta))*cos(ChiVals(IChi+1))            &
                        & *Pot(Itheta,IChi+1)
       EndVec(2,IthVec) = sin(ThetaVal(ITheta))*sin(ChiVals(IChi+1))            &
                        & *Pot(Itheta,IChi+1)
       EndVec(3,IthVec) = Pot(ITheta,IChi+1)*cos(ThetaVal(ITheta))
   ENDDO
ENDDO
 
NVecs = IthVec
NWords = 8*NVecs
WRITE(Out_Unit,*)'NVecs,Nwords', NVecs, NWords
WRITE(Out_Unit,*)'Calling FormVec'
CALL FormVec(BeginVec,EndVec,NVecs,VecList,NWords)
WRITE(Out_Unit,*)'Calling FormElem'
CALL FormElem(BeginVec,EndVec,NTheta,NChi,valmax)
WRITE(Out_Unit,*)'Calling FormPoly'
CALL FormPoly(BeginVec,EndVec,NTheta,NChi)
WRITE(Out_Unit,*)'Calling FormTriangle'
!CALL FormTriangle(ThetaVal,NTheta,ChiVals,NChi,Pot)
WRITE(Out_Unit,*)'Completed writing DATA'

RETURN
ENDSUBROUTINE VECDEF
