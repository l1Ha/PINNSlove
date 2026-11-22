SUBROUTINE SFPlots(Nang, ThetaVal, ChiVals, PhiPlt, SFType)
USE Numeric_Kinds_Module
USE Numbers_Module
USE FileUnits_Module
IMPLICIT NONE
CHARACTER(LEN=6) SFType
INTEGER Iang, Nang
REAL(KIND=SP_Kind) SP_PSI     ! Single Precision for Plotting
REAL(KIND=WP_Kind) valmax, SignAtMax, theta, chi, phi, u, v
REAL(KIND=WP_Kind) ThetaVal(Nang), ChiVals(Nang), PhiPlt(nang)

OPEN(Unit=Plt_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'SfunctionPlots/SfuncType_'//SFType//'.txt')
ValMax = 0.D0
SignAtMax = 0.D0
DO Iang = 1, Nang
   phi = PhiPlt(Iang)
   IF(ValMax<ABS(Phi))THEN
      ValMax=ABS(Phi)
      SignATMax=SIGN(1.d0,Phi)
   ENDIF
ENDDO

WRITE(Out_Unit,*)"ValMax,SignAtMax=",ValMax,SignAtMax
WRITE(Msg_Unit,*)"ValMax,SignAtMax=",ValMax,SignAtMax

DO Iang = 1, Nang
   theta = ThetaVal(Iang)
   chi = ChiVals(Iang)
   phi = PhiPlt(Iang)
   CALL PMap(1, 1, Zero, theta, chi, u, v)
   SP_PSI=phi*SignAtMax/ValMax
   IF(ABS(SP_PSI)<EPSILON(1.0_SP_Kind))SP_PSI=0.0
   WRITE(Plt_Unit,'(3ES12.4)')u, v, SP_Psi
ENDDO
CLOSE(Plt_Unit)

RETURN
END