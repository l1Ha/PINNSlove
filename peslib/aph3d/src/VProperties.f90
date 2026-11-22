SUBROUTINE VProperties(NPoints, r, V_Pot, VInfinity, IArran, jval)
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
INTEGER NPoints, IPoint, IOP(2), IArran, jval, iter
REAL(Kind=WP_Kind) rAtVMin, rAtVMax, ydiff, Vmin, VMax, yintlast
REAL(Kind=WP_Kind) VInfinity, VAsymptoticErr, vlast(3), vnow(3), ylast, ynow, yinterp, vinterp(3)
REAL(Kind=WP_Kind) r(NPoints), V_Pot(NPoints), W(NPoints), A(NPoints), B(NPoints), C(NPoints)

WRITE(Out_Unit,'(" IArran=",I2,"  jval=",I5)')IArran, jval
VAsymptoticErr=VInfinity-V_Pot(NPoints)
WRITE(Out_Unit,'(10x,"RMax",16x,"V(RMax)",14x,"VInfinity",10x,"VInfinity-V(RMax)")')
WRITE(Out_Unit,'(4ES22.14)')r(NPoints),V_Pot(NPoints),VInfinity,VAsymptoticErr

VMax=Maxval(V_Pot)
VMin=Minval(V_Pot)
DO IPoint=1,Npoints
   IF(VMax==V_Pot(IPoint))rAtVMax=r(IPoint)
   IF(VMin==V_Pot(IPoint))rAtVMin=r(IPoint)
ENDDO
WRITE(Out_unit,'(" No Interpolation")')
WRITE(Out_Unit,'(" rAtVMax=",ES20.10," VMax=",ES20.10)')rAtVMax,VMax
WRITE(Out_Unit,'(" rAtVMin=",ES20.10," VMin=",ES20.10)')RAtVMin,VMin

WRITE(Out_Unit,*)
WRITE(Out_Unit,'(" Interpolated Results")')
IOP(1)=5
IOP(2)=5
CALL SPL1D1(NPoints, r, V_Pot, W, IOP, 1, A, B, C)

ylast=r(1)
CALL SPL1D2(NPoints, r, V_Pot, W, 1, ylast, vlast)
DO IPoint=2,NPoints
   ynow=r(IPoint)
   CALL SPL1D2(NPoints, r, V_Pot, W, 1, ynow, vnow)
   IF(ynow>=0.0d0)THEN   !tempmod
      IF((Vlast(1)-VInfinity)*(Vnow(1)-VInfinity)<0.d0)THEN
         ydiff=1.d0
         yintlast=-1.d0
         DO WHILE (ABS(ydiff)>1.d-12)
            yinterp=(Vlast(1)*ynow+VInfinity*ylast-VInfinity*ynow-Vnow(1)*ylast)/(Vlast(1)-Vnow(1))
            CALL SPL1D2(NPoints, r, V_Pot, W, 1, yinterp, vinterp)
            IF((Vinterp(1)-VInfinity)*(Vnow(1)-VInfinity)<0.d0)THEN
               ylast=yinterp
               vlast=vinterp
            ELSE
               ynow=yinterp
               vnow=vinterp
            ENDIF
            ydiff=yinterp-yintlast
            yintlast=yinterp
         ENDDO
         ylast=r(Ipoint-1)
         CALL SPL1D2(NPoints, r, V_Pot, W, 1, ylast, vlast)
         ynow=r(IPoint)
         CALL SPL1D2(NPoints, r, V_Pot, W, 1, ynow, vnow)
         WRITE(Out_Unit,*)"Potential crosses asymptotic limit at",yinterp
      ENDIF
      IF(Vlast(2)*Vnow(2)<0.d0)THEN
         ydiff=1.d0
         yintlast=-1.d0
         DO WHILE (ABS(ydiff)>1.d-12)
            yinterp=(Vlast(2)*ynow-Vnow(2)*ylast)/(Vlast(2)-Vnow(2))
            CALL SPL1D2(NPoints, r, V_Pot, W, 1, yinterp, vinterp)
            IF(Vinterp(2)*Vnow(2)<0.d0)THEN
               ylast=yinterp
               vlast=vinterp
            ELSE
               ynow=yinterp
               vnow=vinterp
            ENDIF
            ydiff=yinterp-yintlast
            yintlast=yinterp
         ENDDO
         ylast=r(Ipoint-1)
         CALL SPL1D2(NPoints, r, V_Pot, W, 1, ylast, vlast)
         ynow=r(IPoint)
         CALL SPL1D2(NPoints, r, V_Pot, W, 1, ynow, vnow)
         IF(vnow(2)>0.d0)THEN
            WRITE(Out_Unit,'(" Minimum of potential at",ES21.14)')yinterp
            WRITE(Out_Unit,'("     Value of the potential at minimum=",ES21.14)')vinterp(1)
            WRITE(Out_Unit,'("     Curvature of the potential at minimum=",ES21.14)')vinterp(3)
         ELSE
            WRITE(Out_Unit,*)" Maximum of potential at",yinterp
            WRITE(Out_Unit,'("     Value of the potential at maximum=",ES21.14)')vinterp(1)
            WRITE(Out_Unit,'("     Curvature of the potential at maximum=",ES21.14)')vinterp(3)
         ENDIF
      ENDIF
      IF(Vlast(3)*Vnow(3)<0.d0)THEN
         ydiff=1.d0
         yintlast=-1.d0
         DO WHILE (ABS(ydiff)>1.d-12)
            yinterp=(Vlast(3)*ynow-Vnow(3)*ylast)/(Vlast(3)-Vnow(3))
            CALL SPL1D2(NPoints, r, V_Pot, W, 1, yinterp, vinterp)
            IF(Vinterp(3)*Vnow(3)<0.d0)THEN
               ylast=yinterp
               vlast=vinterp
            ELSE
               ynow=yinterp
               vnow=vinterp
            ENDIF
            ydiff=yinterp-yintlast
            yintlast=yinterp
         ENDDO
         ylast=r(Ipoint-1)
         CALL SPL1D2(NPoints, r, V_Pot, W, 1, ylast, vlast)
         ynow=r(IPoint)
         CALL SPL1D2(NPoints, r, V_Pot, W, 1, ynow, vnow)
         WRITE(Out_Unit,*)"Inflection point at",yinterp
      ENDIF
      WRITE(VEE_Unit,'(3(es16.8,","),es16.8)')ynow,vnow
   ENDIF
   ylast=ynow
   vlast=vnow
ENDDO

RETURN
END
