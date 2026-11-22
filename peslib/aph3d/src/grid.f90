SUBROUTINE Grid(Rho, S, NS, CapTheta, NT, Chi, Theta, xij)
USE Numeric_Kinds_Module
USE Numbers_Module
USE FileUnits_Module
IMPLICIT NONE
INTEGER k, IS, NS, JT, NT, KM
REAL(Kind=WP_Kind) xmin, xmax, Sin2Chi, Cos2Chi, CapS, CapSDots, Denominator, TanTheta
REAL(Kind=WP_Kind) RHO, S(0:NS), CapTheta(0:NT), Chi(0:NS,0:NT), Theta(0:NS,0:NT), xij(0:7)
KM=1
Xij(0)=  -Pi/6.d0
Xij(1)=   Zero
Xij(2)=   Pi/6.d0
Xij(3)= 2*Pi/6.d0
Xij(4)= 3*Pi/6.d0
Xij(5)= 4*Pi/6.d0
Xij(6)= 5*Pi/6.d0
Xij(7)= 6*Pi/6.d0

xmin=-Pi/12
xmax= Pi/12

OPEN(UNIT=GridNodes_Unit,FILE='Nodes.txt')
OPEN(UNIT=GridDebug_Unit,FILE='Debug.txt')
Chi=-1024.d0

DO IS=0,NS
   WRITE(GridDebug_Unit,*)"IS=",IS
   DO k=1,KM
      WRITE(GridDebug_Unit,*)"---- k=",k
      DO JT=0,NT
         Caps=SQRT(Rho**2-s(IS)**2)
         CapSDots=CapS*s(IS)*Cos(CapTheta(JT))
         Denominator=SQRT((CapS**2-s(IS)**2)**2+(Two*CapSDots)**2)
         Sin2Chi=Two*CapSDots/Denominator
         Cos2Chi=(CapS**2-s(IS)**2)/Denominator
         Chi(IS,JT)=Half*ATAN2(Sin2Chi,Cos2Chi)
         TanTheta=Denominator/(Two*CapS*s(IS)*Sin(CapTheta(JT)))
         Theta(IS,JT)=ATAN(TanTheta)
         !IF(Chi(IS,JT)-xij(k).ge.-(Xij(k)+Xij(k+1))/2.and.Chi(IS,JT)-xij(k).le.-(Xij(k)+Xij(k-1))/2)THEN
            WRITE(GridNodes_Unit,*)Theta(IS,JT),Chi(IS,JT)-xij(k)
            WRITE(GridDebug_Unit,'(4es14.7)')S(IS),CapTheta(JT),Theta(IS,JT),Chi(IS,JT)-xij(k)
         !ENDIF
      ENDDO
   ENDDO
   WRITE(GridNodes_Unit,'(A1)')'&'
ENDDO

!IF(Pi.ne.0)RETURN

DO k=1,KM
   DO JT=0,NT
      DO IS=0,NS
         WRITE(GridNodes_Unit,*)Theta(IS,JT),Chi(IS,JT)-xij(k)
      ENDDO
      WRITE(GridNodes_Unit,'(A1)')'&'
   ENDDO
ENDDO



RETURN
END SUBROUTINE Grid
