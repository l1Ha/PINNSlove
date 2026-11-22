SUBROUTINE sfplot_FEM_main
USE Numeric_Kinds_Module
USE FileUnits_Module
USE FileUnits_Common_Module
USE InputFile_Module
USE Narran_Module
USE Masses_Module
USE Numbers_Module
USE ChiAng_Module
USE Parms_Module
USE Numbers_Module
USE PES_Module
USE Convrsns_Module
IMPLICIT NONE
!*--SFPLOT_FEM_MAIN12
!
!   $RCSfile: MAIN_g3d.f,v $ $Revision: 1.13 $
!   $Date: 89/10/09 14:39:40 $
!   $State: Stable $
!
CHARACTER(LEN=25) :: title
 
LOGICAL :: there, Plot_Sfunc
INTEGER NTheta_Plot, ithmode, NChi_Plot, itheta, ichi, nang, nmode, ithrho, ienergy
REAL(KIND=WP_Kind) :: thetaval(mxang), chivals(mxang), Pot(mxang), energy(10)
REAL(KIND=WP_Kind) :: dtheta, dchi
REAL(KIND=WP_Kind) :: theta8(mxang), chi8(mxang)
REAL(KIND=WP_Kind) :: rho, theta, chi, phi, valmax, a, b, c, phi_s, rhofem
REAL(KIND=WP_Kind) :: phinow(mxang*mxbasis)
DATA energy/0.1D0, 0.2D0, 0.3D0, 0.4D0, 0.5D0, 0.6D0, 0.7D0, 0.8D0, 0.9D0,      &
   & 1.0D0/
 
NAMELIST/sfplot_FEM/ rho, NTheta_Plot, a, b, c, phi_s, ithmode, Plot_Sfunc

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')

INQUIRE(FILE='plot.vecs',EXIST=there)
IF(.NOT.there)THEN
   OPEN(Unit=PlotVecs_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'SfunctionPlots/PlotVecs.txt')
   CLOSE(UNIT=PlotVecs_Unit)
ENDIF
!
! OPEN output files
!

WRITE(Out_Unit,*)'Called    SFPlot_FEM_Main'
CLOSE(Out_Unit)
OPEN(UNIT=Out_Unit,FILE='sfplot_FEM.out',STATUS='unknown')
OPEN(UNIT=Plot_Dat_Unit,FILE='plot.dat',STATUS='unknown')
OPEN(UNIT=MathPoly_Unit,FILE='Math.Poly',STATUS='unknown')
 
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
READ(In_Unit,99001)title
 
WRITE(Out_Unit,99003)title
READ(In_Unit,99001)title
 
WRITE(Out_Unit,99003)title
READ(In_Unit,99002)title
 
WRITE(Out_Unit,99003)title
READ(In_Unit,99002)PES_Name
WRITE(Out_Unit,99003)PES_Name
CLOSE(UNIT=In_Unit)

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
READ(In_Unit,mass)
 
CLOSE(UNIT=In_Unit)
CALL constant
 
CALL massfac(.TRUE.)
 
READ(In_Unit,sfplot_FEM)
 
NChi_Plot = 8*NTheta_Plot - 7
IF(NTheta_Plot*NChi_Plot>mxang)THEN
   WRITE(Out_Unit,*)'NTheta_Plot*NChi_Plot must be less than mxang'
   STOP 'main'
ENDIF

IF(NTheta_Plot<=0 .OR. NTheta_Plot>mxang)THEN
   WRITE(Out_Unit,*)'Error NTheta_Plot =', NTheta_Plot
   STOP 'error'
ENDIF


IF(NChi_Plot<=0 .OR. NChi_Plot>mxang)THEN
   WRITE(Out_Unit,*)'Error NChi_Plot =', NChi_Plot
   STOP 'error'
ENDIF

dtheta = halfpi/(NTheta_Plot-1)
DO itheta = 1, NTheta_Plot
   thetaval(itheta) = dtheta*(itheta-1)
ENDDO

dchi = twopi/(NChi_Plot-1)
DO ichi = 1, NChi_Plot
   chivals(ichi) = dchi*(ichi-1)
ENDDO

DO ichi = 1, NChi_Plot
   DO itheta = 1, NTheta_Plot
       theta8(itheta+(ichi-1)*NTheta_Plot) = dtheta*(itheta-1)
       chi = dchi*(ichi-1)
       chi8(itheta+(ichi-1)*NTheta_Plot) = chi
   ENDDO
ENDDO

rhofem = rho
nang = NChi_Plot*NTheta_Plot
CALL basfem(theta8,chi8,nang,phinow,nmode,rhofem,ithrho)
IF(abs(rho-rhofem)>1.E-7)THEN
   WRITE(Out_Unit,*)'Error rhofem/=rho'
   WRITE(Out_Unit,*)'rhofem,rho', rhofem, rho
   WRITE(Msg_Unit,*)'Error rhofem/=rho'
   WRITE(Msg_Unit,*)'rhofem,rho', rhofem, rho
   STOP 'main'
ENDIF

WRITE(Out_Unit,*)'ithmode=', ithmode
IF(ithmode>nmode .OR. ithmode<=0)THEN
   WRITE(Out_Unit,*)'ithmode must be <= nmode'
   WRITE(Out_Unit,*)'ithmode,nmode=', ithmode, nmode
   STOP 'main'
ENDIF

valmax = 0.D0
DO ichi = 1, NChi_Plot
   DO itheta = 1, NTheta_Plot
       theta = thetaval(itheta)
       chi = chivals(ichi)
       phi = phinow(itheta+(ichi-1)*NTheta_Plot+(ithmode-1)*nang)
       phi = phi*autoev
       IF(phi<-phi_s)THEN
           phi = -phi_s*(1.0+c*log(abs(phi/phi_s)))
       ELSEIF(phi>phi_s)THEN
           phi = phi_s*(1.0+c*log(phi/phi_s))
       ENDIF
       phi = (1.0+b*phi)*a
       IF(phi<1.0E-30)phi = 1.0E-30
       Pot(itheta+(ichi-1)*NTheta_Plot) = phi
       IF(phi>valmax)valmax = phi
   ENDDO
ENDDO

WRITE(Out_Unit,99004)rho
WRITE(Out_Unit,99004)1.1*valmax*cos(chif1(1)), 1.1*valmax*sin(chif1(1)), zero
WRITE(Out_Unit,99004)1.1*valmax*cos(chif1(2)), 1.1*valmax*sin(chif1(2)), zero
WRITE(Out_Unit,99004)1.1*valmax*cos(chif1(3)), 1.1*valmax*sin(chif1(3)), zero
WRITE(Spin_Unit,99004)rho
WRITE(Spin_Unit,99004)1.1*valmax*cos(chif1(1)), 1.1*valmax*sin(chif1(1)), zero
WRITE(Spin_Unit,99004)1.1*valmax*cos(chif1(2)), 1.1*valmax*sin(chif1(2)), zero
WRITE(Spin_Unit,99004)1.1*valmax*cos(chif1(3)), 1.1*valmax*sin(chif1(3)), zero

DO ienergy = 1, 10
   phi = energy(ienergy)
   IF(phi>phi_s)phi = phi_s*(1.0+c*log(phi/phi_s))
   phi = (1.0+b*phi)*a
   WRITE(Spin_Unit,99004)phi, energy(ienergy)
ENDDO
 
CALL VecDef(thetaval,NTheta_Plot,chivals,NChi_Plot,Pot,valmax)

CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/BasisForScatt.txt',Form='FORMATTED',ACCESS='APPEND')
WRITE(Out_Unit,*)'Completed SFPlot_FEM_Main'
WRITE(Out_Unit,*)

STOP

99001 FORMAT(a25)
99002 FORMAT(a25)
99003 FORMAT(1x,a70)
99004 FORMAT(1x,3E16.7)
ENDSUBROUTINE SFPLOT_FEM_MAIN
