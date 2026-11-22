SUBROUTINE sfplot_ABM_main(Nbasis, Nbasiss, Neigmin, rho, nang, ithrho, MegaVal)
USE Numeric_Kinds_Module
USE Numbers_Module
USE FileUnits_Module
USE Narran_Module
USE Masses_Module
USE chiang_Module
USE Parms_Module
USE convrsns_Module
USE converge_Module
USE totj_Module
USE quantb_Module
USE gaussb_Module
USE InputFile_Module
IMPLICIT NONE
CHARACTER(LEN=5) RhoName, MegaName
CHARACTER(LEN=5), PARAMETER:: blank='     '
INTEGER Plot_sf
INTEGER itheta, ichi, ienergy, NTheta_Plot, nang, maxn, maxl, NChi_Plot, scheme, nthetaChi_Plot, nsfn, iang, i, j, iarran, vib
INTEGER :: nbasis(0:mxmega), nbasiss(0:mxmega), neigmin(0:mxmega), ithmode, ithrho, MegaVal, NABM_Sfunc, NSYM_States, NFull_States
LOGICAL :: primitive, Plot_Sfunc
LOGICAL :: there
REAL(KIND=SP_Kind) SP_PSI     ! Single Precision for Plotting
REAL(KIND=WP_Kind) u, v, dchi, dtheta, energy(10)
REAL(KIND=WP_Kind) rho, theta, chi, phi, valmax, a, b, c, phi_s, rhoabm, SignATMax
REAL(KIND=WP_Kind), ALLOCATABLE:: thetaval(:), chivals(:), Phi_Plot(:), Phi_ABM(:,:)
REAL(KIND=WP_Kind), ALLOCATABLE:: vect(:,:), thetaonly(:), chionly(:)
REAL(KIND=WP_Kind), ALLOCATABLE:: Theta_Plot(:), Chi_Plot(:), weight(:)
REAL(KIND=WP_Kind), ALLOCATABLE:: phinow(:,:), sfun(:,:)
REAL(KIND=WP_Kind), ALLOCATABLE:: pn(:,:,:), dpndx(:,:,:)
REAL(KIND=WP_Kind), ALLOCATABLE:: oscil(:,:,:)
REAL(KIND=WP_Kind), ALLOCATABLE:: doscdz(:,:,:), parityfn(:,:)
REAL(KIND=WP_Kind), ALLOCATABLE:: bfaci(:,:), eigen(:,:)
DATA energy/0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00/
 
NAMELIST/sfplot_ABM/ NTheta_Plot, a, b, c, phi_s, ithmode, primitive, Plot_Sfunc

rhoabm=rho
INQUIRE(file=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,exist=there)
IF(there)THEN
   OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
   READ(In_Unit,NML=sfplot_ABM, IOSTAT=plot_sf)
   IF(plot_sf < 0) Plot_SFunc=.False.
   NChi_Plot = 8*NTheta_Plot - 7
   nthetaChi_Plot=NTheta_Plot*NChi_Plot
   nsfn=NBasis(0)
   rho=rhoabm
   CLOSE(In_Unit)
   IF(.NOT.Plot_Sfunc)RETURN
ELSE
   WRITE(Msg_Unit,*)InputDIR(1:LEN(TRIM(InputDIR)))//InputFile, " does not exist."
   Plot_SFunc=.False.
   RETURN
ENDIF

NABM_Sfunc=NEigMin(MegaVal)
NSYM_States=NBasiss(MegaVal)
NFull_States=NBasis(MegaVal)
IF(.NOT.ALLOCATED(ThetaOnly))ALLOCATE(ThetaOnly(NTheta_Plot))
IF(.NOT.ALLOCATED(ChiOnly))ALLOCATE(ChiOnly(NChi_Plot))
IF(.NOT.ALLOCATED(ThetaVal))ALLOCATE(thetaval(nang))
IF(.NOT.ALLOCATED(ChiVals))ALLOCATE(chivals(nang))
IF(.NOT.ALLOCATED(Weight))ALLOCATE(weight(nang))
IF(.NOT.ALLOCATED(Phi_Plot))ALLOCATE(Phi_Plot(nthetaChi_Plot))
IF(.NOT.ALLOCATED(Phi_ABM))ALLOCATE(Phi_ABM(nang,NFull_States))
IF(.NOT.ALLOCATED(PhiNow))ALLOCATE(phinow(nthetaChi_Plot,NFull_States))
IF(.NOT.ALLOCATED(SFUN))ALLOCATE(sfun(nthetaChi_Plot,NABM_Sfunc))
IF(.NOT.ALLOCATED(ParityFn))ALLOCATE(parityfn(nthetaChi_Plot,narran))
IF(.NOT.ALLOCATED(Bfaci))ALLOCATE(bfaci(nthetaChi_Plot,narran))
IF(.NOT.ALLOCATED(pn))ALLOCATE(pn(0:mxl,nthetaChi_Plot,narran))
IF(.NOT.ALLOCATED(dpndx))ALLOCATE(dpndx(0:mxl,nthetaChi_Plot,narran))
IF(.NOT.ALLOCATED(oscil))ALLOCATE(oscil(0:mxl,nthetaChi_Plot,narran))
IF(.NOT.ALLOCATED(doscdz))ALLOCATE(doscdz(0:mxl,nthetaChi_Plot,narran))
IF(.NOT.ALLOCATED(eigen))ALLOCATE(eigen(nsfn,0:mxmega))
IF(.NOT.ALLOCATED(vect))ALLOCATE(vect(NSym_States,NABM_Sfunc))
IF(.NOT.ALLOCATED(Theta_Plot))ALLOCATE(Theta_Plot(nthetaChi_Plot))
IF(.NOT.ALLOCATED(Chi_Plot))ALLOCATE(Chi_Plot(nthetaChi_Plot))

rhoabm = rho
INQUIRE(FILE='plot.vecs',EXIST=there)
IF(.NOT.there)THEN
   OPEN(Unit=PlotVecs_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'SfunctionPlots/PlotVecs.txt')
   CLOSE(UNIT=PlotVecs_Unit)
ENDIF
!
! OPEN output files
!
WRITE(Out_Unit,*)'Called    SFPlot_ABM_Main'
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/SFPlot_ABM_Main.txt')
OPEN(Unit=Plot_txt_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'SfunctionPlots/plot.txt')
OPEN(Unit=MathPoly_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'SfunctionPlots/Math.Poly')
OPEN(Unit=Spin_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'SfunctionPlots/Spin.txt')

rho=rhoabm
IF(abs(rho-rhoabm)>1.E-7)THEN
   WRITE(Out_Unit,*)'Error rhoabm/=rho'
   WRITE(Out_Unit,*)'rhoabm,rho', rhoabm, rho
   WRITE(Msg_Unit,*)'Error rhoabm/=rho'
   WRITE(Msg_Unit,*)'rhoabm,rho', rhoabm, rho
   STOP "SFPlot_ABM"
ENDIF

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
READ(In_Unit,NML=quantum)
CLOSE(In_Unit)
DO iarran=1,narran
   nlegndre(iarran)=jmax(0,iarran)
   DO Vib=MinVib(iarran),MaxVib(iarran)
      jmax(Vib,iarran)=jmax(0,iarran)
   ENDDO
ENDDO
WRITE(Out_Unit,NML=quantum)

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
READ(In_Unit,NML=gauss)
noscil=maxvib+1
CLOSE(In_Unit)
WRITE(Out_Unit,NML=gauss)

IF(NTheta_Plot<=0 .OR. NTheta_Plot>mxang)THEN
   WRITE(Out_Unit,*)'Error NTheta_Plot =', NTheta_Plot
   STOP 'error'
ENDIF

IF(NChi_Plot<=0 .OR. NChi_Plot>mxang)THEN
   WRITE(Out_Unit,*)'Error NChi_Plot =', NChi_Plot
   STOP 'error'
ENDIF

dtheta = halfpi/(NTheta_Plot-1)
dchi= 4.d0*halfpi/(NChi_Plot-1)
DO ichi = 1, NChi_Plot
   chi = dchi*(ichi-1)
   ChiOnly(ichi)=chi
   DO itheta = 1, NTheta_Plot
       Theta=dtheta*(itheta-1)
       ThetaOnly(Itheta)=Theta
       Theta_Plot(itheta+(ichi-1)*NTheta_Plot) = Theta
       Chi_Plot(itheta+(ichi-1)*NTheta_Plot) = chi
   ENDDO
ENDDO

CALL readbas(rho, NFull_States, MegaVal, nang, Phi_ABM, ThetaVal, ChiVals, weight, vect, NABM_Sfunc, eigen(1,MegaVal), NSYM_States, ithrho)

IF(NFull_States>mxbasis)THEN
   WRITE(Out_Unit,*)'mxbasis<NFull_States:', mxbasis, NFull_States
   WRITE(Msg_Unit,*)'mxbasis<NFull_States:', mxbasis, NFull_States
   STOP "SFPlot_ABM"
ENDIF

WRITE(Out_Unit,*)'nang,mxang,NFull_States,mxbasis'
WRITE(Out_Unit,*)nang, mxang, NFull_States, mxbasis
CALL bascal2(nthetaChi_Plot, Theta_Plot, Chi_Plot, NFull_States, phinow, pn, dpndx, oscil, doscdz, parityfn, bfaci)
!  --------------------------------------------------------------
!  DO similarity transform to get overlaps of surface functions.
!  snew=s*vect contains overlaps of surface fcns.
!  -----------------------------------------------------------------
IF(primitive)THEN
   WRITE(Out_Unit,*)'Primitive basis functions being plotted'
ELSE
   WRITE(Out_Unit,*)'Constructing Actual Basis Functions'
   !CALL sgemm_junk(nthetaChi_Plot,NSym_States,NABM_Sfunc,phinow,nthetaChi_Plot,vect,NSym_States,sfun,NSym_States,0,1)
   CALL DGEMM('N','N',nthetaChi_Plot,NABM_Sfunc,NSym_States,1.d0,phinow,nthetaChi_Plot,vect,NSym_States,0.d0,sfun,nthetaChi_Plot)
   CALL vcopy(nthetaChi_Plot*NABM_Sfunc,phinow,sfun)
ENDIF

WRITE(Out_Unit,*)'ithmode=', ithmode
IF(ithmode<=0)THEN
   WRITE(Out_Unit,*)'ithmode must be greater than 0'
   WRITE(Msg_Unit,*)'ithmode must be greater than 0'
   STOP "SFPlot_ABM"
ELSEIF(primitive .AND. ithmode>NFull_States)THEN
   WRITE(Out_Unit,*)'ithmode must be less than NFull_States'
   WRITE(Msg_Unit,*)'ithmode must be less than NFull_States'
   STOP "SFPlot_ABM"
ELSEIF(.NOT.primitive .AND. ithmode>NABM_Sfunc)THEN
   WRITE(Out_Unit,*)'ithmode must be less than NABM_Sfunc'
   WRITE(msg_Unit,*)'ithmode must be less than NABM_Sfunc'
   STOP "SFPlot_ABM"
ENDIF

RhoName =blank
IF(ithrho<10)THEN
   WRITE(RhoName,'(i1)') ithrho
ELSEIF(ithrho<100)THEN
   WRITE(RhoName,'(i2)') ithrho
ELSEIF(ithrho<1000)THEN
   WRITE(RhoName,'(i3)') ithrho
ELSEIF(ithrho<10000)THEN
   WRITE(RhoName,'(i4)') ithrho
ELSEIF(ithrho<100000)THEN
   WRITE(RhoName,'(i5)') ithrho
ELSE
   WRITE(Msg_Unit,*)'ithrho to large=',ithrho
   STOP 'SFPlot_ABM'
ENDIF

MegaName =blank
IF(MegaVal<10)THEN
   WRITE(MegaName,'(i1)') MegaVal
ELSEIF(MegaVal<100)THEN
   WRITE(MegaName,'(i2)') MegaVal
ELSEIF(MegaVal<1000)THEN
   WRITE(MegaName,'(i3)') MegaVal
ELSEIF(MegaVal<10000)THEN
   WRITE(MegaName,'(i4)') MegaVal
ELSEIF(MegaVal<100000)THEN
   WRITE(MegaName,'(i5)') MegaVal
ELSE
   WRITE(Msg_Unit,*)'MegaVal to large=',MegaVal
   STOP 'SFPlot_ABM'
ENDIF

OPEN(Unit=Plt_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'SfunctionPlots/Sfunc_ABM_'//Rhoname(1:LEN(TRIM(RhoName)))//'_'//Meganame(1:LEN(TRIM(MegaName)))//'.txt')
valmax = 0.D0
SignAtMax = 0.D0
DO ichi = 1, NChi_Plot
   DO itheta = 1, NTheta_Plot
       phi = phinow(itheta+(ichi-1)*NTheta_Plot,ithmode)
       IF(ValMax<ABS(Phi))THEN
          ValMax=ABS(Phi)
          SignATMax=SIGN(1.d0,Phi)
       ENDIF
   ENDDO
ENDDO

DO ichi = 1, NChi_Plot
   DO itheta = 1, NTheta_Plot
       theta = thetaonly(itheta)
       chi = chionly(ichi)
       phi = phinow(itheta+(ichi-1)*NTheta_Plot,ithmode)
       CALL PMap(1,1,Zero,theta,chi,u,v)
       SP_PSI=phi*SignAtMax/ValMax
       IF(ABS(SP_PSI)<EPSILON(1.0_SP_Kind))SP_PSI=0.0
       WRITE(Plt_Unit,'(3ES12.4)')u,v,SP_PSI
       Phi_Plot(itheta+(ichi-1)*NTheta_Plot) = phi*SignAtMax/ValMax
   ENDDO
ENDDO
CLOSE(Plt_Unit)

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
 
CALL VecDef(thetaonly,NTheta_Plot,chionly,NChi_Plot,Phi_Plot,valmax)

DEALLOCATE(thetaval)
DEALLOCATE(chivals)
DEALLOCATE(Phi_Plot)
DEALLOCATE(Theta_Plot)
DEALLOCATE(Chi_Plot)
DEALLOCATE(weight)
DEALLOCATE(phinow)
DEALLOCATE(sfun)
DEALLOCATE(parityfn)
DEALLOCATE(bfaci)
DEALLOCATE(pn)
DEALLOCATE(dpndx)
DEALLOCATE(oscil)
DEALLOCATE(doscdz)
DEALLOCATE(eigen)
DEALLOCATE(vect)

CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/BasisForScatt.txt',Form='FORMATTED',ACCESS='APPEND')
WRITE(Out_Unit,*)'Completed SFPlot_ABM_Main'
WRITE(Out_Unit,*)

RETURN
99001 FORMAT(a25)
99002 FORMAT(a25)
99003 FORMAT(1x,a70)
 
99004 FORMAT(1x,3E16.7)
ENDSUBROUTINE SFPLOT_ABM_MAIN
