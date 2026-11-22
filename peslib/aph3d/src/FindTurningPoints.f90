 SUBROUTINE FindTurningPoints(NStates, Eig, Rvals, V, LeftTP, RightTP, Method, jrot, jvalue, NthElec)
USE Numeric_Kinds_Module
USE Numbers_Module
USE FileUnits_Common_Module
USE FileUnits_OneDim_Module
USE OneDim_Parms_Module
USE CriticalDist_Module
USE QState_Module
USE DiatomicPot_Module
USE Convrsns_Module
IMPLICIT NONE
CHARACTER(LEN=11) Method
CHARACTER(LEN=1) elecname
CHARACTER(LEN=5) jvalname
INTEGER(KIND=IW_Kind) NStates, IthState, I, jvalue, jrot(0:NStates), NthElec(0:NStates), nu_state, iloc(1), jthiter
REAL(KIND=WP_Kind)  DistMin(1), DistBest(1), DistMax(1), VMinimum(1), VMaximum(1), VBest(1), Vm(0:10,0:100)
REAL(KIND=WP_Kind)  r, hstep, vp(0:nmax), vpp(0:nmax), mu, we_au(0:100), wexe_au(1:100), weye_au(2:100), v0, v1, v2, tolfit
REAL(KIND=WP_Kind)  Eig(0:NStates), Rvals(0:Nmax), V(0:Nmax), LeftTP(0:NStates), RightTP(0:NStates)
REAL(KIND=WP_Kind), PARAMETER:: Tol=1.d-6
REAL(KIND=WP_Kind), PARAMETER:: r0=Zero

WRITE(*,*) "NStates=",NStates
jval=jvalue
! Find the left turning points
hstep=(ri-r0)/nmax
nu_state=0
DO IthState=0,NStates-1
   DO i=0,nmax
      r=i*hstep
      jval=jrot(IthState)
      elec=NthElec(IthState)
      CALL PotDiatomic(r, hstep, v0, v1, v2, jval, mu)
      v(i)=v0
      vp(i)=v1
      vpp(i)=v2
   ENDDO
   Vm(elec,jval)=MINVAL(V)
   iloc=MINLOC(ABS(V-Eig(IthState)))-1
   IF(NthElec(IthState).eq.elec)THEN
   DO I=1,Nmax
      IF((Eig(IthState)-V(I))*(Eig(IthState)-V(I-1))<=Zero)THEN
         DistMin=Rvals(I-1)
         DistMax=Rvals(I)
         DistBest=(DistMin+DistMax)/Two
         VMinimum=MIN(V(I),V(I-1))
         VMaximum=MAX(V(I),V(I-1))
         CALL AKIMA (3,Nmax+1, Rvals, V, 1, DistBest, VBest)
         EXIT
      ENDIF
   ENDDO
   DO jthiter=1,50
   !DO WHILE(ABS((DistMax(1)-DistMin(1)))>=Tol)
   !DO WHILE(ABS((DistMax(1)-DistMin(1))/Eig(IthState))>=Tol)
      IF(VBest(1)-Eig(IthState)>=Zero)THEN
         DistMin=DistBest
         DistMax=DistMax
         DistBest=(DistMin+DistMax)/Two
         VMinimum=VMinimum
         VMaximum=VBest
      ELSEIF(VBest(1)-Eig(IthState)<=Zero)THEN
         DistMin=DistMin
         DistMax=DistBest
         DistBest=(DistMin+DistMax)/Two
         VMinimum=VBest
         VMaximum=VMaximum
      ELSE
         EXIT
      ENDIF
         nu_state=nu_state+1
      CALL AKIMA (3,Nmax+1,Rvals,V,1,DistBest,VBest)
   ENDDO
   LeftTP(IthState)=DistBest(1)
   ENDIF
ENDDO

!RETURN  !TmpModGregParker

!Find the right turning points
DO IthState=0,NStates-1
   DO i=0,nmax
      r=i*hstep
      jval=jrot(IthState)
      elec=NthElec(IthState)
      CALL PotDiatomic(r, hstep, v0, v1, v2, jval, mu)
      v(i)=v0
      vp(i)=v1
      vpp(i)=v2
   ENDDO
   Vm(elec,jval)=MINVAL(V)
   iloc=MINLOC(ABS(V-Eig(IthState)))-1
   IF(iloc(1)<10)THEN
      write(*,*)iloc
   ENDIF
   IF(NthElec(IthState).eq.elec)THEN
   DO I=Nmax,iloc(1),-1
      IF((Eig(IthState)-V(I))*(Eig(IthState)-V(I-1))<=Zero)THEN
         DistMin=Rvals(I-1)
         DistMax=Rvals(I)
         DistBest=(DistMin+DistMax)/Two
         VMinimum=MIN(V(I),V(I-1))
         VMaximum=MAX(V(I),V(I-1))
         CALL AKIMA (3,Nmax+1,Rvals,V,1,DistBest,VBest)
         EXIT
      ENDIF
   ENDDO
   
   DO jthiter=1,50
   !DO WHILE(ABS((DistMax(1)-DistMin(1))/Eig(IthState))>=Tol)
      IF(VBest(1)-Eig(IthState)<=Zero)THEN
         DistMin=DistBest
         DistMax=DistMax
         DistBest=(DistMin+DistMax)/Two
         VMinimum=VMinimum
         VMaximum=VBest
      ELSEIF(VBest(1)-Eig(IthState)>=Zero)THEN
         DistMin=DistMin
         DistMax=DistBest
         DistBest=(DistMin+DistMax)/Two
         VMinimum=VBest
         VMaximum=VMaximum
      ELSE
         EXIT
      ENDIF
      CALL AKIMA (3,Nmax+1,Rvals,V,1,DistBest,VBest)
   ENDDO
   RightTP(IthState)=DistBest(1)
   ENDIF
ENDDO


DO jval=jmin,jmax
   IF(jval.lt.10)THEN
      WRITE(jvalname,'(i1)') jval
   ELSEIF(jval.lt.100)THEN
      WRITE(jvalname,'(i2)') jval
   ELSEIF(jval.lt.1000)THEN
      WRITE(jvalname,'(i3)') jval
   ELSEIF(jval.lt.10000)THEN
      WRITE(jvalname,'(i4)') jval
   ELSEIF(jval.lt.100000)THEN
      WRITE(jvalname,'(i5)') jval
   ELSE
      WRITE(Msg_Unit,*)'jval to large=',jval
      STOP 'FindTurningPoints'
   ENDIF
   WRITE(elecname,'(i1)') elec
   vmin(elec)=Vm(elec,jval)

   OPEN(Unit=Eig_Plot_Unit,File=TRIM(OutDIR)//"DiatomicPlots/"//TRIM(DiatomicPot)//"_EigPlot"//TRIM(Method)//"e"//elecname//"j"//jvalname(1:LEN(TRIM(jvalname)))//".csv")
   !Print Heading
   WRITE(Eig_Plot_Unit,'(12(A15,","))')"nvib","Turning Point","Eigenenergy", "we(au)", "wexe(au)", "weye(au)", "we(ev)", "wexe(eV)", "weye(eV)", "we(cm-1)", "wexe(cm-1)", "weye(cm-1)"
   ! Print rmin and corresponding vmin
   WRITE(Eig_Plot_Unit,'(A15,",",2(ES15.7,","))')"-.5", rmin(elec), vmin(elec)
   ! Print the Right Turning points and corresponding eigenenergy
   IthState=0
   nu_state=0
   we_au(0)=2*(Eig(0)-vmin(elec))
   DO IthState=1,NStates-1
      IF(jrot(IthState)==jval.and.elec==NthElec(IthState))THEN
         nu_state=nu_state+1
         we_au(nu_state)=Eig(nu_state)-Eig(nu_state-1)
      ENDIF
   ENDDO
   nu_state=0
   DO IthState=1,NStates-1
      IF(jrot(IthState)==jval.and.elec==NthElec(IthState))THEN
         nu_state=nu_state+1
         wexe_au(nu_state)=(we_au(nu_state)-we_au(nu_state-1))/two
      ENDIF
   ENDDO
   nu_state=1
   DO IthState=2,NStates-1
      IF(jrot(IthState)==jval.and.elec==NthElec(IthState))THEN
         nu_state=nu_state+1
         weye_au(nu_state)=(wexe_au(nu_state)-wexe_au(nu_state-1))/six
      ENDIF
   ENDDO
   nu_state=0
   WRITE(Eig_Plot_Unit,'(I15,",",3(ES15.7,","),TR15,",",TR15,",",(ES15.7,","),TR15,",",TR15,",",(ES15.7,","))')nu_state, LeftTP(0), Eig(0), 2*(Eig(0)-vmin(elec)), 2*(Eig(0)-vmin(elec))*autoev, 2*(Eig(0)-vmin(elec))/cmm1toau
   DO IthState=1,NStates-1
      IF(jrot(IthState)==jval.and.elec==NthElec(IthState))THEN
         nu_state=nu_state+1
         IF(nu_state>1)THEN
            WRITE(Eig_Plot_Unit,'(I15,",",11(ES15.7,","))')nu_state, LeftTP(IthState), Eig(IthState), we_au(nu_state), wexe_au(nu_state), weye_au(nu_state), we_au(nu_state)*autoev, wexe_au(nu_state)*autoev, weye_au(nu_state)*autoev, we_au(nu_state)/cmm1toau, wexe_au(nu_state)/cmm1toau, weye_au(nu_state)/cmm1toau
            ELSEIF(nu_State>0)THEN
            WRITE(Eig_Plot_Unit,'(I15,",",4(ES15.7,","),TR15,",",2(ES15.7,","),TR15,",",2(ES15.7,","))')nu_state, LeftTP(IthState), Eig(IthState), we_au(nu_state), wexe_au(nu_state), we_au(nu_state)*autoev, wexe_au(nu_state)*autoev, we_au(nu_state)/cmm1toau, wexe_au(nu_state)/cmm1toau
         ELSE
            WRITE(Eig_Plot_Unit,'(I15,",",3(ES15.7,","),TR15,",",TR15,",",(ES15.7,","),TR15,",",TR15,",",(ES15.7,","))')nu_state, LeftTP(IthState), Eig(IthState), we_au(nu_state), we_au(nu_state)*autoev, we_au(nu_state)/cmm1toau
         ENDIF
      ENDIF
   ENDDO
   
   ! Print rmin and corresponding vmin
   WRITE(Eig_Plot_Unit,*)"&"
   WRITE(Eig_Plot_Unit,'(A15,",",2(ES15.7,","))')"-.5", rmin(elec), vmin(elec)
   ! Print the Right Turning points and corresponding eigenenergy
   IthState=0
   nu_state=0
   WRITE(Eig_Plot_Unit,'(I15,",",3(ES15.7,","),TR15,",",TR15,",",(ES15.7,","),TR15,",",TR15,",",(ES15.7,","))')nu_state, RightTP(0), Eig(0), 2*(Eig(0)-vmin(elec)), 2*(Eig(0)-vmin(elec))*autoev, 2*(Eig(0)-vmin(elec))/cmm1toau
   DO IthState=1,NStates-1
      IF(jrot(IthState)==jval.and.elec==NthElec(IthState))THEN
         nu_state=nu_state+1
         IF(nu_state>1)THEN
            WRITE(Eig_Plot_Unit,'(I15,",",11(ES15.7,","))')nu_state, RightTP(IthState), Eig(IthState), we_au(nu_state), wexe_au(nu_state), weye_au(nu_state), we_au(nu_state)*autoev, wexe_au(nu_state)*autoev, weye_au(nu_state)*autoev, we_au(nu_state)/cmm1toau, wexe_au(nu_state)/cmm1toau, weye_au(nu_state)/cmm1toau
            ELSEIF(nu_State>0)THEN
            WRITE(Eig_Plot_Unit,'(I15,",",4(ES15.7,","),TR15,",",2(ES15.7,","),TR15,",",2(ES15.7,","))')nu_state, RightTP(IthState), Eig(IthState), we_au(nu_state), wexe_au(nu_state), we_au(nu_state)*autoev, wexe_au(nu_state)*autoev, we_au(nu_state)/cmm1toau, wexe_au(nu_state)/cmm1toau
         ELSE
            WRITE(Eig_Plot_Unit,'(I15,",",3(ES15.7,","),TR15,",",TR15,",",(ES15.7,","),TR15,",",TR15,",",(ES15.7,","))')nu_state, RightTP(IthState), Eig(IthState), we_au(nu_state), we_au(nu_state)*autoev, we_au(nu_state)/cmm1toau
         ENDIF
      ENDIF
   ENDDO
   FLUSH(Eig_Plot_Unit)
   CLOSE(Unit=Eig_Plot_Unit)
ENDDO

RETURN
ENDSUBROUTINE FindTurningPoints