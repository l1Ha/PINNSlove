SUBROUTINE QuantumStates(Jtot,jmax)
IMPLICIT NONE
LOGICAL :: TestAllPosibilities=.False.
INTEGER :: Jtot, jmax, test
LOGICAL :: Parity      = .False.
LOGICAL :: PosPar      = .True.
LOGICAL :: Homonuclear = .False.
LOGICAL :: Evenj       = .True.

CALL NumberOfStates(Jtot, jmax)
IF(TestAllPosibilities)THEN
  DO test=1,9
     IF(Test==1)THEN
        Parity=.False.
        Homonuclear=.False.
     ELSEIF(Test==2)THEN
        Parity=.True.
        PosPar=.True.
        Homonuclear=.False.
     ELSEIF(Test==3)THEN
        Parity=.True.
        PosPar=.False.
        Homonuclear=.False.
     ELSEIF(Test==4)THEN
        Parity=.False.
        Homonuclear=.True.
        Evenj=.True.
     ELSEIF(Test==5)THEN
        Parity=.False.
        Homonuclear=.True.
        Evenj=.False.
     ELSEIF(Test==6)THEN
        Parity=.True.
        PosPar=.True.
        Homonuclear=.True.
        Evenj=.True.
     ELSEIF(Test==7)THEN
        Parity=.True.
        PosPar=.True.
        Homonuclear=.True.
        Evenj=.False.
     ELSEIF(Test==8)THEN
        Parity=.True.
        PosPar=.False.
        Homonuclear=.True.
        Evenj=.True.
     ELSE
        Parity=.True.
        PosPar=.False.
        Homonuclear=.True.
        Evenj=.False.
     ENDIF
     CAll j_and_l_States(Jtot, jmax, Parity, PosPar, HomoNuclear, Evenj)
  ENDDO
ELSE
   CAll j_and_l_States(Jtot, jmax, Parity, PosPar, HomoNuclear, Evenj)
ENDIF

RETURN
ENDSUBROUTINE QuantumStates
   
SUBROUTINE j_and_l_States(Jtot, jmax, Parity, PosPar, HomoNuclear, Evenj)
USE FileUnits_Common_Module
IMPLICIT NONE
INTEGER :: Jtot, jmax, QS_Unit
INTEGER Ntot, NPosPar, NNegPar, NEvenj, NOddj, NPosParEvenj, NPosParOddj, NNegParEvenj, NNegParOddj
INTEGER j,l
LOGICAL Parity, PosPar, Homonuclear, Evenj
QS_Unit=Out_Unit


WRITE(QS_Unit,*)"--------------------------------------------------"
WRITE(QS_Unit,*)
WRITE(QS_Unit,*)"Parity=", Parity," PosPar=", PosPar," Homonuclear=", Homonuclear," Evenj=",Evenj
WRITE(QS_Unit,*)
IF(.not.Parity)WRITE(QS_Unit,*)"Both even and odd parity states"
IF(Parity.and.PosPar)WRITE(QS_Unit,*)"Positive parity states"
IF(Parity.and..not.PosPar)WRITE(QS_Unit,*)"Odd parity states"

IF(.not.Homonuclear)WRITE(QS_Unit,*)"Both even and odd rotational states"
IF(Homonuclear.and.Evenj)WRITE(QS_Unit,*)"Even rotational states"
IF(Homonuclear.and..not.Evenj)WRITE(QS_Unit,*)"Odd rotational states"

WRITE(QS_Unit,*)
WRITE(QS_Unit,'(2(A,I5))')" Jtot=", Jtot, " jmax=", jmax
Ntot=0           ! Total number of states
NPosPar=0        ! Number of positive parity states
NNegPar=0        ! Number of negative parity states
NEvenj=0         ! Number of even rotational states
NOddj=0          ! Number of odd rotational states
NPosParEvenj=0   ! Number of even rotational states with positive parity
NPosParOddj=0    ! Number of odd rotational states with positive parity
NNegParEvenj=0   ! Number of even rotational states with negaive parity
NNegParOddj=0    ! Number of odd rotational states with odd parity
DO j=0,jmax
   DO l=ABS(Jtot-j),Jtot+j
      Ntot=Ntot+1
      IF(MOD(j+l,2)==0)NPosPar=NPosPar+1 ! Positiuve parity
      IF(MOD(j+l,2)==1)NNegPar=NNegPar+1 ! Positiuve parity
      IF(MOD(j,2)==0)NEvenj=Nevenj+1     ! Even rotational states
      IF(MOD(j,2)==1)Noddj=Noddj+1       ! Even rotational states
      IF(MOD(j+l,2)==0.and.MOD(j,2)==0)NPosParEvenj=NPosParEvenj+1  ! Positiuve parity and even j
      IF(MOD(j+l,2)==0.and.MOD(j,2)==1)NPosParOddj=NPosParOddj+1    ! Positiuve parity and odd j
      IF(MOD(j+l,2)==1.and.MOD(j,2)==0)NNegParEvenj=NNegParEvenj+1  ! Negative parity and even j
      IF(MOD(j+l,2)==1.and.MOD(j,2)==1)NNegParOddj=NNegParOddj+1    ! Negative parity and odd j
      IF(.not.Parity.and..not.Homonuclear)THEN
         WRITE(QS_Unit,'(4(A,I5))')"Ntot=",Ntot," j=",j," l=",l
      ELSEIF(Parity.and..not.Homonuclear)THEN
         IF(PosPar.and.MOD(j+l,2)==0)WRITE(QS_Unit,'(4(A,I5))')"NPosPar=",NPosPar," j=",j," l=",l
         IF(.not.PosPar.and.MOD(j+l,2)==1)WRITE(QS_Unit,'(4(A,I5))')"NNegPar=",NNegPar," j=",j," l=",l 
      ELSEIF(.not.Parity.and.Homonuclear)THEN
         IF(Evenj.and.MOD(j,2)==0)WRITE(QS_Unit,'(4(A,I5))')"NEvenj=",NEvenj," j=",j," l=",l
         IF(.not.Evenj.and.MOD(j,2)==1)WRITE(QS_Unit,'(4(A,I5))')"NOddj=",NOddj," j=",j," l=",l
      ELSEIF(Parity.and.Homonuclear)THEN
         IF(PosPar.and.MOD(j+l,2)==0)THEN
            IF(Evenj.and.MOD(j,2)==0)WRITE(QS_Unit,'(4(A,I5))')"NPosParEvenj=",NPosParEvenj," j=",j," l=",l
            IF(.not.Evenj.and.MOD(j,2)==1)WRITE(QS_Unit,'(4(A,I5))')"NPosParOddj=",NPosParOddj," j=",j," l=",l
         ELSEIF(.not.PosPar.and.MOD(j+l,2)==1)THEN
            IF(Evenj.and.MOD(j,2)==0)WRITE(QS_Unit,'(4(A,I5))')"NNegParEvenj=",NNegParEvenj," j=",j," l=",l
            IF(.not.Evenj.and.MOD(j,2)==1)WRITE(QS_Unit,'(4(A,I5))')"NNegParOddj=",NNegParOddj," j=",j," l=",l
         ENDIF
      ENDIF
   ENDDO
ENDDO

IF(Ntot/=NPosPar+NNegPar)Stop "Ntot/=NPosPar+NNegPar"
IF(Ntot/=NEvenj+NOddj)Stop "Ntot/=NEvenj+NOddj"
IF(Ntot/=NPosParEvenj+NPosParOddj+NNegParEvenj+NNegParOddj)Stop "Ntot/=NPosParEvenj+NPosParOddj+NNegParEvenj+NNegParOddj"

WRITE(QS_Unit,*)"--------------------------------------------------"

RETURN
ENDSUBROUTINE j_and_l_States
   
SUBROUTINE NumberOfStates(Jtot, jmax)
USE FileUnits_Common_Module
IMPLICIT NONE
INTEGER :: Jtot, jmax, QS_Unit
INTEGER Ntot, NPosPar, NNegPar, NEvenj, NOddj, NPosParEvenj, NPosParOddj, NNegParEvenj, NNegParOddj
INTEGER j,l
QS_Unit=Out_Unit

Ntot=0           ! Total number of states
NPosPar=0        ! Number of positive parity states
NNegPar=0        ! Number of negative parity states
NEvenj=0         ! Number of even rotational states
NOddj=0          ! Number of odd rotational states
NPosParEvenj=0   ! Number of even rotational states with positive parity
NPosParOddj=0    ! Number of odd rotational states with positive parity
NNegParEvenj=0   ! Number of even rotational states with negaive parity
NNegParOddj=0    ! Number of odd rotational states with odd parity
DO j=0,jmax
   DO l=ABS(Jtot-j),Jtot+j
      Ntot=Ntot+1
      IF(MOD(j+l,2)==0)NPosPar=NPosPar+1 ! Positiuve parity
      IF(MOD(j+l,2)==1)NNegPar=NNegPar+1 ! Positiuve parity
      IF(MOD(j,2)==0)NEvenj=Nevenj+1     ! Even rotational states
      IF(MOD(j,2)==1)Noddj=Noddj+1       ! Even rotational states
      IF(MOD(j+l,2)==0.and.MOD(j,2)==0)NPosParEvenj=NPosParEvenj+1  ! Positiuve parity and even j
      IF(MOD(j+l,2)==0.and.MOD(j,2)==1)NPosParOddj=NPosParOddj+1    ! Positiuve parity and odd j
      IF(MOD(j+l,2)==1.and.MOD(j,2)==0)NNegParEvenj=NNegParEvenj+1  ! Negative parity and even j
      IF(MOD(j+l,2)==1.and.MOD(j,2)==1)NNegParOddj=NNegParOddj+1    ! Negative parity and odd j
   ENDDO
ENDDO
WRITE(QS_Unit,*)
WRITE(QS_Unit,'(A,I5)')" Jtot=",jtot
WRITE(QS_Unit,'(A,I5)')" jmax=",jmax
WRITE(QS_Unit,'(A,I5)')" Ntot=",Ntot

WRITE(QS_Unit,*)
WRITE(QS_Unit,'(A,I5)')" NPosPar=",NPosPar
WRITE(QS_Unit,'(A,I5)')" NNegPar=",NNegPar
WRITE(QS_Unit,'(A,I5)')" NPosPar+NNegPar=",NPosPar+NNegPar

WRITE(QS_Unit,*)
WRITE(QS_Unit,'(A,I5)')" NEvenj=",NEvenj
WRITE(QS_Unit,'(A,I5)')" NOddj=",NOddj
WRITE(QS_Unit,'(A,I5)')" NEvenj+NOddj=",NEvenj+NOddj

WRITE(QS_Unit,*)
WRITE(QS_Unit,'(A,I5)')" NPosParEvenj",NPosParEvenj
WRITE(QS_Unit,'(A,I5)')" NPosParOddj=",NPosParOddj
WRITE(QS_Unit,'(A,I5)')" NNegParEvenj=",NNegParEvenj
WRITE(QS_Unit,'(A,I5)')" NNegParOddj=",NNegParOddj
WRITE(QS_Unit,'(A,I5)')" NPosParEvenj+NPosParOddj+NNegParEvenj+NNegParOddj=",NPosParEvenj+NPosParOddj+NNegParEvenj+NNegParOddj

IF(Ntot/=NPosPar+NNegPar)Stop "Ntot/=NPosPar+NNegPar"
IF(Ntot/=NEvenj+NOddj)Stop "Ntot/=NEvenj+NOddj"
IF(Ntot/=NPosParEvenj+NPosParOddj+NNegParEvenj+NNegParOddj)Stop "Ntot/=NPosParEvenj+NPosParOddj+NNegParEvenj+NNegParOddj"

   
RETURN
ENDSUBROUTINE NumberOfStates
   