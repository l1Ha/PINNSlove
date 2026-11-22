SUBROUTINE calc_pot2 (nrho, ntheta, nchi, rho_val, theta_val, chi_val, vpot, morph, morph_val)
USE Numeric_Kinds_Module
USE Numbers_Module
USE FileUnits_Module
USE PES_MODULE
USE Masses_Module
USE InputFile_Module
IMPLICIT NONE
LOGICAL :: const_rho, const_chi, const_theta, vmin_rho, vmin_theta, There
LOGICAL :: vmin_chi, morph
LOGICAL, Parameter :: Prt_Masses=.False.
INTEGER :: nrho, nchi, ntheta, irho, itheta, ichi, icase, nkeyfr, kfrcnt, I, IERR
REAL(Kind=WP_Kind) :: rho_min, rho_max, rho, chi,theta, upot, drho
REAL(Kind=WP_Kind) :: dchi, theta_min, theta_max, chi_min, chi_max
REAL(Kind=WP_Kind) :: vmin, radtodeg, morph_val
REAL(Kind=WP_Kind) :: rho_val (nrho), theta_val (ntheta), chi_val (nchi)
REAL(Kind=WP_Kind) :: vpot (nrho, ntheta, nchi)
REAL(Kind=WP_Kind), ALLOCATABLE:: temp(:,:)

NAMELIST/pot2/ rho_min, rho_max, const_rho, rho, vmin_rho, const_chi, chi, vmin_chi, const_theta, theta, vmin_theta
ALLOCATE(temp(ntheta,ntheta))
WRITE(Out_Unit, * )
WRITE(Out_Unit, * ) 'Begin(Calc_Pot2) Calculate Potential'

! Atomic Masses
INQUIRE (File=InputDIR(1:LEN(TRIM(InputDIR)))//TRIM(InputFile), exist = there)
IF(there)THEN
   REWIND In_Unit
   OPEN(Unit=In_Unit,File=TRIM(InputDIR)//TRIM(InputFile))
      READ(In_Unit,  NML=AtomicMasses, IOSTAT=IERR)                            ! Atomic Masses
      IF(IERR==-1)THEN
         IERR=0
         REWIND In_Unit
         READ(In_Unit, NML=Atoms, IOSTAT=IERR)
         IF(Prt_Masses)THEN
            WRITE(Out_Unit,NML=Atoms)
            WRITE(Msg_Unit,NML=Atoms)
            WRITE(Out_Unit,*)
         ENDIF
         DO I=1,3
            CALL AtomicWeights(AtomicSymbol(I), MassNumber(I), AtomicNumber(I), Mass(I), Abundance(I), AtomicWeight(I), Notes(I))
            IF(Prt_Masses)THEN
            WRITE(Out_Unit, '(A,I3)')'Atomic Number = ', AtomicNumber(I)
            WRITE(Out_Unit, '(A,A3)')'Atomic Symbol = ', TRIM(AtomicSymbol(I))
            WRITE(Out_Unit, '(A,I3)')'Mass Number = '  , MassNumber(I)
            WRITE(Out_Unit, '(A,1PES23.15)')'Relative Atomic Mass = '  , Mass(I)
            WRITE(Out_Unit, '(A,1PES23.15)')'Isotopic Composition = '  , Abundance(I)
            WRITE(Out_Unit, '(A,1PES23.15)')'Standard Atomic Weight = ', AtomicWeight(I)
            WRITE(Out_Unit, '(A,A)')'Notes = ', TRIM(Notes(I))
            WRITE(Out_Unit,*)
            ENDIF
         ENDDO
      ELSEIF (IERR/=0)THEN
         WRITE(Msg_Unit,*)' Missing a necessary namelist: Atoms or AtomicMasses'
      ENDIF
      WRITE(Out_Unit,NML=AtomicMasses)
   CLOSE(In_Unit)
ELSE
   WRITE(Out_Unit,  * ) 'Data file must exist: ', InputDIR(1:LEN(TRIM(InputDIR)))//InputFile
   WRITE(Msg_Unit,  * ) 'Data file must exist: ', InputDIR(1:LEN(TRIM(InputDIR)))//InputFile
   STOP "Calc_Pot2"
ENDIF


OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile)
READ(In_Unit,NML=pot2, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'Pot2')
      rho_min=1.d0
      rho_max=10.d0
      const_rho=.false.
      rho=9.5d0
      vmin_rho=.false.
      const_theta=.false.
      theta=Pi/two
      vmin_theta=.false.
      const_chi=.false.
      chi=0.0
      vmin_chi=.false.
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Pot2'
      WRITE(Msg_Unit,*)'ERROR with Namelist Pot2'
      STOP 'Calc_Pot2: Pot2'
   ENDIF
CLOSE(In_Unit)
WRITE(Out_Unit,NML=pot2)

CALL constant

CALL massfac (.true.)

IF(const_rho)THEN
   nrho = 1
   IF(morph)THEN
      rho = morph_val
   ENDIF
   rho_val (1) = rho
ELSE
   drho = (rho_max - rho_min) / (nrho - 1)
   DO irho = 1, nrho
      rho_val (irho) = (irho - 1) * drho + rho_min
   ENDDO
ENDIF
theta_min = 0.d0
theta_max = 2.d0 * atan (1.d0)
IF(const_theta)THEN
   ntheta = 1
   IF(morph)THEN
      theta = morph_val
   ENDIF
   theta_val (1) = theta
ELSE
!--------------------------------------
!set new theta values
!--------------------------------------
CALL DVR_kinetic(ntheta,temp,theta_val,One)
ENDIF
chi_min = 0.d0
chi_max = 8.d0 * atan (1.d0)
IF(const_chi)THEN
   nchi = 1
   IF(morph)THEN
      chi = morph_val
   ENDIF
   chi_val (1) = chi
ELSE
   dchi = (chi_max - chi_min) / (nchi - 1)
   DO ichi = 1, nchi
      chi_val (ichi) = (ichi - 1) * dchi
   ENDDO

ENDIF

kfrcnt = 0
DO ichi = 1, nchi
   chi = chi_val (ichi)
   DO itheta = 1, ntheta
      theta = theta_val (itheta)
      DO irho = 1, nrho
         rho = rho_val (irho)
         nkeyfr = nchi * ntheta * nrho
         kfrcnt = kfrcnt + 1
         CALL potaph2 (rho, theta, chi, upot, nkeyfr, kfrcnt)
         vpot (irho, itheta, ichi) = upot
      ENDDO
   ENDDO
ENDDO
IF(kfrcnt/=nkeyfr)THEN
   STOP 'kfrcnt_calc_pot.f'

ENDIF
IF(vmin_rho)THEN
   DO itheta = 1, ntheta
      DO ichi = 1, nchi
         vmin = Ten**(RANGE(One)-10)
         DO irho = 1, nrho
            IF(vpot (irho, itheta, ichi) <vmin)THEN
               vmin = vpot (irho, itheta, ichi)
            ENDIF
         ENDDO
         DO irho = 1, nrho
            vpot (irho, itheta, ichi) = vmin
         ENDDO
      ENDDO
   ENDDO
   nrho = 1
ENDIF
IF(vmin_theta)THEN
   DO irho = 1, nrho
      DO ichi = 1, nchi
         vmin = Ten**(RANGE(One)-10)
         DO itheta = 1, ntheta
            IF(vpot (irho, itheta, ichi) <vmin)THEN
               vmin = vpot (irho, itheta, ichi)
            ENDIF
         ENDDO
         DO itheta = 1, ntheta
            vpot (irho, itheta, ichi) = vmin
         ENDDO
      ENDDO
   ENDDO
   ntheta = 1

ENDIF
IF(vmin_chi)THEN
   DO irho = 1, nrho
      DO itheta = 1, ntheta
         vmin = Ten**(RANGE(One)-10)
         DO ichi = 1, nchi
            IF(vpot (irho, itheta, ichi) <vmin)THEN
               vmin = vpot (irho, itheta, ichi)
            ENDIF
         ENDDO
         DO ichi = 1, nchi
            vpot (irho, itheta, ichi) = vmin
         ENDDO
      ENDDO
   ENDDO
   nchi = 1
ENDIF
OPEN(Unit=PotDat_Bin_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/PotDat.bin')
OPEN(Unit=Pot_Graph_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'GraphicsOut/pot_graph.csv')
OPEN(Unit=Pot_Bin_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Pot.bin', form = 'unformatted')
WRITE(Out_Unit, '(1x,A5,3x,i5)') 'nrho=', nrho
WRITE(Out_Unit, '(1x,A7,1x,i5)') 'ntheta=', ntheta
WRITE(Out_Unit, '(1x,A5,3x,i5)') 'nchi=', nchi
WRITE(PotDat_Bin_Unit, * ) nrho, ntheta, nchi
WRITE(Pot_Bin_Unit) PES_Name
WRITE(Pot_Bin_Unit) mass(1), mass(2), mass(3), usys2
WRITE(Pot_Bin_Unit) rho_min, rho_max
WRITE(Pot_Bin_Unit) nrho, const_rho, rho, vmin_rho
WRITE(Pot_Bin_Unit) nchi, const_chi, chi, vmin_chi
WRITE(Pot_Bin_Unit) ntheta, const_theta, theta, vmin_theta
radtodeg = fortyfive / atan (1.d0)
IF(nrho>1.AND.ntheta>1.AND.nchi>1)THEN
   icase = 1
ELSEIF (nrho>1.AND.ntheta>1.AND.nchi==1)THEN
   icase = 2
ELSEIF (nrho>1.AND.ntheta==1.AND.nchi>1)THEN
   icase = 3
ELSEIF (nrho==1.AND.ntheta>1.AND.nchi>1)THEN
   icase = 4
ELSEIF (nrho>1.AND.ntheta==1.AND.nchi==1)THEN
   icase = 5
ELSEIF (nrho==1.AND.ntheta>1.AND.nchi==1)THEN
   icase = 6
ELSEIF (nrho==1.AND.ntheta==1.AND.nchi>1)THEN
   icase = 7
ENDIF
WRITE(Out_Unit, '(14x,A3,20x,A5,22x,A3,15x,A16)') 'Rho', 'Theta', 'Chi', 'V(Rho,Theta,Chi)'
DO irho = 1, nrho
   DO itheta = 1, ntheta
      DO ichi = 1, nchi
         WRITE(Out_Unit, 20) rho_val (irho), theta_val (itheta), chi_val (ichi), vpot (irho, itheta, ichi)
         WRITE(PotDat_Bin_Unit, 30) rho_val (irho), theta_val (itheta), chi_val (ichi),vpot (irho, itheta, ichi)
         IF(icase==1)THEN
            WRITE(Pot_Graph_Unit, 40) rho_val (irho), theta_val (itheta) * radtodeg, chi_val (ichi) * radtodeg, vpot (irho, itheta, ichi)
         ELSEIF (icase==2)THEN
            WRITE(Pot_Graph_Unit, 40) rho_val (irho), theta_val (itheta) * radtodeg, vpot (irho, itheta, ichi)
         ELSEIF (icase==3)THEN
            WRITE(Pot_Graph_Unit, 40) rho_val (irho), chi_val (ichi) * radtodeg, vpot (irho, itheta, ichi)
         ELSEIF (icase==4)THEN
            WRITE(Pot_Graph_Unit, 40) theta_val (itheta), chi_val (ichi) * radtodeg,  vpot (irho, itheta, ichi)
         ELSEIF (icase==5)THEN
            WRITE(Pot_Graph_Unit, 40) rho_val (irho), vpot (irho, itheta, ichi)
         ELSEIF (icase==6)THEN
            WRITE(Pot_Graph_Unit, 40) theta_val (itheta), vpot (irho, itheta, ichi)
         ELSEIF (icase==7)THEN
            WRITE(Pot_Graph_Unit, 40) chi_val (ichi), vpot (irho, itheta, ichi)
         ENDIF
         WRITE(Pot_Bin_Unit) rho_val (irho), theta_val (itheta), chi_val (ichi), vpot (irho, itheta, ichi)
      ENDDO
   ENDDO
ENDDO
CLOSE(unit=Pot_Bin_Unit)
CLOSE(unit=Pot_Graph_Unit)

CLOSE(unit=PotDat_Bin_Unit)
   20 FORMAT(1x,4d25.15)
   30 FORMAT(1x,4d25.15)
   40 FORMAT(1x,4(1pe15.7','))
WRITE(Out_Unit, * ) 'Begin(Calc_Pot2)'

WRITE(Out_Unit, * )

DEALLOCATE(Temp)
RETURN
ENDSUBROUTINE calc_pot2
