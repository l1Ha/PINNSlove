SUBROUTINE Kmatrx_Main(Nopen)
!
! This program calculates the K_Matrix and/or derivative of
! the K_Matrix with respect to the total energy from the R_matrix.
!
!
! modified by m. braunstein, april, 1991, to include
! photodissociation option will take R matrix, gamma vector
! and find the real and imaginary photodissociation
! amplitudes
!
! gamma = inhomogeneity
! phor = real part of photodissociation amplitude
! phoi = imaginary part of photodissociation amplitude
!
!------------------------------------------------------------------
USE FileUnits_MODULE
USE FileNames_Module
USE QAsy_Numbers_Module
USE Masses_Module
USE parms_MODULE
!USE nstate_MODULE
USE Narran_Module
USE chltot_MODULE
USE NumNuj_Module
USE VFunc_Module
USE Storage_Module
USE Regins_Module
USE Energy_Module      ! PORT 2026-10-03: eref (energy reference)
USE Region_Module
USE Oops_Module
USE dip_MODULE
USE engpro_MODULE
USE TotalEng_Module, ONLY: ETot
USE Oops_Module
USE Masses_Module !, ONLY: SysReducedMass=>USys2, AtMass=>AtomicMass
USE InputFile_MODULE
USE parms_MODULE
USE qnumsp_Module
USE GaussQuady_Module
USE cherm1_Module
USE cherm2_Module
USE EDeriv_Module
USE Jtot_Module
USE Qall_Module
IMPLICIT NONE
SAVE
LOGICAL calc_kmat, more, there
LOGICAL, PARAMETER :: Full=.False.
CHARACTER(LEN=5) ename
CHARACTER(LEN=20) TypeOfMat, approx, proj, symm
INTEGER, PARAMETER :: inunit1 =80
INTEGER, PARAMETER :: inunit2 =81
INTEGER, PARAMETER :: inunit3 =82
INTEGER, PARAMETER :: outunit3=85
INTEGER, PARAMETER :: outunit4=86
INTEGER, PARAMETER :: outunit1=83
INTEGER, PARAMETER :: outunit2=84
INTEGER ichan,  nopen, ii, ithmat, lnblnk, njacobi, ndelves, IErr, i, nopen_last/-1/
INTEGER, ALLOCATABLE ::   elect(:),  nvib(:),   jrot(:)
INTEGER, ALLOCATABLE ::    lorb(:), chanl(:), nbasis(:)
REAL(KIND=WP_Kind), ALLOCATABLE :: energy(:), wvec_kmat(:)
REAL(KIND=WP_Kind), ALLOCATABLE ::   xksq_kmat(:),   xk(:)
REAL(KIND=WP_Kind), ALLOCATABLE ::  R_mat(:,:), RE_mat(:,:)
REAL(KIND=WP_Kind), ALLOCATABLE ::  K_mat(:,:), KE_mat(:,:)
REAL(KIND=WP_Kind), ALLOCATABLE ::  A_mat(:,:),  B_mat(:,:)
REAL(KIND=WP_Kind), ALLOCATABLE ::  E_mat(:,:),  F_mat(:,:)
REAL(KIND=WP_Kind), ALLOCATABLE :: AE_mat(:,:), BE_mat(:,:)
REAL(KIND=WP_Kind), ALLOCATABLE :: EE_mat(:,:), FE_mat(:,:)
REAL(KIND=WP_Kind), ALLOCATABLE ::   temp(:,:)
REAL(KIND=WP_Kind), ALLOCATABLE ::  temp1(:,:)
REAL(KIND=WP_Kind), ALLOCATABLE ::  temp2(:,:)
REAL(KIND=WP_Kind), ALLOCATABLE :: gamma(:), phor(:), phoi(:)

NAMELIST/kmatrix/ calc_kmat
!    calc_kmat  = .true. ==> Calculate K_Matrix.
!    EDeriv = .true. ==> Calculate derivative of K_Matrix with
!                           respect to the total energy.
WRITE(Out_Unit,*)'Called    KMatrix_Main'
NAtoms=3
!ALLOCATE(AtomicSymbol(NAtoms))  ! Atomic symbol
!ALLOCATE(Notes(NAtoms))         ! Notes
!ALLOCATE(AtomicNumber(NAtoms))  ! Atomic number
!ALLOCATE(MassNumber(NAtoms))    ! Mass number
!ALLOCATE(Abundance(NAtoms))     ! Isotopic abundance
!ALLOCATE(AtomicWeight(NAtoms))  ! Atomic Weight
!ALLOCATE(AtomicMass(NAtoms))    ! Atomic mass
CLOSE(Out_Unit)
!USys2=SysReducedMass
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/KMatrx_Main.txt',form='formatted', status='unknown')

irdindep = .false.
iwrindep = .true.
INQUIRE(file=InputDIR(1:LEN(TRIM(InputDIR)))//TRIM(InputFile),exist=there)
IF(there)THEN
   OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
   READ(In_Unit,NML=kmatrix, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'kmatrix')
      calc_kmat=.True.
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist kmatrix'
      WRITE(Msg_Unit,*)'ERROR with Namelist kmatrix'
      STOP 'Kmatrix_Main: kmatrix'
   ENDIF
   CLOSE(Unit=In_Unit, status='keep')
   WRITE(Out_Unit,NML=kmatrix)
ELSE
   WRITE(Msg_Unit,*)'Error: InputFile does not exist'
   STOP 'Stopping Inside Kmatrx_Main'
ENDIF

WRITE(*,*)"iregion=", iregion,"  ioops=", ioops, "Loc5 why Basis.bin?"
OPEN(Unit=bas_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Basis.bin', form='unformatted', status='unknown')

!-----------------------------------------------------------------------
! Establish all conversion factors and physical constants
!-----------------------------------------------------------------------
CALL constant
!-----------------------------------------------------------------------
! Read in all input data
!-----------------------------------------------------------------------
CALL readaph_ph
maxosc=maxhermt

ALLOCATE( elect(NVibRot),  nvib(NVibRot),   jrot(NVibRot) )
ALLOCATE(  lorb(NVibRot), chanl(NVibRot), nbasis(NVibRot) )
ALLOCATE(  R_mat(NVibRot, NVibRot), RE_mat(NVibRot, NVibRot) )
ALLOCATE(  K_mat(NVibRot, NVibRot), KE_mat(NVibRot, NVibRot) )
ALLOCATE(  A_mat(NVibRot, NVibRot),  B_mat(NVibRot, NVibRot) )
ALLOCATE(  E_mat(NVibRot, NVibRot),  F_mat(NVibRot, NVibRot) )
ALLOCATE( AE_mat(NVibRot, NVibRot), BE_mat(NVibRot, NVibRot) )
ALLOCATE( EE_mat(NVibRot, NVibRot), FE_mat(NVibRot, NVibRot) )
ALLOCATE(   temp(NVibRot, NVibRot) )
ALLOCATE(  temp1(NVibRot, NVibRot) )
ALLOCATE(  temp2(NVibRot, NVibRot) )
ALLOCATE(energy(NVibRot), wvec_kmat(NVibRot) )
ALLOCATE(  xksq_kmat(NVibRot),   xk(NVibRot) )
ALLOCATE( gamma(NVibRot), phor(NVibRot), phoi(NVibRot) )
IF(.not.ALLOCATED(chlast))ALLOCATE(chlast(maxosc*maxosc,0:mxl,narran))
IF(.not.ALLOCATED(numax))ALLOCATE(numax(nvibrot), numin(nvibrot))
IF(.not.ALLOCATED(xpth))ALLOCATE(xpth(maxhermt,narran),wpth(maxhermt,narran))

IF(.NOT.allocated(chinuj))ALLOCATE(chinuj(nvbrthrt))
IF(.not.ALLOCATED(chnow))ALLOCATE(chnow(maxosc*maxosc,0:mxl,narran))

IF(lphoto)THEN
   INQUIRE(file=OutDIR(1:LEN(TRIM(OutDIR)))//'R_matdph.bin', exist=there)
   IF(there)THEN
      OPEN(Unit=inunit3,File=OutDIR(1:LEN(TRIM(OutDIR)))//'R_Matdph.bin',form='unformatted', status='old')
   ELSE
      WRITE(Msg_Unit,*)'Error: No input file by that name'
      STOP 'main_kmatrix'
   ENDIF
   OPEN(Unit=outunit3,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Smatphr.bin', form='unformatted', status='unknown')
   OPEN(Unit=outunit4,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Smatphi.bin', form='unformatted', status='unknown')
ENDIF

rhoj = finish
! PORT (2026-10-11): when the long-range van der Waals propagation hook has
! carried the R-matrix out to APH3D_RF, the asymptotic matching must be
! evaluated consistently at that SAME hyperradius: rhoj drives strans (the
! Jacobi separation), the Morse coordinate map (yg), the Jacobian factor
! (bfaci) and the Bessel argument (z = k*rhoj*cos(theta_d)) inside
! Delvsf_New/upsiln.  Previously only the Bessel argument was moved to the
! larger radius (a hybrid convention) while every geometric factor stayed at
! `finish`.
BLOCK
   CHARACTER(LEN=32) :: env_rf_rhoj
   INTEGER :: stat_rf_rhoj
   REAL(Kind=WP_Kind) :: rf_rhoj
   CALL get_environment_variable('APH3D_RF', env_rf_rhoj, status=stat_rf_rhoj)
   IF (stat_rf_rhoj == 0) THEN
      READ(env_rf_rhoj, *) rf_rhoj
      IF (rf_rhoj > rhoj) rhoj = rf_rhoj
   END IF
END BLOCK
WRITE(Msg_Unit,*)'rhoj=',rhoj

!
! Read in URU_matrix
!

INQUIRE(file=OutDIR(1:LEN(TRIM(OutDIR)))//URU_File, exist=there)
IF(there)THEN
    OPEN(Unit=inunit1,File=OutDIR(1:LEN(TRIM(OutDIR)))//URU_File, form='unformatted', status='old')
ELSE
   WRITE(Msg_Unit,*)URU_File
   WRITE(Msg_Unit,*)'Error: No input file by that name'
   WRITE(Out_Unit,*)URU_File
   WRITE(Out_Unit,*)'Error: No input file by that name'
   STOP 'Stopping in Kmatrx_Main'
ENDIF

IF(EDeriv)THEN
   INQUIRE(file=OutDIR(1:LEN(TRIM(OutDIR)))//UREU_File(:LEN(TRIM(UREU_File))), exist=there)
   IF(there)THEN
      OPEN(Unit=inunit2,File=OutDIR(1:LEN(TRIM(OutDIR)))//UREU_File(:LEN(TRIM(UREU_File))), &
           form='unformatted',status='old')
   ELSE
      WRITE(Msg_Unit,*) URU_File(:lnblnk(URU_File))
      WRITE(Msg_Unit,*)'Error: No input file by that name'
      STOP 'Stopping in Kmatrx_Main'
   ENDIF
ENDIF

IF(calc_kmat)THEN
   OPEN(Unit=outunit1,File=OutDIR(1:LEN(TRIM(OutDIR)))//K_File(:lnblnk(K_File)), form='unformatted',status='unknown')
ENDIF
IF(EDeriv)THEN
   OPEN(Unit=outunit2,File=OutDIR(1:LEN(TRIM(OutDIR)))//KE_File(:lnblnk(KE_File)), form='unformatted', status='unknown')
ENDIF

WRITE(Msg_Unit,*)
WRITE(Msg_Unit,*)'Read from file= ',URU_File,' Write to file= ',K_File
WRITE(Out_Unit,*)'Read from file= ',URU_File,' Write to file= ',K_File

ithmat = 0
more = .true.
Write(Out_Unit,*)
Write(Out_Unit,*)
DO WHILE (more)
   iregion='delves '
   WRITE(Out_Unit,*)'Starting jth_Energy = ',ithmat+1
   WRITE(Out_Unit,*)'Total Energy = ', Etot
   WRITE(Out_Unit,'(A,ES21.14)')' usys2=',usys2
   WRITE(Out_Unit,*)'ndelves=',ndelves
   !WRITE(Msg_Unit,*)'Inside Kmatrx_Main Reading URU_Matrix'
   WRITE(Out_Unit,*)'Inside Kmatrx_Main Reading URU_Matrix'
   TypeOfMat = 'URU_Matrix'
   CALL Read_Matrix2(R_mat, energy, chanl, elect, nvib, jrot,      &
                     lorb, ndelves, Etot, ithmat, TypeOfMat, jtot, &
                     proj, approx, symm, inunit1, Out_Unit,        &
                     NVibRot, more)
   IF(EDeriv)THEN
      TypeOfMat = 'UREU_Matrix'
      ithmat = ithmat - 1
      CALL Read_Matrix2(RE_mat, energy, chanl, elect, nvib, jrot,      &
                        lorb, ndelves, Etot, ithmat, TypeOfMat, jtot,  &
                        proj, approx, symm, inunit2, Out_Unit,         &
                        NVibRot, more)
   ENDIF
   IF(.NOT.more)THEN
      ! PORT (2026-10-03): flush/close the K-matrix and URU files here.  They
      ! are written on outunit1 and read back later (Asymptotic_Analysis) on a
      ! different unit; without the CLOSE, files smaller than the Fortran I/O
      ! buffer (e.g. a 2x2 cold-collision K matrix) are still unflushed when
      ! the reader opens them and the read hits EOF.
      CLOSE(Unit=inunit1, status='keep')
      CLOSE(Unit=outunit1, status='keep')
      IF(EDeriv)THEN
         CLOSE(Unit=inunit2, status='keep')
         CLOSE(Unit=outunit2, status='keep')
      ENDIF
      !STOP 'complete'
      GOTO 991
   ENDIF
   !
   ! read in gamma vector
   !
   IF(lphoto)THEN
      READ(inunit3)(gamma(ii),ii=1,ndelves)
   ENDIF

   nopen = 0
   BLOCK
      CHARACTER(LEN=32) :: env_log
      INTEGER :: stat_log
      REAL(Kind=WP_Kind) :: eshift
      eshift = eref
      CALL get_environment_variable('APH3D_LOG_ENERGY', env_log, status=stat_log)
      IF (stat_log == 0 .AND. TRIM(ADJUSTL(env_log)) == '1') eshift = 0.0d0
      DO ichan = 1, ndelves
         xksq_kmat(ichan)=usys2*(Etot-(energy(ichan)-eshift))
         IF(Etot>(energy(ichan)-eshift))nopen = nopen + 1
         wvec_kmat(ichan)=SQRT(ABS(xksq_kmat(ichan)))
      ENDDO
   END BLOCK
   
   IF(Full)THEN
   WRITE(Out_Unit,*)'Inside Kmatrx_Main: The following are Delves states at last rho'
   Call Quant_Out (NDelves, chanl, elect, nvib, jrot, lorb, energy, xksq_kmat,'KMatrx_Main_1')
   ENDIF


   !--------------------------------------------------------------------
   ! Obtain the asymptotic Jacobi basis functions.
   !--------------------------------------------------------------------
   WRITE(Out_Unit,*)' G E N E R A T E    J A C O B I    B A S I S'

   ioops=.false.
   iregion='jacobi '
   WRITE(Out_Unit,*)'In Kmatrx_Main Calculating Jacobi Basis Functions'
   CALL JacBasis (nchanl, A_mat, B_mat, xksq3, wvec_kmat, endaph, chinuj, nnuj, njacobi)
   WRITE(Out_Unit,*)'completed JacBasis'
   BLOCK
      CHARACTER(LEN=32) :: env_log
      INTEGER :: stat_log
      REAL(Kind=WP_Kind) :: eshift
      eshift = eref
      CALL get_environment_variable('APH3D_LOG_ENERGY', env_log, status=stat_log)
      IF (stat_log == 0 .AND. TRIM(ADJUSTL(env_log)) == '1') eshift = 0.0d0
      IF(eshift/=0.d0)THEN
         DO ichan = 1, njacobi
            xksq3(ichan) = xksq3(ichan) + usys2*eshift
         ENDDO
         WRITE(Out_Unit,'(A,ES14.6)')' Kmatrx_Main: xksq3 shifted by usys2*eshift = ', usys2*eshift
      ENDIF
   END BLOCK
   
   nopen = 0
   DO ichan = 1, njacobi
      Energy3_Jacobi(ichan)=Etot-xksq3(ichan)/usys2
      IF(Etot>Energy3_Jacobi(ichan))nopen = nopen + 1
      chanl(ichan) = kchan(ichan)
      elect(ichan) = elect3(ichan)
      nvib(ichan)  = mvib3(ichan)
      jrot(ichan)  = jrot3(ichan)
      lorb(ichan)  = lorb3(ichan)
   ENDDO
   
   IF(Full)THEN
   Write(Out_Unit,*)"Nopen_Max=",Nopen
   WRITE(Out_Unit,*)'In Kmatrx_Main After Calling JacBasis: NJacobi=', NJacobi
   CALL Quant_Out(NJacobi, kchan, elect3, mvib3, jrot3, lorb3, energy3_Jacobi, xksq3,'KMatrx_Main_2')
   ENDIF
   
   !
   ! Calculate the K_Matrix.
   !
   WRITE(Out_Unit,*)'jth_Energy = ',ithmat
   WRITE(Out_Unit,*)'Total Energy = ', Etot
   WRITE(Out_Unit,'(A,ES21.14)')' usys2=',usys2
   WRITE(Out_Unit,*)'ndelves=',ndelves
   kenergy = kenergy + 1

   WRITE(Out_Unit,*)'nopen=',nopen
   IF(FULL)THEN
      WRITE(Out_Unit,*)'In Kmatrx_Main Just Before Calling Kmatrx_Calc: NJacobi=', NJacobi
      CALL Quant_Out(ndelves, chanl, elect, nvib, jrot, lorb, energy3_Jacobi, xksq3,'KMatrx_Main_3')
   ENDIF
   
   CALL Kmatrx_Calc (ndelves, R_mat, RE_mat, K_mat, KE_mat, xk,  &
                 chinuj, wvec_kmat, nbasis, usys2, Etot, jtot,    &
                 chanl, elect, nvib, jrot, lorb, xksq3,        &
                 energy3_Jacobi, nchanl, narran, temp, endaph,       &
                 Out_Unit, A_mat, B_mat, E_mat, F_mat,       &
                 AE_mat, BE_mat, EE_mat, FE_mat,             &
                 nopen, calc_kmat, EDeriv,njacobi,           &
                 gamma, phor, phoi, temp1,temp2)
   
   REWIND bas_unit
   irdindep = .false.
   iwrindep = .false.
   !WRITE(Msg_Unit,*)'Kmatrx_Main writing'
   WRITE(Out_Unit,*)'Kmatrx_Main writing'
   IF(calc_kmat)THEN
      TypeOfMat='K_Matrix'
      IF(.not.ALLOCATED(QNumbr))ALLOCATE(QNumbr(NQuant,njacobi))
      IF(.not.ALLOCATED(QLabels))ALLOCATE(QLabels(NQuant))
      IF(.not.ALLOCATED(ConLabels))ALLOCATE(ConLabels(NConserved))
      IF(.not.ALLOCATED(Conserved))ALLOCATE(Conserved(NConserved))
      IF(.not.ALLOCATED(Asymptenergies))ALLOCATE(Asymptenergies(njacobi))
      IF(.not.ALLOCATED(wvec))ALLOCATE(wvec(njacobi))
      IF(.not.ALLOCATED(xksq))ALLOCATE(xksq(njacobi))
      QLabels=['chanl','elect','nvib ','jrot ','lorb ']
      ConLabels=['Jtot  ','Parity','C2v   ','A1    ']
      Conserved=[jtot,parity,2,1]
      DO i=1,njacobi
          Asymptenergies(i)=energy3_Jacobi(i)
          xksq(i)=usys2*(Etot-energy3_Jacobi(i))
          wvec(i)=SQRT(ABS(xksq(i)))   !tmpmod gregparker abs shouldn't be necessary...
          QNumbr(1,i)=chanl(i)
          QNumbr(2,i)=elect(i)
          QNumbr(3,i)=nvib(i)
          QNumbr(4,i)=jrot(i)
          QNumbr(5,i)=lorb(i)
      ENDDO
      CALL Write_Matrix2(K_Mat, nopen, etot, ithmat, TypeOfMat, outunit1)
      WRITE(Out_Unit,*)'K_Matrix Open-Open Block'
      CALL MxOut(K_Mat, nopen, nopen)
   ENDIF
   IF(EDeriv)THEN
      TypeOfMat='KE_Matrix'
      CALL Write_Matrix2(K_Mat, nopen, etot, ithmat, TypeOfMat, outunit2)
      WRITE(Out_Unit,*)'Energy Derivative of the K_Matrix'
      CALL MxOut(KE_Mat, nopen, nopen)
   ENDIF
   IF(nopen/=nopen_last)THEN
      WRITE(0,*)'Ithmat=',ithmat,'  nopen=',nopen
      nopen_last=nopen
   ENDIF
   !DEALLOCATE(QNumbr) 
   !DEALLOCATE(QLabels)
   DEALLOCATE(ConLabels)
   !DEALLOCATE(Conserved)
   !DEALLOCATE(Asymptenergies)
   !DEALLOCATE(wvec)
   !DEALLOCATE(xksq)
   IF(lphoto)THEN
      !
      ! scale these by 1/sqrt(2mu) - normalization
      !
      DO ii=1,nopen
         phor(ii)=phor(ii)/sqrt(usys2)
         phoi(ii)=phoi(ii)/sqrt(usys2)
      ENDDO
      WRITE(outunit3)(phor(ii),ii=1,nopen)
      WRITE(Out_Unit,*)'phor',(phor(ii),ii=1,nopen)
      WRITE(outunit4)(phoi(ii),ii=1,nopen)
      WRITE(Out_Unit,*)'phoi',(phoi(ii),ii=1,nopen)
   ENDIF
ENDDO

991 CONTINUE
    
OPEN(Unit=NOpenMAX_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/NOpenMax.bin', form='unformatted', status='unknown')
WRITE(NOpenMAX_Unit)NOpen   !"At this point MOpenMax=Nopen"
CLOSE(NOpenMAX_Unit)


!DEALLOCATE(numax,numin)
!DEALLOCATE(xpth)
!DEALLOCATE(wpth)
!DEALLOCATE(chinuj)
!DEALLOCATE(chlast)
!DEALLOCATE(chnow)
DEALLOCATE(R_mat)
DEALLOCATE(RE_mat)
DEALLOCATE(K_mat)
DEALLOCATE(KE_mat)
!DEALLOCATE(energy)
!DEALLOCATE(wvec_kmat)
!DEALLOCATE(temp)
DEALLOCATE(A_mat)
DEALLOCATE(B_mat)
DEALLOCATE(E_mat)
DEALLOCATE(F_mat)
!DEALLOCATE(gamma)
!DEALLOCATE(phor)
!DEALLOCATE(phoi)
!DEALLOCATE(temp1)
!DEALLOCATE(temp2)
DEALLOCATE(AE_mat)
DEALLOCATE(BE_mat)
DEALLOCATE(EE_mat)
DEALLOCATE(FE_mat)
!DEALLOCATE(xksq)
!DEALLOCATE(xk)

!DEALLOCATE(AtomicSymbol)
!DEALLOCATE(Notes)
!DEALLOCATE(AtomicNumber)
!DEALLOCATE(MassNumber)
!DEALLOCATE(Abundance)
!DEALLOCATE(AtomicWeight)
!DEALLOCATE(AtomicMass)


WRITE(Out_Unit,*)'K_Matrix Calculation Complete'
CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/BasisForScatt.txt',Form='FORMATTED',ACCESS='APPEND')
WRITE(Out_Unit,*)'Completed Kmatrx_Main'
WRITE(Out_Unit,*)

RETURN
ENDSUBROUTINE Kmatrx_Main
