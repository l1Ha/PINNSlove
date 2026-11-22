SUBROUTINE UMatPh
!
! This program transforms the R_Matrix in the APH basis to the
! R_Matrix in the delves basis.
!
! Modified by M. Braunstein, Spring 1991, to include photodissociation
! Transforms Gamma_APH matrix from aph to delves (half transform)
!
USE numeric_kinds_Module
USE dip_Module
USE FileUnits_Module
USE FileNames_Module
USE InputFile_Module
USE Energy_Module
USE EDeriv_Module
USE Masses_Module
USE TotalEng_Module
USE APHchl_Module, ONLY: Neigmin_Sum, Neigmin_Save
IMPLICIT NONE
INTEGER naph, ndelves, IErr
INTEGER i, ithmat, n, lnblnk, jtot, nchanl, kenergy, DtValues(8)
INTEGER, PARAMETER:: inunit1=80
INTEGER, PARAMETER:: inunit2=81
INTEGER, PARAMETER:: inunit4=92
INTEGER, PARAMETER:: outunit1=83
INTEGER, PARAMETER:: outunit2=84
INTEGER, PARAMETER:: outunit3=85
INTEGER, ALLOCATABLE :: IDummy(:)
INTEGER, ALLOCATABLE :: elect3(:), nvib3(:), jrot3(:), lorb3(:), chanl3(:)
INTEGER lam_n(0:1000), lamfst, lamlst
LOGICAL calc_rmat, more, there
CHARACTER(len=5) CurZone
CHARACTER(len=8) Today
CHARACTER(len=10)Hour
CHARACTER(len=6) iprt
CHARACTER(len=20) TypeOfRmat, approx, proj, symm
CHARACTER(len=20) R_file
CHARACTER(len=21) RE_file
CHARACTER(LEN=200) File1, File2
CHARACTER(len=5) ename
REAL(Kind=WP_Kind), ALLOCATABLE ::  Energy_APH(:), R_mat_APH(:,:), RE_mat_APH(:,:),  umat(:, :)
REAL(Kind=WP_Kind), ALLOCATABLE ::  Energy_Delves(:), R_mat_Delves(:,:), RE_mat_Delves(:,:)
REAL(Kind=WP_Kind), ALLOCATABLE ::  energy3(:), Gamma_APH(:), Gamma_Delves(:), xksq3(:)
NAMELIST/umat_trans/ calc_rmat, iprt
NAMELIST/basis_n/ lam_n, lamfst, lamlst
NAMELIST/dipole2/ jf,prf,lamf,mf,jd,prd,lamd,md,ji,pri,lami, mi,nsfuncu,nsfuncl,lphoto,iaxis,enddip

WRITE(Out_Unit,*)'Called    Umatph'
CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/UMatph.txt',form='formatted',status='unknown')

CALL Date_And_Time(today, hour, curzone, dtvalues)
CALL Print_Date_And_Time(Today, Hour, CurZone, DtValues, Out_Unit)
! 
! Read in Data file to determine input and output file names.
! Read in Options.
!    calc_rmat  = .true. ==> Transform R_Matrix.
!    EDeriv = .true. ==> Transform derivative of R_Matrix with
!                           respect to the total energy.
lphoto = .false.
INQUIRE (file=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,exist=there)
IF(there)THEN
   OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
   READ(In_Unit,NML=Totenergy)
   CLOSE(Unit=In_Unit, status='keep')
   WRITE(Out_Unit,NML=TotEnergy)
   
   OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
   READ(In_Unit,NML=umat_trans, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'umat_trans')
      Calc_RMat=.true.
      Iprt='medium'
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist umat_trans'
      WRITE(Msg_Unit,*)'ERROR with Namelist umat_trans'
      STOP 'UMatPh: umat_trans'
   ENDIF
   CLOSE(Unit=In_Unit, status='keep')
   WRITE(Out_Unit,NML=umat_trans)
ELSE
   WRITE(Msg_Unit,*)'Error: umat_trans does not exist'
   STOP 'Umatph: umat_trans'
ENDIF

INQUIRE (file=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,exist=there)
IF(there)THEN
   Lam_N=0
   OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
   READ(In_Unit,NML=basis_n,IOSTAT=IERR)
   IF(IERR==-1)THEN
      WRITE(Out_Unit,*)'ERROR missing Namelist Basis_n'
      WRITE(Msg_Unit,*)'ERROR missing Namelist Basis_n'
      STOP 'UmatPh: Basis_n'
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Basis_n'
      WRITE(Msg_Unit,*)'ERROR with Namelist Basis_n'
      STOP 'UMatPh: Basis_n'
   ENDIF
   WRITE(Out_Unit,NML=basis_n)
   CLOSE(Unit=In_Unit)
   lam_N=Neigmin_Save
   NChanl=SUM(Lam_N)
   NAPH=NChanl
ELSE
   WRITE(Msg_Unit,*)'Error: InputFile does not exist', InputDIR(1:LEN(TRIM(InputDIR)))//InputFile
   WRITE(Out_Unit,*)'Error: InputFile does not exist', InputDIR(1:LEN(TRIM(InputDIR)))//InputFile
   STOP 'Umatph: Basis_n'
ENDIF

!
INQUIRE (file=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,exist=there)
IF(there)THEN
   OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
   READ(In_Unit,NML=dipole2,IOSTAT=IERR)
   IF(IERR==-1)THEN
      LPhoto=.False.
      CALL NameList_Default(Out_Unit, 'Dipole2')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Dipole2'
      WRITE(Msg_Unit,*)'ERROR with Namelist Dipole2'
      STOP 'UMatPh: Dipole2'
   ENDIF
   CLOSE(Unit=In_Unit, status='keep')
   WRITE(Out_Unit,NML=dipole2)
ELSE
   WRITE(Msg_Unit,*)'Error: InputFile does not exist', InputDIR(1:LEN(TRIM(InputDIR)))//InputFile
   WRITE(Out_Unit,*)'Error: InputFile does not exist', InputDIR(1:LEN(TRIM(InputDIR)))//InputFile
   STOP 'Umatph: Dipole2'
ENDIF
!
! Open input and output files.
!
INQUIRE (file=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Umat.bin', exist=there)
IF(there)THEN
   OPEN(Unit=umat_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Umat.bin', form='unformatted',status='old')
   READ(umat_unit)naph, ndelves
   WRITE(Msg_Unit,*)'naph=',naph,'   ndelves=',ndelves
   WRITE(Out_Unit,*)'naph=',naph,'   ndelves=',ndelves
   IF(naph>nchanl)THEN
      naph = nchanl
      WRITE(Out_Unit,*)'Reducing naph'
      WRITE(Msg_Unit,*)'naph=',naph,'   ndelves=',ndelves
      WRITE(Out_Unit,*)'naph=',naph,'   ndelves=',ndelves
   ELSEIF(naph==nchanl)THEN
      WRITE(Msg_Unit,*)'naph=',naph,'   ndelves=',ndelves
      WRITE(Out_Unit,*)'naph=',naph,'   ndelves=',ndelves  
   ELSE
      WRITE(Out_Unit,*)'nchanl too large',nchanl,' naph=',naph
      STOP 'UmatPh: nchanl too large'
   ENDIF


!---------------------------------------------------------------
   ALLOCATE(R_mat_APH(naph,naph))
   ALLOCATE(R_mat_Delves(ndelves,ndelves))
   ALLOCATE(umat(ndelves, naph))
   ALLOCATE(Gamma_APH(naph))
   ALLOCATE(Gamma_Delves(ndelves))
   ALLOCATE(IDummy(naph))
   ALLOCATE(Energy_APH(naph))
   ALLOCATE(xksq3(ndelves))
   ALLOCATE(energy3(ndelves))
   ALLOCATE(elect3(ndelves))
   ALLOCATE(nvib3(ndelves))
   ALLOCATE(jrot3(ndelves))
   ALLOCATE(lorb3(ndelves))
   ALLOCATE(chanl3(ndelves))
   IF(EDeriv)THEN
      ALLOCATE(RE_mat_APH(naph,naph))
      ALLOCATE(RE_mat_Delves(ndelves,ndelves))
   ENDIF
!---------------------------------------------------------------

   CALL Read_Umat(naph, ndelves, umat, chanl3, elect3, nvib3, jrot3, lorb3, energy3, umat_unit)
   WRITE(Out_Unit,*)'U_Matrix'
   CALL MxOut(umat, ndelves, naph)
   CLOSE(umat_unit)
ELSE
   WRITE(Msg_Unit,*)'Error:  1 No input file by that name'
   STOP 'Umatph: No input file by that name'
ENDIF


OPEN(Unit=Rmat_APH_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Rmat_APH.bin',  form='unformatted',status='unknown')
OPEN(Unit=REmat_APH_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/REmat_APH.bin', form='unformatted',status='unknown')
WRITE(Msg_Unit,*)'Read from= ',R_File,'  Write to= ',URU_File
WRITE(Out_Unit,*)'Read from= ',R_File,'  Write to= ',URU_File
DO kenergy=1,nenergy
   IF(kenergy<10)THEN
         WRITE(ename,'(i1)') kenergy
      ELSEIF(kenergy<100)THEN
         WRITE(ename,'(i2)') kenergy
      ELSEIF(kenergy<1000)THEN
         WRITE(ename,'(i3)') kenergy
      ELSEIF(kenergy<10000)THEN
         WRITE(ename,'(i4)') kenergy
      ELSEIF(kenergy<100000)THEN
         WRITE(ename,'(i5)') kenergy
      ELSE
         WRITE(Msg_Unit,*)'kenergy to large=',kenergy
         STOP 'Umatph: kenergy to large'
      ENDIF

      R_file = 'working/R'//ename(:lnblnk(ename))
      RE_file = 'working/RE'//ename(:lnblnk(ename))

IF(calc_rmat)THEN
   File1=OutDIR(1:LEN(TRIM(OutDIR)))//R_File
   INQUIRE (file=OutDIR(1:LEN(TRIM(OutDIR)))//R_File, exist=there)
   IF(there)THEN
      OPEN(Unit=inunit1,File=OutDIR(1:LEN(TRIM(OutDIR)))//R_File, form='unformatted', status='old')
   ELSE
      WRITE(Msg_Unit,*)'Error: 2 No input file by that name'
      STOP 'Umatph: No input file by that name'
   ENDIF
ENDIF

IF(EDeriv)THEN
   File2=OutDIR(1:LEN(TRIM(OutDIR)))//RE_File(1:LEN(TRIM(R_File)))
   INQUIRE (file=OutDIR(1:LEN(TRIM(OutDIR)))//RE_File(1:LEN(TRIM(RE_File))), exist=there)
   IF(there)THEN
      OPEN(Unit=inunit2,File=OutDIR(1:LEN(TRIM(OutDIR)))//RE_File(1:LEN(TRIM(RE_File))),form='unformatted', status='old')
   ELSE
      WRITE(Msg_Unit,*)'Error:  3 No input file by that name'
      STOP 'Umatph: No input file by that name'
   ENDIF
ENDIF
!
IF(lphoto)THEN
   WRITE(Out_Unit,*)'Doing photodissociation '
   INQUIRE (file=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/R_Matrixph.bin', exist=there)
   IF(there)THEN
      OPEN(Unit=inunit4,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/R_Matrixph.bin', form='unformatted', status='old')
   ELSE
      WRITE(Msg_Unit,*)'Error: 4 No input file by that name'
      STOP 'Umatph: No input file by that name'
   ENDIF
ENDIF
!

IF(calc_rmat.AND.kenergy==1)THEN
   OPEN(Unit=outunit1,File=OutDIR(1:LEN(TRIM(OutDIR)))//URU_File(:lnblnk(URU_File)),form='unformatted',status='REPLACE')
ENDIF

IF(EDeriv)THEN
   OPEN(Unit=outunit2,File=OutDIR(1:LEN(TRIM(OutDIR)))//UREU_File(:lnblnk(UREU_File)), form='unformatted', status='unknown')
ENDIF

!
IF(lphoto)THEN
   OPEN(Unit=outunit3,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/R_matdph.bin', form='unformatted', status='unknown')
ENDIF


!
! Read in R_matrix
!
IF(kenergy==1)ithmat = 0
ithmat=kenergy-1
more = .true.
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)
DO WHILE (more)
   WRITE(Out_Unit,*)'umatph: Reading APH R_Matrix'
   IF(calc_rmat)THEN
      TypeOfRmat = 'R_Matrix'
      CALL Read_Matrix2(R_mat_APH, Energy_APH, IDummy, IDummy, IDummy, IDummy, IDummy, n, Etot, ithmat, TypeOfRmat, jtot, &
                        proj, approx, symm, inunit1, Out_Unit, naph, more)
      WRITE(Out_Unit,'(2A3, i5, A,i5,2A)')'P2', " n=", n, " ithmat=", ithmat, " TypeOfRmat=", TypeOfRmat
      WRITE(Out_Unit,'(A,I5,7A,I5,A,I5,A,L)')" Jtot=", jtot," proj=", proj, " approx=", approx, " symm=", symm," inunit1=", inunit1, " out_unit=",Out_Unit," more=",more
   ENDIF

   IF(EDeriv)THEN
      TypeOfRmat = 'RE_Matrix'
      ithmat = ithmat - 1
      CALL Read_Matrix2(RE_mat_APH, Energy_APH, IDummy, IDummy, IDummy, IDummy, IDummy, n, Etot, ithmat, TypeOfRmat, &
                        jtot, proj, approx, symm, inunit2, Out_Unit, naph, more)
      WRITE(Out_Unit,'(2A3, i5, A,i5,2A)')'P4', " n=", n, " ithmat=", ithmat, " TypeOfRmat=", TypeOfRmat
      WRITE(Out_Unit,'(A,I5,7A,I5,A,I5,A,L)')" Jtot=", jtot," proj=", proj, " approx=", approx, " symm=", symm," inunit1=", inunit1, " out_unit=",Out_Unit," more=",more
   ENDIF

   IF(.NOT.more)THEN
      IF(calc_rmat)THEN
         CLOSE(Unit=inunit1, status='keep')
         IF(kenergy==nenergy)CLOSE(Unit=outunit1, status='keep')
      ENDIF
      IF(EDeriv)THEN
         CLOSE(Unit=inunit2, status='keep')
         IF(kenergy==nenergy)CLOSE(Unit=outunit2, status='keep')
      ENDIF
!              STOP 'complete'
         GOTO 990
   ENDIF
!------------------------------------------------
! Read Gamma_APH matrix
!
   IF(lphoto)THEN
         READ(inunit4)(Gamma_APH(i),i=1,naph)
         WRITE(Out_Unit,*)'gamma1',(Gamma_APH(i),i=1,naph)
   ENDIF

!
! Transform the R_Matrix.
!
   WRITE(Out_Unit,*)'jth_Energy = ',ithmat
   WRITE(Out_Unit,*)'Total Energy = ', Etot

   IF(calc_rmat)THEN
      IF(iprt=='medium')THEN
         WRITE(Out_Unit,*)'R_Matrix(APH)'
         CALL MxOut(R_mat_APH, naph, naph)
      ELSEIF(iprt=='full')THEN
         WRITE(Out_Unit,*)'R_Matrix(APH)'
         CALL MxOutD(R_mat_APH, naph, naph, 0)
      ENDIF
      
      !GregParker  TEMPMOD ?
      TypeOfRmat='R_Mat_APH'
      CALL Write_Matrix2Old(R_mat_APH, Energy_APH, IDummy, IDummy, IDummy, IDummy, IDummy, n, Etot, ithmat, &
                         TypeOfRmat, jtot, proj, approx, symm, Rmat_APH_Unit)
      
      CALL transf(umat, R_mat_APH, R_mat_Delves, naph, ndelves)
      
      IF(iprt=='medium')THEN 
         WRITE(Out_Unit,*)'Transformed R_Matrix(Delves)' 
         CALL MxOut(R_mat_Delves, ndelves, ndelves) 
      ELSEIF(iprt=='full')THEN 
         WRITE(Out_Unit,*)'Transformed R_Matrix(Delves)'
         CALL MxOutD(R_mat_Delves, ndelves, ndelves, 0) 
      ENDIF 
   ENDIF

   IF(EDeriv)THEN
      IF(iprt=='medium')THEN 
         WRITE(Out_Unit,*)'Energy derivative of the R_Matrix(APH)' 
         CALL MxOut(RE_mat_APH, naph, naph) 
      ELSEIF(iprt=='full')THEN 
         WRITE(Out_Unit,*)'Energy derivative of the R_Matrix(APH)'
         CALL MxOutD(RE_mat_APH, naph, naph, 0) 
      ENDIF
      !GregParker  TEMPMOD?
      TypeOfRmat='RE_Mat_APH'
      CALL Write_Matrix2Old(RE_mat_APH, Energy_APH, IDummy, IDummy, IDummy, IDummy, IDummy, n, Etot, ithmat, &
                         TypeOfRmat, jtot, proj, approx, symm, REmat_APH_Unit)

      CALL transf(umat, RE_mat_APH, RE_mat_Delves, naph, ndelves)
      
      IF(iprt=='medium')THEN 
         WRITE(Out_Unit,*)'Transformed Energy derivative ', 'of the R_Matrix(Delves)' 
         CALL MxOut(RE_mat_Delves, ndelves, ndelves) 
      ELSEIF(iprt=='full')THEN 
         WRITE(Out_Unit,*)'Transformed Energy derivative ', 'of the R_Matrix(Delves)'
         CALL MxOutD(RE_mat_Delves, ndelves, ndelves, 0) 
      ENDIF 
   ENDIF
!-------------------------------------------------
! Transfrom Gamma_APH matrix from aph to delves
!
   IF(lphoto)THEN
      CALL tgam(ndelves, naph, umat, Gamma_APH, Gamma_Delves)
   ENDIF
!--------------------------------------------------


   DO i=1,ndelves
      xksq3(i)=usys2*(Etot-energy3(i))
   ENDDO

!
! Write Matrices for asymptotic analysis.
!
   WRITE(Out_Unit,*)'umatph: Writing Delves R_matrix'
   IF(calc_rmat)THEN
      IF(iprt/='none  ')THEN
      CALL Quant_Out(ndelves, chanl3, elect3, nvib3, jrot3, lorb3, energy3, xksq3,'UMatph_1')
      ENDIF
      TypeOfRmat='URU_Matrix'
      WRITE(Out_Unit,*)'ndelves=',ndelves
      CALL Write_Matrix2Old(R_mat_Delves, energy3, chanl3, elect3, nvib3, jrot3, lorb3, ndelves, Etot, ithmat, &
                         TypeOfRmat, jtot, proj, approx, symm, outunit1)
   ENDIF

   IF(EDeriv)THEN
      TypeOfRmat='UREU_Matrix'
      CALL Write_Matrix2Old(RE_mat_Delves, energy3, chanl3, elect3, nvib3, jrot3, lorb3, ndelves, Etot, ithmat, &
                         TypeOfRmat, jtot, proj, approx, symm, outunit2)
   ENDIF
!
   IF(lphoto)THEN
      WRITE(outunit3)(Gamma_Delves(i),i=1,ndelves)
   ENDIF
   
   WRITE(Out_Unit,*) 'i, chanl3(i), elect3(i), nvib3(i), jrot3(i), lorb3(i)'
   DO i=1,ndelves
     WRITE(Out_Unit,'(6i3)')i, chanl3(i), elect3(i), nvib3(i), jrot3(i), lorb3(i)
   ENDDO

ENDDO
!---------------------------------------------------------------


990 CONTINUE
ENDDO


DEALLOCATE(R_mat_APH,R_Mat_Delves,Energy_APH, umat, energy3)
DEALLOCATE(Gamma_APH,Gamma_Delves)
DEALLOCATE(IDummy, elect3, nvib3, jrot3, lorb3, chanl3)
IF(EDeriv)DEALLOCATE(RE_mat_APH, RE_mat_Delves)

!---------------------------------------------------------------
WRITE(Out_Unit,*)'URU Transformation is now complete'
CALL Date_And_Time(today, hour, curzone, dtvalues)
CALL Print_Date_And_Time(Today, Hour, CurZone, DtValues, Out_Unit)

CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/BasisForScatt.txt',Form='FORMATTED',ACCESS='APPEND')
WRITE(Out_Unit,*)'Completed Umatph'
WRITE(Out_Unit,*)

RETURN
ENDSUBROUTINE UMatph

SUBROUTINE tgam(ndelves, naph, umat, Gamma_APH, Gamma_Delves)
!this routine is called by:
!         Umatphph
!this routine calls
!         
!-----------------------------------------------------------------------
USE fileunits_Module
USE numeric_kinds_Module
IMPLICIT NONE

INTEGER ndelves, naph, idelve, iaph
REAL(Kind=WP_Kind) umat(ndelves,naph), Gamma_APH(naph), Gamma_Delves(ndelves),  sum
WRITE(Out_Unit,*)'begin tgam'
DO  idelve=1,ndelves
   sum=0.0
   DO  iaph=1,naph
      sum=sum+umat(idelve,iaph)*Gamma_APH(iaph)
!           WRITE(Out_Unit,*)'gamma2',umat(idelve,iaph),Gamma_APH(iaph)
ENDDO
   Gamma_Delves(idelve)=sum
   IF(idelve<=10)THEN
      WRITE(Out_Unit,*)'gamma2',naph,idelve,Gamma_Delves(idelve)
   ENDIF
ENDDO
WRITE(Out_Unit,*)'finish tgam'

RETURN
ENDSUBROUTINE Tgam
