SUBROUTINE READINFBR
!------------------------------------------------
!READ in DVR parameter
!------------------------------------------------
USE DVRNL_Module
USE InputFile_Module
USE System_Module
USE Dipole_Module
USE FileUnits_Module
USE PES_MODULE
INTEGER i
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,     form='formatted',   status='old')

!
!  -------------------------------------------------------------------
!     READ in input
!  -------------------------------------------------------------------
!
!
READ(In_Unit,NML=rstart)
WRITE(Out_Unit,NML=rstart)
READ(In_Unit,NML=dipole)
WRITE(Out_Unit,NML=dipole)
READ(In_Unit,NML=system)
WRITE(Out_Unit,NML=system)
READ(In_Unit,NML=sysu)
WRITE(Out_Unit,NML=sysu)


!
 1004 FORMAT(/'upper state input echo'/)
  995 FORMAT(/'lower state input echo'/)
  996 FORMAT(/'input echo data for both states'/)
   11 FORMAT(a25)
!
!  -------------------------------------------------------------------
! WRITE out input to unit 6
!  -------------------------------------------------------------------
!
WRITE(Out_Unit,rstart)
WRITE(Out_Unit,996)
WRITE(Out_Unit,'(a,3(e22.14))') ' mass = ',(mass(i),i=1,3)
WRITE(Out_Unit,'(a,3(i22))') ' nthesec = ',(nthesec(i),i=1,Num_DVR_Sectors)
WRITE(Out_Unit,'(a,3(i22))') ' nchisec = ',(nchisec(i),i=1,Num_DVR_Sectors)
WRITE(Out_Unit,'(a,3(e22.14))') ' rendsec = ',(rendsec(i),i=1,Num_DVR_Sectors)
WRITE(Out_Unit,*)  'peigv = ',peigv
WRITE(Out_Unit,*)  'debug = ',debug
WRITE(Out_Unit,*)  'medium = ',medium
WRITE(Out_Unit,'(a,3(e22.14))') ' potmsec = ',(potmsec(i),i=1,Num_DVR_Sectors)
WRITE(Out_Unit,'(a,3(e22.14))') ' pevmsec = ',(pevmsec(i),i=1,Num_DVR_Sectors)
WRITE(Out_Unit,'(a,3(i22))') ' ltheta = ',ltheta
WRITE(Out_Unit,'(a,3(e22.14))') ' rho1 = ',rho1
WRITE(Out_Unit,*)  'matelem = ',matelem
WRITE(Out_Unit,*)  'h3sys = ',h3sys
WRITE(Out_Unit,'(a,3(e22.14))') ' rhofbr1 = ',rhofbr1
WRITE(Out_Unit,'(a,3(e22.14))') ' rhofbr2 = ',rhofbr2
WRITE(Out_Unit,'(a,3(i22))') ' nfbr = ',nfbr
WRITE(Out_Unit,'(a,3(i22))') ' nsfunc = ',nsfunc
WRITE(Out_Unit,'(a,3(i22))') ' nrow = ',nrow
WRITE(Out_Unit,'(a,3(i22))') ' ncol = ',ncol
WRITE(Out_Unit,'(a,3(i22))') ' ngood1 = ',ngood1
WRITE(Out_Unit,'(a,3(i22))') ' ngood2 = ',ngood2
WRITE(Out_Unit,'(a,3(i22))') ' Num_DVR_Sectors = ',Num_DVR_Sectors
WRITE(Out_Unit,'(a,3(i22))') ' initsec = ',initsec
WRITE(Out_Unit,'(a,3(i22))') ' nrhomx = ',nrhomx
WRITE(Out_Unit,'(a,3(e22.14))') ' fmpotma = ',fmpotma
WRITE(Out_Unit,'(a,3(e22.14))') ' deltar = ',deltar
WRITE(Out_Unit,*)  ' logspace = ',logspace
!
WRITE(Out_Unit,1004)
WRITE(Out_Unit,'(a)')pesnamu
WRITE(Out_Unit,*)  'jevenu = ',jevenu
WRITE(Out_Unit,'(a,3(i22))') ' nrhos1u = ', nrhos1u
WRITE(Out_Unit,'(a,3(e22.14))') ' rhos1u = ', (rhos1u(i),i=1,nrhos1u)
WRITE(Out_Unit,'(a,3(e22.14))') ' pevms1u = ',(pevms1u(i),i=1,nrhos1u)
WRITE(Out_Unit,'(a,3(i22))') ' iparu = ',iparu
WRITE(Out_Unit,'(a,3(i22))') ' jtotu= ',jtotu
WRITE(Out_Unit,'(a,3(i22))') ' megau= ',megau
WRITE(Out_Unit,*) ' cstestu = ',cstestu
WRITE(Out_Unit,*)  ' pointu = ',pointu

CLOSE(Unit=In_Unit)


RETURN
END
