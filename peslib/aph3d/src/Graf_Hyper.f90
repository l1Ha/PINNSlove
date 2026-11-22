SUBROUTINE Graf_Hyper
USE Numbers_Module
USE FileUnits_Module
USE Masses_Module
USE MassFactor2_Module
USE PES_MODULE
USE InputFile_Module
USE Narran_Module
USE BasisRegions_Module
USE Contour_Module
USE Convrsns_Module
!
!  program gen_contour.
!
IMPLICIT NONE
LOGICAL :: there
LOGICAL, PARAMETER:: debg=.False.
!
!  calculates potential over a grid of points to be used in making
!  contour plots.
!  this version calls surface SUBROUTINE in same FORMAT as aph3d, and
!  needs it loaded as a binary.
!
CHARACTER(LEN=4) SecLabel
CHARACTER(LEN=35) VeeLabel
CHARACTER(LEN=36) CsvLabel, MTVLabel
CHARACTER(LEN=80) SubTitle
INTEGER :: i, j, k, nchi, ISector, IERR
REAL(Kind=WP_Kind) :: chi(3)
REAL(Kind=WP_Kind) :: xi, dx, yi, dy, rr, theta, sinth, vv
REAL(Kind=WP_Kind) :: chii, xp (3)
REAL(Kind=WP_Kind), PARAMETER:: PhiChi=1.d0
REAL(KIND=WP_Kind), ALLOCATABLE::v (:, :)
REAL :: xx(361), yy(361)


WRITE(Out_Unit,*)'Called    Graf_Hyper'
CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/Graf_Hyper.txt',Form='FORMATTED')

DO ISector=1,NSectors
rho=Basis_Dist(ISector)
vmax = Ten**(RANGE(One)-10)
IF(ISector<10)THEN
   WRITE(SecLabel,'(4I1)')0,0,0,ISector
ELSEIF(ISector<100)THEN
   WRITE(SecLabel,'(2I1,I2)')0,0,ISector
ELSEIF(ISector<1000)THEN
   WRITE(SecLabel,'(I1,I3)')0,ISector
ELSEIF(ISector<10000)THEN
   WRITE(SecLabel,'(I4)')ISector
ELSE
   WRITE(Out_Unit,*)'ISector is greater than 9999:',ISector
   STOP 'Graf_Hyper'
ENDIF
VeeLabel='PlotVee/Hyper/VeeAtSector'//SecLabel//'.wg'
CsvLabel='PlotVee/Hyper/VeeAtSector'//SecLabel//'.csv'
MTVLabel='PlotVee/Hyper/VeeAtSector'//SecLabel//'.mtv'
WRITE(Out_Unit,'("VeeLabel=",A)')TRIM(VeeLabel)
WRITE(Out_Unit,'("CsvLabel=",A)')TRIM(CsvLabel)
!
! Check for existence of the contour namelist
!
INQUIRE (file =InputDIR(1:LEN(TRIM(InputDIR)))//InputFile, exist = there)
IF(there)THEN
   OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile)
   READ(In_Unit,NML=contour, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'Contour')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Contour3'
      WRITE(Msg_Unit,*)'ERROR with Namelist Contour4'
      STOP 'Graf_Hyper: Contour'
   ENDIF
   CLOSE(In_Unit)
   WRITE(Out_Unit,NML=contour)
ELSE
WRITE(Out_Unit,  * ) 'DATA must exist', InputDIR(1:LEN(TRIM(InputDIR)))//InputFile
   STOP 'Graf_Hyper'
ENDIF
!
! ALLOCATE memory
!
ALLOCATE (v (nx, ny) )
!
!  READ the potential surface name
!
WRITE(Out_Unit, '(a25)') PES_Name

OPEN(Unit=Vee_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//TRIM(VeeLabel))

!
!  fix rho and THEN plot the theta, chi surface as a projection.
!
WRITE(Out_Unit, '("rho=",f10.5)') rho
xi = - 1.d0
dx = 2.d0 / (nx - 1)
WRITE(Out_Unit, * ) 'nx,xi,dx=', nx, xi, dx
yi = - 1.d0
dy = 2.d0 / (ny - 1)
WRITE(Out_Unit, * ) 'ny,yi,dy=', ny, yi, dy

DO i = 1, nx
   xp (1) = xi + (i - 1) * dx
   DO j = 1, ny
      xp (2) = yi + (j - 1) * dy !  xp(1) and xp(2) are parallel to the x and y axes, resp., but cover the surface of the hyperspher
      rr = sqrt (xp (1) **2 + xp (2) **2)   !  aph coordinates
      theta = 2. * atan (rr)  !  stereographic projection
      sinth = sin (theta)
      chi (1) = acos (xp (1) / rr)   !  this version uses chi as its variable.
      chi (1) = sign (chi (1), xp (2) )
      chi (2) = chi (1) - chij (2, 1)
      chi (3) = chi (1) - chij (3, 1)
      DO k = 1, 3   !  get scaled and unscaled interparticle distances.
         s (k) = rho * sqrt (0.5 - 0.5 * sinth * cos (PhiChi * chi (k) ) )
         r (k) = d (k) * s (k)
      ENDDO
      CALL surface (vv, r)
      v (i, j) = autoev * vv   !  convert from au to ev
      IF(v (i, j) >vmax) v (i, j) = vmax
   ENDDO
ENDDO

IF(Debg)THEN
   WRITE(Out_Unit, '(" $readxl",/,1x,i5,2g13.5,/," $end")') nx, xi, dx
   WRITE(Out_Unit, '(" $readyl",/,1x,i5,2g13.5,/," $end",/,"0$DATA")') ny, yi, dy
   WRITE(Out_Unit, '(10g12.4)') ( (v (i, j) , j = 1, ny) , i = 1, nx)
   WRITE(Out_Unit, '(" $end")')
ENDIF

nchi = 361   !  compute a circle and WRITE it out
DO i = 1, nchi
   chii = (i - 1) * pi / 180.
   xx (i) = cos (chii)
   yy (i) = sin (chii)
ENDDO

WRITE(Vee_Unit, '(" $readxy",/,1x,i5)') nchi
WRITE(Vee_Unit, '(10g12.4)') (xx (i) , yy (i) , i = 1, nchi)
WRITE(Vee_Unit, '(" $contin")')
WRITE(Vee_Unit, '(" $readxl",/,1x,i5,2g13.5,/," $end")') nx, xi,dx
WRITE(Vee_Unit, '(" $readyl",/,1x,i5,2g13.5,/," $end",/,"0$DATA")') ny, yi, dy
WRITE(Vee_Unit, '(10g12.4)') ( (v (i, j) , j = 1, ny) , i = 1, nx)
WRITE(Vee_Unit, '(" $end")')

OPEN(Unit=CSV_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//TRIM(CsvLabel))
DO i = 1, nx
   xp (1) = xi + (i - 1) * dx
   DO j = 1, ny
      xp (2) = yi + (j - 1) * dy
      WRITE(CSV_Unit, '(3(e14.7,","))') xp (1) , xp (2) , v (i, j)
   ENDDO
ENDDO
CLOSE(Vee_Unit)
CLOSE(CSV_Unit)

WRITE(SubTitle,'(A12,A5,ES23.15,A1)') '%subtitle= "',' rho=', rho, '"'
Call MTV_Output(nx-1, -One, One, -One, One, v, 'hyp', MTVLabel, SubTitle)

DEALLOCATE (v)
ENDDO

CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/BasisForScatt.txt',Form='FORMATTED',ACCESS='APPEND')
WRITE(Out_Unit,*)'Completed Graf_Hyper'
WRITE(Out_Unit,*)

RETURN
ENDSUBROUTINE Graf_Hyper
