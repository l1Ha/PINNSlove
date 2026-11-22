SUBROUTINE Graf_Collinear(RhoMax)
USE Numbers_Module
USE FileUnits_Module
USE Masses_Module
USE MassFactor2_Module
USE PES_MODULE
USE InputFile_Module
USE Contour_Module
USE Convrsns_Module
!
!  program Graf_Collinear
!
IMPLICIT NONE
LOGICAL :: there
LOGICAL, PARAMETER:: debg=.False.
!
!  calculates potential over a grid of points to be used in making contour plots.
!  this version calls surface SUBROUTINE in same FORMAT as aph3d, and
!  needs it loaded as a binary.
!
CHARACTER(LEN=34) MTVLabel
CHARACTER(LEN=80) SubTitle
INTEGER :: i, j, k, nchi, IERR
REAL(Kind=WP_Kind) :: chi(3)
REAL(Kind=WP_Kind) :: rhomax, xi, dx, yi, dy, rr, theta, sinth, vv
REAL(Kind=WP_Kind) :: chii, xp (3)
REAL(Kind=WP_Kind), PARAMETER:: PhiChi=1.d0
REAL(KIND=WP_Kind), ALLOCATABLE::v (:, :)
REAL :: xx (361), yy (361)



WRITE(Out_Unit,*)'Called    Graf_Collinear'
CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/Graf_Collinear.txt',Form='FORMATTED')

vmax = Ten**(RANGE(One)-10)
!
! Check for existence of file .contour
!
INQUIRE (file = InputDIR(1:LEN(TRIM(InputDIR)))//InputFile, exist = there)
IF(there)THEN
   OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile)
   READ(In_Unit,NML=Contour, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'Contour')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Contour1'
      WRITE(Msg_Unit,*)'ERROR with Namelist Contour2'
      STOP 'Graf_Collinear: Contour'
   ENDIF
   CLOSE(In_Unit)
   WRITE(Out_Unit,NML=contour)
ELSE
WRITE(Out_Unit,  * ) 'DATA file must exist: ',InputDIR(1:LEN(TRIM(InputDIR)))//InputFile
   STOP 'Graf_Collinear'
ENDIF
!
! ALLOCATE memory
!
ALLOCATE(v(nx, ny))

OPEN(Unit=Vee_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'PlotVee/Collinear/Vee_Collinear.wg')

!
!  fix rho and THEN plot the theta, chi surface as a projection.
!
WRITE(Out_Unit, '("rhomax=",f10.5)') rhomax
xi = - 1.d0
dx = 2.d0 / (nx - 1)
WRITE(Out_Unit, * ) 'nx,xi,dx=', nx, xi, dx
yi = - 1.d0
dy = 2.d0 / (ny - 1)
WRITE(Out_Unit, * ) 'ny,yi,dy=', ny, yi, dy
DO i = 1, nx
   xp (1) = xi + (i - 1) * dx
   xp (1) = xp (1) * rhomax
   DO j = 1, ny
      xp (2) = yi + (j - 1) * dy !xp(1) and xp(2) are parallel to the x and y axes, resp., but cover the surface of the hypersphere.
      xp (2) = xp (2) * rhomax
      rr = sqrt (xp (1) **2 + xp (2) **2) !  aph coordinates
      theta = pi / 2.d0  !  stereographic projection
      sinth = sin (theta)
      chi (1) = atan2 (xp (2), xp (1) ) !  this version uses chi as its variable.
      chi (2) = chi (1) - chij (2, 1)
      chi (3) = chi (1) - chij (3, 1)
      DO k = 1, 3   !  get scaled and unscaled interparticle distances.
         s (k) = rr * sqrt (0.5d0 - 0.5d0 * sinth * cos (PhiChi * chi (k) ) )
         r (k) = d (k) * s (k)
      ENDDO
      CALL surface (vv, r)
      v (i, j) = autoev * vv !  convert from au to ev
      IF(v (i, j) >vmax) v (i, j) = vmax
   ENDDO
ENDDO
IF(Debg)THEN
   WRITE(Out_Unit, '(" $readxl",/,1x,i5,2g13.5,/," $end")') nx, xi, dx
   WRITE(Out_Unit, '(" $readyl",/,1x,i5,2g13.5,/," $end",/,"0$DATA")') ny, yi, dy
   WRITE(Out_Unit, '(10g12.4)') ( (v (i, j) , j = 1, ny) , i = 1,nx)
   WRITE(Out_Unit, '(" $end")')
ENDIF
!
!  compute a circle and WRITE it out
!
nchi = 361
DO i = 1, nchi
   chii = (i - 1) * pi / 180.
   xx (i) = cos (chii)
   yy (i) = sin (chii)
ENDDO

WRITE(Vee_Unit, '(" $readxy",/,1x,i5)') nchi
WRITE(Vee_Unit, '(10g12.4)') (xx (i) , yy (i) , i = 1, nchi)
WRITE(Vee_Unit, '(" $CONTINUE")')
WRITE(Vee_Unit, '(" $readxl",/,1x,i5,2g13.5,/,"0$END")') nx, xi,dx
WRITE(Vee_Unit, '(" $readyl",/,1x,i5,2g13.5,/,"0$END",/," $data")') ny, yi, dy
WRITE(Vee_Unit, '(10g12.4)') ( (v (i, j) , j = 1, ny) , i = 1, nx)
WRITE(Vee_Unit, '(" $end")')

OPEN(Unit=CSV_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'PlotVee/Collinear/Vee_Collinear.csv')
DO i = 1, nx
   xp (1) = xi + (i - 1) * dx
   xp (1) = xp (1) * rhomax
   DO j = 1, ny
      xp (2) = yi + (j - 1) * dy
      xp (2) = xp (2) * rhomax
      WRITE(CSV_Unit, '(3(e14.7,","))') xp (1) , xp (2) , v (i, j)
   ENDDO
ENDDO

CLOSE(Vee_Unit)
CLOSE(CSV_Unit)

MTVLabel="PlotVee/Collinear/VeeCollinear.mtv"
WRITE(SubTitle,'(A26)') '%subtitle= Collinear Plane"'
Call MTV_Output(nx-1, -rhomax, rhomax, -rhomax, rhomax, v, 'pla', MTVLabel, SubTitle)

DEALLOCATE(v)

CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/BasisForScatt.txt',Form='FORMATTED',ACCESS='APPEND')
WRITE(Out_Unit,*)'Completed Graf_Collinear'
WRITE(Out_Unit,*)

RETURN
ENDSUBROUTINE Graf_Collinear
