SUBROUTINE Graf_Planar(RhoMax)
USE Numbers_Module
USE FileUnits_Module
USE DVRMasses_Module
USE PES_MODULE
USE InputFile_Module
USE Masses_Module
USE MassFactor2_Module
USE Contour_Module
USE Convrsns_Module
!
!  program Graf_Planar
!
IMPLICIT NONE
LOGICAL :: there
LOGICAL, PARAMETER:: debg=.False.
!
!  calculates potential over a grid of points to be used in making contour plots.
!  this version calls surface SUBROUTINE in same FORMAT as aph3d, and
!  needs it loaded as a binary.
!
CHARACTER(LEN=1)  ChiLabel
CHARACTER(LEN=34) VeeLabel
CHARACTER(LEN=35) CsvLabel, MTVLabel
CHARACTER(LEN=80) SubTitle
INTEGER :: i, j, k, ichi, IERR
REAL(Kind=WP_Kind) :: chi(3)
REAL(Kind=WP_Kind) :: rhomax, xi, dx, yi, dy, rr, theta, sinth, vv, chival
REAL(Kind=WP_Kind) :: xp (3)
REAL(Kind=WP_Kind), PARAMETER:: PhiChi=1.d0
REAL(KIND=WP_Kind), ALLOCATABLE::v(:, :)


WRITE(Out_Unit,*)'Called    Graf_Planar'
CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/Graf_Planar.txt',Form='FORMATTED')
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
      WRITE(Out_Unit,*)'ERROR with Namelist Contour5'
      WRITE(Msg_Unit,*)'ERROR with Namelist Contour6'
      STOP 'Graf_Planar: Contour'
   ENDIF
   CLOSE(In_Unit)
   WRITE(Out_Unit,NML=contour)
ELSE
WRITE(Out_Unit,  * ) 'DATA file must exist: ',InputDIR(1:LEN(TRIM(InputDIR)))//InputFile
   STOP 'Graf_Planar'
ENDIF
!
! ALLOCATE memory
!
ALLOCATE(v(nx, ny))

OPEN(Unit=Vee_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'PlotVee/Planar/Vee_Planar.wg')

!
!  fix rho and THEN plot the theta, chi surface as a projection.
!
WRITE(Out_Unit, '("rhomax=",f10.5)') rhomax
xi = 0.d0
dx = 1.d0 / (nx - 1)
WRITE(Out_Unit, * ) 'nx,xi,dx=', nx, xi, dx
yi = 0.d0
dy = 1.d0 / (ny - 1)
WRITE(Out_Unit, * ) 'ny,yi,dy=', ny, yi, dy
DO ichi=1,6
   IF(ichi==1)THEN
      chival=0.d0
   ELSEIF(ichi==2)THEN
      chival=(pi+chij(1,3))/2
   ELSEIF(ichi==3)THEN
      chival=pi+chij(1,3)
   ELSEIF(ichi==4)THEN
      chival=(pi+chij(1,3)+chij(1,2))/2
   ELSEIF(ichi==5)THEN
      chival=chij(1,2)
   ELSEIF(ichi==6)THEN
      chival=(pi+chij(1,2))/2
   ELSE
      chival=pi
   ENDIF
   WRITE(Out_Unit,*)'Plotting Plane at constant Chi=',chival
   WRITE(ChiLabel,'(I1)')ichi
   VeeLabel='PlotVee/ChiPlane/VeeAtChiPlane'//ChiLabel//'.wg'
   CsvLabel='PlotVee/ChiPlane/VeeAtChiPlane'//ChiLabel//'.csv'
   MTVLabel='PlotVee/ChiPlane/VeeAtChiPlane'//ChiLabel//'.mtv'
   DO i = 1, nx
      xp (1) = xi + (i - 1) * dx
      xp (1) = xp (1) * rhomax
      DO j = 1, ny
         xp (2) = yi + (j - 1) * dy !xp(1) and xp(2) are parallel to the x and y axes, resp., but cover the surface of the hypersphe
         xp (2) = xp (2) * rhomax
         rr = sqrt (xp (1) **2 + xp (2) **2) !  aph coordinates
         theta = atan2(xp(2),xp(1))
         IF(i==1.or.i==nx)THEN
            IF(j==1.or.j==ny)THEN
               WRITE(Out_Unit,'("x=",f10.5," y=",f10.5," rho=",f10.5," theta=", f10.5," chi=",2f10.5)')&
               xp(1),xp(2),rr,theta,chival,chival*180/pi
            ENDIF
         ENDIF
         sinth = sin (theta)
         chi (1) = chival !  this version uses chi as its variable.
         chi (2) = chi (1) - chij (2, 1)
         chi (3) = chi (1) - chij (3, 1)
         DO k = 1, 3   !  get scaled and unscaled interparticle distances.
            s (k) = rr * sqrt (0.5d0 - 0.5d0 * sinth * cos (PhiChi * chi(k) ) )
            r (k) = d (k) * s (k)
         ENDDO
         IF(i==1.or.i==nx)THEN
            IF(j==1.or.j==ny)THEN
               WRITE(Out_Unit,*)s
               WRITE(Out_Unit,*)r
            ENDIF
         ENDIF
         CALL surface (vv, r)
         v(i, j) = autoev * vv !  convert from au to ev
         IF(v(i, j) >vmax) v(i, j) = vmax
      ENDDO
   ENDDO

   OPEN(Unit=Vee_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//TRIM(VeeLabel))
   WRITE(Vee_Unit, '(" $readxl",/,1x,i5,2g13.5,/," $end")') nx, xi, dx
   WRITE(Vee_Unit, '(" $readyl",/,1x,i5,2g13.5,/," $end",/,"0$DATA")') ny, yi, dy
   WRITE(Vee_Unit, '(10g12.4)') ( (v(i, j) , j = 1, ny) , i = 1,nx)
   WRITE(Vee_Unit, '(" $end")')
   CLOSE(Vee_Unit)

   OPEN(Unit=CSV_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//TRIM(CsvLabel))
   DO i = 1, nx
      xp (1) = xi + (i - 1) * dx
      xp (1) = xp (1) * rhomax
      DO j = 1, ny
         xp (2) = yi + (j - 1) * dy
         xp (2) = xp (2) * rhomax
         WRITE(CSV_Unit, '(3(e14.7,","))') xp (1) , xp (2) , v(i, j)
      ENDDO
   ENDDO
   CLOSE(CSV_Unit)

   WRITE(SubTitle,'(A12,A5,ES23.15,A6,F6.2,A8)') '%subtitle= "',' Chi=', chival, ' (Rad)', chival*180/Pi, ' (Deg)"'
   WRITE(Out_Unit,'(A12,A5,ES23.15,A1)') '%subtitle= "',' Chi=', chival, '"'
   Call MTV_Output(nx-1, Zero, rhomax, Zero, rhomax, v, 'pla', MTVLabel, SubTitle)
ENDDO

DEALLOCATE(v)

CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/BasisForScatt.txt',Form='FORMATTED',ACCESS='APPEND')
WRITE(Out_Unit,*)'Completed Graf_Planar'
WRITE(Out_Unit,*)

RETURN
ENDSUBROUTINE Graf_Planar
