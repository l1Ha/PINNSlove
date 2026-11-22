SUBROUTINE Graf_Surm
USE Numbers_Module
USE FileUnits_Module
USE Masses_Module
USE MassFactor2_Module
USE PES_MODULE
USE InputFile_Module
USE BasisRegions_Module
USE Boundary_Module
USE Contour_Module  !, ONLY: nx, ny, vmax, rho, iaph
USE Convrsns_Module
IMPLICIT NONE
!
!  calculates potential over a grid of points to be used in making
!  contour plots.
!  looks for minimum reaction path.
!  this version calls surface SUBROUTINE in same FORMAT as aph3d, and
!  needs it loaded as a binary.
LOGICAL :: there
LOGICAL, PARAMETER:: Debg=.False.
INTEGER i, irho, j, k, nchi, IERR
CHARACTER(LEN=28) MTVLabel
CHARACTER(LEN=80) SubTitle
REAL(Kind=WP_Kind) chi(3)
REAL(Kind=WP_Kind) xp (3), xx (361), yy (361)
REAL(Kind=WP_Kind) rlittle(3)
REAL(Kind=WP_Kind) chii, dx, dy, rr, sinth, theta, vv, Rmax, xi, yi
REAL(Kind=WP_Kind), PARAMETER:: PhiChi=1.d0
REAL(KIND=WP_Kind), ALLOCATABLE:: v(:,:), Rhom(:,:)


WRITE(Out_Unit,*)'Called    Graf_Surm'
CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/Graf_Surm.txt',Form='FORMATTED')

INQUIRE (file =InputDIR(1:LEN(TRIM(InputDIR)))//InputFile, exist = there)
IF(there)THEN
   OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile)
   READ(In_Unit,NML=contour, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'Contour')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Contour7'
      WRITE(Msg_Unit,*)'ERROR with Namelist Contour8'
      STOP 'Graf_Surm: Contour'
   ENDIF
   CLOSE(In_Unit)
   WRITE(Out_Unit,NML=contour)
ELSE
WRITE(Out_Unit,  * ) 'DATA must exist', InputDIR(1:LEN(TRIM(InputDIR)))//InputFile
   STOP 'Graf_Surm'
ENDIF
ALLOCATE(V(nx,ny),Rhom(nx,ny))

rmax = 1.d0

xi = - 1.d0
dx = 2.d0 / (nx - 1)
yi = - 1.d0
dy = 2.d0 / (ny - 1)
v=1.d35
DO i = 1, nx
   xp (1) = xi + (i - 1) * dx
   DO j = 1, ny
      xp (2) = yi + (j - 1) * dy  !  xp(1) and xp(2) are parallel to the x and y axes, resp., but coverthe surface of the hyperspher
      rr = sqrt (xp (1) **2 + xp (2) **2)  !  aph coordinates
      rhom (i, j) = 10.0
      v(i,j)=vmax
      theta = 2. * atan (rr)   !  stereographic projection
      sinth = sin (theta)
      chi (1) = acos (xp (1) / rr)   !  this version uses chi as its variable.
      chi (1) = sign (chi (1), xp (2) )
      chi (2) = chi (1) - chij (2, 1)
      chi (3) = chi (1) - chij (3, 1)
      DO irho = 1, NSectors    !  start loop on rho to look for min at the given theta and chi
         rho = Basis_Dist(irho)
         DO k = 1, 3           !  get scaled and unscaled interparticle distances.
            s (k) = rho * sqrt (0.5 - 0.5 * sinth * cos (PhiChi * chi(k) ) )
            rlittle (k) = d (k) * s (k)
         ENDDO
         CALL surface (vv, rlittle)  !  ***calculate potential energy***
         vv = autoev * vv  !  convert from au to ev
         IF(vv<v (i, j) )THEN
            v (i, j) = vv
            rhom (i, j) = rho
         ENDIF
      ENDDO
   ENDDO
ENDDO

IF(debg)THEN
   WRITE(Out_Unit, 64) nx, xi, dx
   WRITE(Out_Unit, 66) ny, yi, dy
   WRITE(Out_Unit, '(1x,10g12.4)') ( (v (i, j), j = 1, ny), i = 1, nx)
ENDIF

nchi = 361
DO i = 1, nchi   !  compute a circle and WRITE it out
   chii = (i - 1) * pi / 180.
   xx (i) = cos (chii)
   yy (i) = sin (chii)
ENDDO

OPEN(Unit=Vee_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'PlotVee/Surm/VeeSurm.wg')
WRITE(Vee_Unit, 68) nchi
WRITE(Vee_Unit, '(1x,10g12.4)') (xx (i), yy (i), i = 1, nchi)
WRITE(Vee_Unit, '( " $CONTINUE")')
WRITE(Vee_Unit, 64) nx, xi, dx
WRITE(Vee_Unit, 66) ny, yi, dy
WRITE(Vee_Unit, '(1x,10g12.4)') ( (v (i, j), j = 1, ny), i = 1, nx)
WRITE(Vee_Unit, '( " $END" )')
CLOSE(Vee_Unit)

OPEN(Unit=Vee_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'PlotVee/Surm/RhoSurm.wg')
WRITE(Vee_Unit, 68) nchi
WRITE(Vee_Unit, '(1x,10g12.4)') (xx (i), yy (i), i = 1, nchi)
WRITE(Vee_Unit, '( " $CONTINUE")')
WRITE(Vee_Unit, 64) nx, xi, dx
WRITE(Vee_Unit, 66) ny, yi, dy
WRITE(Vee_Unit, '(1x,10g12.4)') ( (rhom (i, j), j = 1, ny), i = 1, nx)
WRITE(Vee_Unit, '( " $END" )')
CLOSE(Vee_Unit)

OPEN(Unit=CSV_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'PlotVee/Surm/VeeSurm.csv')
DO i = 1, nx
   xp (1) = xi + (i - 1) * dx
   DO j = 1, ny
      xp (2) = yi + (j - 1) * dy
      WRITE(CSV_Unit, '(3(e14.7,","))') xp (1) , xp (2) , v (i, j)
   ENDDO
ENDDO
CLOSE(CSV_Unit)

OPEN(Unit=CSV_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'PlotVee/Surm/RhoSurm.csv')
DO i = 1, nx
   xp (1) = xi + (i - 1) * dx
   DO j = 1, ny
      xp (2) = yi + (j - 1) * dy
      WRITE(CSV_Unit, '(3(e14.7,","))') xp (1) , xp (2) , Rhom(i, j)
   ENDDO
ENDDO
CLOSE(CSV_Unit)

MTVLabel='PlotVee/Surm/VeeSurm.mtv'
WRITE(SubTitle,'(A12,A8,ES23.15,A8,ES23.15,A1)') '%subtitle= "',' rhomin=', Basis_Dist(1),' rhomax=',Basis_Dist(NSectors+1), '"'
Call MTV_Output(nx-1, xi, -xi, yi, -yi, v,    'hyp', MTVLabel, SubTitle)
MTVLabel='PlotVee/Surm/RhoSurm.mtv'
Call MTV_Output(nx-1, xi, -xi, yi, -yi, Rhom, 'hyp', MTVLabel, SubTitle)

DEALLOCATE(V,Rhom)

CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/BasisForScatt.txt',Form='FORMATTED',ACCESS='APPEND')
WRITE(Out_Unit,*)'Completed Graf_Surm'
WRITE(Out_Unit,*)

RETURN
   52 FORMAT(i5,2f10.4)
   61 FORMAT(/,'m(i),mu,mm=',5ES15.7)
   62 FORMAT(/,'d(i)=',3f10.5)
   63 FORMAT(1x,3f10.5)
   64 FORMAT( ' $readxl',/, 1x,i5,2g13.5, /, ' $end' )
   66 FORMAT( ' $readyl',/, 1x,i5,2g13.5, /, ' $end', /, '0$DATA' )
   68 FORMAT( ' $readxy',/,1x,i5)
  177 FORMAT(1000(e14.7,','))
  631 FORMAT(/,'   chij')
  632 FORMAT(/,'nrho,rhomin,deltarho1,deltarho2=',i5,3f10.4)
ENDSUBROUTINE Graf_surm
