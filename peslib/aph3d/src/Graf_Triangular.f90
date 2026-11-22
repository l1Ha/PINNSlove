SUBROUTINE Graf_Triangular
USE Numeric_Kinds_Module
USE Numbers_Module
USE FileUnits_Module
USE PES_MODULE
USE Masses_Module
USE BasisRegions_Module
USE Boundary_Module
USE InputFile_Module
USE Contour_Module
USE Convrsns_Module
IMPLICIT NONE
LOGICAL :: maskit, debug, there, dtemp
PARAMETER (maskit = .false., debug = .false.)
INTEGER :: r2s, IERR
REAL(Kind=WP_Kind) :: vpot, dsq3
CHARACTER(LEN=4) SecLabel
CHARACTER(LEN=40) VeeLabel
CHARACTER(LEN=40) CsvLabel, MTVLabel
CHARACTER(LEN=80) SubTitle
REAL(Kind=WP_Kind) :: m0, m1, m2
REAL(Kind=WP_Kind) :: s3, s1, s2
INTEGER ::  nxy, ISector, len_trim, vtype
REAL(KIND=WP_Kind), ALLOCATABLE:: surf (:,:)
REAL(Kind=WP_Kind) :: delta1, delta2
REAL(Kind=WP_Kind) :: x10b, x02b, x10, x20, x02
REAL(Kind=WP_Kind) :: mu, r (3)
REAL(Kind=WP_Kind) :: d0, d1, d2
REAL(Kind=WP_Kind) :: x, y, th, perim, yc
INTEGER :: ii, jj
REAL(Kind=WP_Kind) :: xx, yy, deltac
DIMENSION xx (1000), yy (1000)
REAL(Kind=WP_Kind) :: xval, yval


WRITE(Out_Unit,*)'Called    Graf_Triangular'
CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/Graf_Triangular.txt',Form='FORMATTED')

INQUIRE (file =InputDIR(1:LEN(TRIM(InputDIR)))//InputFile, exist = there)
IF(there)THEN
   OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile)
   READ(In_Unit,NML=contour, IOSTAT=IERR)
   IF(IERR==-1)THEN
      CALL NameList_Default(Out_Unit, 'Contour')
   ELSEIF(IERR/=0)THEN
      WRITE(Out_Unit,*)'ERROR with Namelist Contour9'
      WRITE(Msg_Unit,*)'ERROR with Namelist Contour10'
      STOP 'Graf_Triangular: Contour'
   ENDIF
   CLOSE(In_Unit)
   WRITE(Out_Unit,NML=contour)
ELSE
WRITE(Out_Unit,  * ) 'DATA must exist', InputDIR(1:LEN(TRIM(InputDIR)))//InputFile
   STOP 'Graf_Triangular'
ENDIF


nxy=nx
ALLOCATE(Surf(nxy,nxy))

m0=mass(1)
m1=mass(2)
m2=mass(3)
!Mtot = m0 + m1 + m2
mu = SQRT (m0 * m1 * m2 / Mtot)
D0 = sqrt ( (M0 / mu) * (1.0d0 - M0 / Mtot) )
D1 = sqrt ( (M1 / mu) * (1.0d0 - M1 / Mtot) )
D2 = sqrt ( (M2 / mu) * (1.0d0 - M2 / Mtot) )
X10 = - acos ( - mu / (D0 * D1 * M2) )
X10b = asin ( - 1.0d0 / (D0 * D1) )
X02 = acos ( - mu / (D2 * D0 * M1) )
X02b = asin ( - 1.0d0 / (D2 * D0) )
IF(X02>X02b)THEN
   X02 = - 1.0d0 * abs (X02)
ELSE
   X02 = - 1.0d0 * abs (X02b)
ENDIF
X20 = - 1.0d0 * X02

WRITE(Out_Unit, * ) 'Subroutine Graf_Triangular'
WRITE(Out_Unit, * ) '---------------------------'
WRITE(Out_Unit, * ) 'PES= ', PES_Name
WRITE(Out_Unit, * ) 'Mass Constants:'
WRITE(Out_Unit, * ) 'Mass= ', M0, M1, M2
WRITE(Out_Unit, * ) 'D0= ', d0, ' D1= ', d1, ' D2= ', d2
WRITE(Out_Unit,  * ) 'X10= ', x10, '  X20= ', x20
WRITE(Out_Unit, * ) 'mu= ', mu, ' Mtot= ', Mtot
WRITE(Out_Unit, * ) 'nxy= ', nxy
WRITE(Out_Unit, * ) '============================'

surf=vmax
dsq3 = SQRT (3d0)
DO ISector = 1, NSectors
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
      STOP 'Graf_Triangular'
   ENDIF
   VeeLabel='PlotVee/Triangular/VeeTriangular'//SecLabel//'.wg'
   CsvLabel='PlotVee/Triangular/VeeTriangular'//SecLabel//'.csv'
   MTVLabel='PlotVee/Triangular/VeeTriangular'//SecLabel//'.mtv'
   WRITE(Out_Unit,'("VeeLabel=",A30)')VeeLabel
   WRITE(Out_Unit,'("CsvLabel=",A31)')CsvLabel
!   perim = 2*Basis_Dist(ISector)
   perim= Basis_Dist(ISector)*(6.d0/Sqrt(7.d0)+2.d0)/2.d0
   xval = perim / (2d0 * SQRT (2d0) )
   yval = perim / (SQRT (6d0) )
   s1 = perim / SQRT (3d0)
   delta1 = 2.0d0 * xval / (nxy - 1)
   delta2 = 3.0d0 * yval / 2.d0 / (nxy - 1)
   DO ii = 1, nxy
      x = - xval + DBLE (ii - 1) * delta1
      s2 = x
      DO jj = 1, nxy
         y = - yval + delta2 * DBLE (jj - 1)
         IF(x<0d0)THEN
            yc = y / 2d0 + x * dsq3
         ELSE
            yc = y / 2d0 - x * dsq3
         ENDIF
         s3 = y
         r (1) = (s1 + SQRT (2.0d0) * s3) / dsq3  *D0
         r (2) = (s1 + SQRT (1.5d0) * s2 - SQRT (0.5d0) * s3) / dsq3  *D1
         r (3) = (s1 - SQRT (1.5d0) * s2 - SQRT (0.5d0) * s3) / dsq3  *D2
         IF(r(1)<0.7d0 .or. r(2)<0.7d0 .or. r(3)<0.7d0)THEN
            cycle
         ENDIF
         CALL surface (vpot, r)
         vpot=autoev * vpot !  convert from au to ev
         IF(ISector==1.AND.vpot<0.d0)THEN
            dtemp=.true.
         ELSE
            dtemp=.false.
         ENDIF
         IF(dtemp)THEN
            WRITE(Out_Unit, '(a9,3f16.10)') 'x,y=     ', x, y
            WRITE(Out_Unit, '(a9,3f16.10)') 's2,s3,s1=', s2, s3, s1
            WRITE(Out_Unit, '(a9,3f16.10)') 'r1,r2,r3=', r (1) , r (2) , r (3)
            WRITE(Out_Unit, '(a9,3f16.10)') 'vpot=    ', vpot
            WRITE(Out_Unit, '(a9,3f16.10)') 'p,p-sumri', perim, perim - r (1) - r (2) - r (3)
            WRITE(Out_Unit, * ) '----------------------------------'
         ENDIF
         IF(vpot>vmax) vpot = vmax
         surf (ii, jj) = vpot
      ENDDO
   ENDDO

   WRITE(Out_Unit,'("VeeLabel=",A)')TRIM(VeeLabel)
   WRITE(Out_Unit,'("CsvLabel=",A)')TRIM(CsvLabel)
   WRITE(Out_Unit,'("MTVLabel=",A)')TRIM(MTVLabel)

   OPEN(Unit=CSV_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//TRIM(CsvLabel))
   DO ii = 1, nx
      x = - xval + DBLE (ii - 1) * delta1
      DO jj = 1, ny
         y = - yval + delta2 * DBLE (jj - 1)
         WRITE(CSV_Unit, '(3(e14.7,","))')x, y, surf(ii, jj)
      ENDDO
   ENDDO
   CLOSE(CSV_Unit)

   OPEN(Unit=Vee_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//TRIM(VeeLabel))
   WRITE(Vee_Unit, '(" $readxl",/,1x,i5,2g13.5,/," $end")') nx, -xval, delta1
   WRITE(Vee_Unit, '(" $readyl",/,1x,i5,2g13.5,/," $end",/," $data")') ny, -yval, delta2
   WRITE(Vee_Unit, '(10g12.4)') ( (surf(ii, jj) , jj = 1, ny) , ii = 1,nx)
   WRITE(Vee_Unit, '(" $end")')
   CLOSE(Vee_Unit)

   WRITE(SubTitle,'(A12,A11,ES23.15,A8,F8.5,A8,F8.5,A1)') '%subtitle= "',&
   ' Perimeter=', Perim,' Rhomin=',Sqrt(7.d0)*Perim/6,' Rhomax=',Perim/2, '"'
   CALL MTV_Output(nxy-1, -xval, xval, -yval, -yval+(nxy - 1)*delta2, surf, &
   'spl', MTVLabel, SubTitle)
ENDDO
DEALLOCATE(surf)

CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/BasisForScatt.txt',Form='FORMATTED',ACCESS='APPEND')
WRITE(Out_Unit,*)'Completed Graf_Triangular'
WRITE(Out_Unit,*)

RETURN
ENDSUBROUTINE Graf_Triangular
