SUBROUTINE mtv_output (nxy, xmin, xmax, ymin, ymax, surf, pltype, MTVLabel, SubTitle)
USE FileUnits_Module
USE Numeric_Kinds_Module
USE Numbers_Module
USE PES_MODULE
IMPLICIT NONE
INTEGER :: nxy1, nxy2, nxy3, nxy4
LOGICAL :: maskit
CHARACTER(LEN=3) :: pltype
CHARACTER(LEN=*) SubTitle
CHARACTER(LEN=*) MTVLabel
INTEGER :: nxy
REAL(KIND=WP_Kind) :: surf (nxy + 1, nxy + 1)
REAL(Kind=WP_Kind) :: vmax, vmin, vpot
REAL(Kind=WP_Kind) :: x, y, th, xi, deltac
INTEGER :: ii, jj
INTEGER :: len_trim
REAL(Kind=WP_Kind) :: xx (0:nxy), yy (0:nxy)
REAL(Kind=WP_Kind) :: cmin, cmax, cminsave, cmaxsave
REAL(Kind=WP_Kind) :: xmin, xmax, ymin, ymax

WRITE(Out_Unit,*)
WRITE(Out_Unit,*)"Start of MTV_Output"
cmin=MinVal(surf)
cmax=MaxVal(surf)
cminsave = cmin
cmaxsave = cmax

WRITE(Out_Unit, * ) 'nxy= ', nxy
WRITE(Out_Unit, * ) 'cmin= ', cmin, ' cmax= ', cmax

vmin=minval(surf)
vmax=maxval(surf)
IF(pltype=='hyp')THEN
   maskit = .true.
   deltac = 2.0d0 * PI / DBLE (nxy)
   th = PI / 2.0d0
   WRITE(Out_Unit, * ) 'Creating Circle DATA'
   DO ii = 0, nxy
   xi = DBLE (ii) * deltac
   xx (ii) = tan (th / 2.0d0) * cos (xi)
   yy (ii) = tan (th / 2.0d0) * sin (xi)
   ENDDO
ELSEIF (pltype=='spl')THEN
   maskit = .true.
   xx (1) = xmin
   yy (1) = ymax
   xx (2) = 0d0
   yy (2) = ymin
   xx (3) = xmax
   yy (3) = ymax
   xx (4) = xmin
   yy (4) = ymax
ELSEIF (pltype=='pla')THEN
   maskit = .false.
   xx (1) = xmin
   yy (1) = ymax
   xx (2) = 0d0
   yy (2) = ymin
   xx (3) = xmax
   yy (3) = ymax
   xx (4) = xmin
   yy (4) = ymax
ELSE
   STOP 'wrong type'
ENDIF

IF(vmax- vmin<cmaxsave-cminsave)THEN
   cmin = vmin
   cmax = vmax
ELSEIF (vmax<cmaxsave)THEN
   cmin = vmax- (cmaxsave-cminsave)
   IF(cmin<vmin) cmin = vmin
   cmax = vmax
ELSEIF (vmin>cminsave)THEN
   cmax = vmin + (cmaxsave-cminsave)
   cmin = vmin
ELSE
   WRITE(Out_Unit, * ) 'cmin=', cmin
   WRITE(Out_Unit, * ) 'vmin=', vmin
   WRITE(Out_Unit, * ) 'cmax=', cmax
   WRITE(Out_Unit, * ) 'vmax=', vmax
!             STOP 'error'
ENDIF

IF(cmin>cmax)THEN
   STOP 'cmin>cmax'
ENDIF

OPEN(unit =MTV_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//MTVLabel, status = 'unknown')

WRITE(MTV_Unit, * ) '$DATA=CONTOUR'
IF(pltype=='hyp')THEN
   WRITE(MTV_Unit, '(A12,A25,A1)' ) '%toplabel= "', PES_Name,'"'
ELSE
   WRITE(MTV_Unit, '(A12,A25,A1)' ) '%toplabel= "', PES_Name,'"'
ENDIF

WRITE(MTV_Unit, '(A80)')SubTitle
WRITE(MTV_Unit, '(a5,i5,a6,e16.4,a6,e16.4)') '%nx= ', nxy, ' xmin=',xmin, ' xmax=', xmax
WRITE(MTV_Unit, '(a5,i5,a6,e16.4,a6,e16.4)') '%ny= ', nxy, ' ymin=',ymin, ' ymax=', ymax
WRITE(MTV_Unit, * ) '%xlabel= "x"'
WRITE(MTV_Unit, * ) '%ylabel= "y"'
WRITE(MTV_Unit, * ) '%xticks=10 yticks=10'
WRITE(MTV_Unit, * ) '%nsteps=50'
WRITE(MTV_Unit, * ) '%contfill=true'
WRITE(MTV_Unit, * ) '%interpolate=1'
WRITE(MTV_Unit, * ) '%cmin=', cmin
WRITE(MTV_Unit, * ) '%cmax=', cmax

DO jj = 1, nxy
   DO ii = 1, nxy
      WRITE(MTV_Unit, * ) surf (ii, jj)
   ENDDO
   WRITE(MTV_Unit, * )
ENDDO

IF(pltype=='hyp')THEN   !  HYP type
   DO ii = 0, nxy - 1
      WRITE(MTV_Unit, 500)'@line lc=4 lw=3 x1=',xx(ii),' y1=',yy(ii),' z1=',Zero,' x2=',xx(ii+1),' y2=',yy(ii+1),' z2=',Zero
   ENDDO
   WRITE(MTV_Unit, * )
   IF(maskit)THEN
      !FIRST MASK
      WRITE(MTV_Unit, * ) '$DATA=CURVE2D'
      WRITE(MTV_Unit,  * ) '%fillcolor=-1  filltype=1 linetype=0'
      nxy1 = nxy / 4
      nxy2 = nxy / 2
      nxy3 = 3 * nxy / 4
      nxy4 = nxy
      DO ii = 0, nxy1
      WRITE(MTV_Unit, '(2f14.9)') xx (ii) , yy (ii)
      ENDDO
      IF(abs (xx (nxy1) ) >1.d-8)THEN
         WRITE(MTV_Unit, '(f14.9,a5)') xx (nxy1) , ' 1.0'
      ENDIF
      WRITE(MTV_Unit, * ) ' 1.0 1.0'
      IF(abs (yy (0) - 1.d0) >1.d-8)THEN
         WRITE(MTV_Unit, '(a5,f14.9)') ' 1.0', yy (0)
      ENDIF
      WRITE(MTV_Unit, '(2f14.9)') xx (0) , yy (0)
      WRITE(MTV_Unit, * )
      !SECOND MASK
      WRITE(MTV_Unit, * ) '$DATA=CURVE2D'
      WRITE(MTV_Unit, * ) '%fillcolor=-1 filltype=1 linetype=0'
      DO ii = nxy1, nxy2
      WRITE(MTV_Unit, '(2f14.9)') xx (ii) , yy (ii)
      ENDDO
      IF(abs (yy (nxy2) ) >1.d-8)THEN
         WRITE(MTV_Unit, '(a5,f14.9)') '-1.0', yy (nxy2)
      ENDIF
      WRITE(MTV_Unit, * ) '-1.0 1.0'
      IF(abs (xx (nxy1) - 1.d0) >1.d-8)THEN
         WRITE(MTV_Unit, '(f14.9,a5)') xx (nxy1) , ' 1.0'
      ENDIF
      WRITE(MTV_Unit, '(2f14.9)') xx (nxy1) , yy (nxy1)
      WRITE(MTV_Unit, * )
      !THIRD MASK
      WRITE(MTV_Unit, * ) '$DATA=CURVE2D'
      WRITE(MTV_Unit, * ) '%fillcolor=-1 filltype=1 linetype=0'
      DO ii = nxy2, nxy3
      WRITE(MTV_Unit, '(2f14.9)') xx (ii) , yy (ii)
      ENDDO
      IF(abs (xx (nxy3) ) >1.d-8)THEN
         WRITE(MTV_Unit, '(f14.9,a5)') xx (nxy3) , ' -1.0'
      ENDIF
      WRITE(MTV_Unit, * ) '-1.0 -1.0'
      IF(abs (yy (nxy2) ) >1.d-8)THEN
         WRITE(MTV_Unit, '(a5,f14.9)') '-1.0', yy (nxy2)
      ENDIF
      WRITE(MTV_Unit, '(2f14.9)') xx (nxy2) , yy (nxy2)
      WRITE(MTV_Unit, * )
      !FOURTH MASK
      WRITE(MTV_Unit, * ) '$DATA=CURVE2D'
      WRITE(MTV_Unit, * ) '%fillcolor=-1 filltype=1 linetype=0'
      DO ii = nxy3, nxy4
      WRITE(MTV_Unit, '(2f14.9)') xx (ii) , yy (ii)
      ENDDO
      IF(abs (yy (nxy4) + 1.d0) >1.d-8)THEN
         WRITE(MTV_Unit, '(a5,f14.9)') ' 1.0', yy (nxy4)
      ENDIF
      WRITE(MTV_Unit, * ) ' 1.0 -1.0'
      IF(abs (xx (nxy3) ) >1.e-8)THEN
         WRITE(MTV_Unit, '(f14.9,a5)') xx (nxy3) , ' -1.0'
      ENDIF
      WRITE(MTV_Unit, '(2f14.9)') xx (nxy3) , yy (nxy3)
      WRITE(MTV_Unit, * )
   ENDIF
ELSEIF(PlType/='pla')THEN      !  TRian type
   DO ii = 1, 3
      !WRITE( *, * ) xx (ii), yy (ii)
      WRITE(MTV_Unit,500) '@line lc=4 lw=4 x1=',xx(ii),' y1=',yy(ii),' z1=', Zero,' x2=',xx(ii+1),' y2=',yy(ii+ 1),' z2=',Zero
   ENDDO
   WRITE(MTV_Unit,  * ) '@point x1=0.0 y1=0.0 z1=0.0 ', 'mt=13 mc=5 ft=1'
   IF(maskit)THEN
      WRITE(MTV_Unit, * ) '$DATA=CURVE2D'
      WRITE(MTV_Unit, * ) '%fillcolor=-1 filltype=1'
      WRITE(MTV_Unit, 600) xx(1), yy(1)
      WRITE(MTV_Unit, 600) xx(2), yy(2)
      WRITE(MTV_Unit, 600) xx(1), yy(2)
      WRITE(MTV_Unit, 600) xx(1), yy(1)
      WRITE(MTV_Unit, * )
      WRITE(MTV_Unit, * ) '%fillcolor=-1 filltype=1'
      WRITE(MTV_Unit, 600) xx(2), yy(2)
      WRITE(MTV_Unit, 600) xx(3), yy(3)
      WRITE(MTV_Unit, 600) xx(3), yy(2)
      WRITE(MTV_Unit, 600) xx(2), yy(2)
      WRITE(MTV_Unit, * )
   ENDIF
ENDIF
WRITE(MTV_Unit, * ) '$END'
CLOSE(MTV_Unit)

WRITE(Out_Unit,*)"End of MTV_Output"
RETURN
  500 FORMAT(A19,f9.4,A4,f9.4,A4,f9.4,A4,f9.4,A4,f9.4,A4,f9.4)
  600 FORMAT(2f9.4)
ENDSUBROUTINE mtv_output
