SUBROUTINE laguer (xpt, wht, npt, kase)
USE FileUnits_Module
USE scalag_Module
IMPLICIT NONE
REAL(Kind=WP_Kind) xpt(*), wht(*), fac, hpol(0:200), sum, abs
INTEGER i, npt, kase
!---------------------------------------------------------------------
!
!  purpose
!     supplies the Gauss_Laguerre abscissas and weights
!
!  description of PARAMETERs
!     npt-number of points desired
!     xpt-resultant array containing the Gauss_Laguerre abscissas.
!     wht-resultant array containing the Gauss_Laguerre weights
!     kase=0 IF integrating exp(-x)*f(x)*dx
!     kase=1 IF integrating f(x)*dx
!
!  SUBROUTINE and functions used
!     NONE
!
!---------------------------------------------------------------------
IF(npt>2000)THEN   ! PORT: enlarged from 200
   WRITE(Out_Unit,*)'npt>200', npt
   STOP 'laguer'
ENDIF

CALL lagzro (npt, xpt)

fac = 1.d0/npt

DO 20 i = 1, npt
   CALL laguer_fun(xpt(i), npt, hpol, npt+1)
   wht(i) =  xpt(i)*(fac/(hpol(npt)-hpol(npt-1)))**2
   20 CONTINUE

!
! Check for Accuracy
!
sum = 0.d0
DO 34 i=1,npt
  sum = sum + wht(i)
   34 CONTINUE

IF( ABS(sum-1.d0)>1.d-12)THEN
    WRITE(Out_Unit,*)'Error in weights'
    STOP 'laguer'
ENDIF

sum =0.d0
DO 35 i =1,npt
   sum = sum +wht(i)*xpt(i)
   35 CONTINUE

IF( ABS(sum-1.d0)>1.d-12)THEN
    WRITE(Out_Unit,*)'Error in weights'
    STOP 'laguer'
ENDIF

DO 40 i=1,npt
   xpt(i)=xpt(i)/scalagr
   wht(i)=wht(i)/scalagr
   IF(kase==1)wht(i)=wht(i)*exp(scalagr*xpt(i))
   40 CONTINUE

RETURN
END
