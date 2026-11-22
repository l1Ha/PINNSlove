SUBROUTINE Hermit (npt, xpt, wht, kase)
USE Numeric_Kinds_Module
USE FileUnits_Module
USE Numbers_Module
IMPLICIT NONE
REAL(Kind=WP_Kind) xpt(npt), wht(npt), fac, hpol(201)
INTEGER i, npt, kase
!---------------------------------------------------------------------
!
!  purpose
!     supplies the Gauss_Hermite abscissas and weights
!
!  description of parameters
!     npt-number of points desired
!     xpt-resultant array containing the Gauss_Hermite abscissas.
!     wht-resultant array containing the Gauss_Hermite weights
!     kase=0 IF integrating exp(-x**2)*f(x)*dx
!     kase=1 IF integrating f(x)*dx
!---------------------------------------------------------------------
IF(npt>2000)THEN   ! PORT: enlarged from 200
   WRITE(Out_Unit,*)'npt>200', npt
   STOP 'hermit'
ENDIF

CALL hrmzro (npt, xpt)

fac = 1.d0
DO i = 1, npt
   fac = fac*sqrt(REAL(i,WP_Kind))
ENDDO
fac = fac*sqrt(2.d0**(npt-1))*pi**0.25d0/npt

DO i = 1, (npt+1)/2
   CALL hermit_fun(xpt(i), npt-1, hpol, npt)
   wht(i) =  (fac/hpol(npt))**2
ENDDO

DO i = 1, (npt+1)/2
   xpt(npt-i+1) = xpt((npt+1)/2-i+1)
   wht(npt-i+1) = wht((npt+1)/2-i+1)
ENDDO

DO i = 1, npt/2
   xpt(i) = -xpt(npt+1-i)
   wht(i) =  wht(npt+1-i)
ENDDO

IF(kase==1)THEN
   DO i=1,npt
      wht(i)=wht(i)*exp(xpt(i)**2)
   ENDDO
ENDIF

RETURN
ENDSUBROUTINE Hermit
