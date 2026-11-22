SUBROUTINE Glegen (npt, xpt, wht, a, b)
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
INTEGER i, npt
REAL(Kind=WP_Kind) xpt(npt), wht(npt), pol(0:2000), a, b, t1, t2   ! PORT: pol enlarged from 0:200
!---------------------------------------------------------------------
!
!  purpose
!     supplies the Gauss_Legendre abscissas and weights
!
!  description of parameters
!     npt-number of points desired
!     xpt-resultant array containing the Gauss_Legendre abscissas.
!     wht-resultant array containing the Gauss_Legendre weights
!---------------------------------------------------------------------
IF(npt>2000)THEN
   WRITE(Out_Unit,*)'npt>200', npt
   STOP 'glegen'
ENDIF

CALL GlgnZro (npt, xpt)

DO i = 1, (npt+1)/2
   CALL Pn(xpt(i), npt, pol)
   wht(i) = 2.d0*(1.d0-xpt(i)**2)/(npt*(-xpt(i)*pol(npt)+pol(npt-1)) )**2
ENDDO

DO i = 1, (npt+1)/2
   xpt(npt-i+1) = xpt((npt+1)/2-i+1)
   wht(npt-i+1) = wht((npt+1)/2-i+1)
ENDDO

DO i = 1, npt/2
   xpt(i) = -xpt(npt+1-i)
   wht(i) =  wht(npt+1-i)
ENDDO

t1=(b-a)/2.0d0
t2=(b+a)/2.0d0
DO i=1,npt
   xpt(i)=xpt(i)*t1+t2
   wht(i)=wht(i)*t1
ENDDO


RETURN
ENDSUBROUTINE Glegen
