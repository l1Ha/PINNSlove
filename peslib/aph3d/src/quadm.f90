SUBROUTINE quadm (nel, xmx, cm, xx)
Use FileUnits_Module
USE pass_Module
USE todim_Module
USE quad_Module
USE Numbers_Module
!            I N P U T    A R G U M E N T S
!  nel
!  xmx
!  cm
!  xx
IMPLICIT NONE
INTEGER nel, iintp, i, j, k, lx, ly, kl, ks
REAL(Kind=WP_Kind) phi, r, s, wt, det, theta, fac
REAL(Kind=WP_Kind) cm(171), xmx(27), d(18), xx(2,9), h(9), p(2,9),  xj(2,2), c(4,4)
EXTERNAL funct2, ststl
!-----------------------------------------------------------------------
! Initialize to zero.
!-----------------------------------------------------------------------
iintp=0
DO i=1,171
   cm(i)=zero
ENDDO
DO i=1,18
   xmx(i)=zero
ENDDO
!-----------------------------------------------------------------------
! DO a 3-point Gauss_Legendre quadrature in each direction.
!-----------------------------------------------------------------------
DO lx=1,nquad
   r=xg(lx)
   DO ly=1,nquad
      s=xg(ly)
      wt=wgt(lx)*wgt(ly)
!-----------------------------------------------------------------------
!     find interpolation functions and jacobian
!-----------------------------------------------------------------------
      CALL funct2 (r, s, h, p, xj, det, xx, iintp)
!-----------------------------------------------------------------------
!     compute the radius at point (r,s)
!-----------------------------------------------------------------------
      theta=zero
      phi =zero
      DO k=1,9
         phi = phi +h(k)*xx(2,k)
         theta=theta+h(k)*xx(1,k)
      ENDDO
      de=sin(2.d0*theta)
      CALL ststl(c, theta, phi)
      de = de*gtot*gtot
      fac=wt*det*de
      DO i=1,9
         d(2*i-1)=h(i)
         d(2*i)=h(i)
      ENDDO
      kl=1
      DO i=1,18,2
         DO j=i,18,2
            cm(kl)=cm(kl)+d(i)*d(j)*fac
            kl=kl+2
         ENDDO
         kl=kl+18-i
      ENDDO
   ENDDO
ENDDO
kl=1
DO i=1,18,2
   ks=kl+18-i+1
   DO j=i,18,2
      cm(ks)=cm(kl)
      ks=ks+2
      kl=kl+2
   ENDDO
   kl=kl+18-i
ENDDO
RETURN
END
