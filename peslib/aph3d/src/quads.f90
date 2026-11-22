SUBROUTINE quads (nel, b, s, yz)
USE FileUnits_Module
USE todim_Module
USE quad_Module
!
!   $RCSfile: quads.f,v $ $Revision: 1.15 $
!   $Date: 89/11/16 11:03:20 $
!   $State: Stable $
!


!            P U R P O S E   O F    S U B R O U T I N E


!            I N P U T    A R G U M E N T S
!
!    nel
!    b
!    s
!    yz


IMPLICIT NONE



!            I N T E G E R S
INTEGER lx, ly, kl, ks, i, ii, j, k, nel

!            R E A L S
REAL(Kind=WP_Kind) b, s, yz, wt, det, e1, e2, xbar, theta, phi, fac, db, d

!            D I M E N S I O N S
DIMENSION b(4,18), s(378), yz(*), db(4),   d(4,4)

!            E X T E R N A L S
EXTERNAL deriq, ststl
!-----------------------------------------------------------------------
! DO a 3-point Gauss_Legendre quadrature in each direction
!-----------------------------------------------------------------------
DO lx=1,nquad
   e1=xg(lx)
   DO ly=1,nquad
      e2=xg(ly)
      wt=wgt(lx)*wgt(ly)
      CALL deriq (nel, yz, b, det, e1, e2, xbar, theta, phi)
      CALL ststl (d, theta, phi)
      fac=wt*xbar*det
      kl=1
      DO j=1,18,2
         DO k=1,3
            db(k)=d(k,1)*b(1,j)+d(k,3)*b(3,j)
            db(k)=db(k)*fac
         ENDDO
         DO i=j,18,2
            s(kl)=s(kl)+b(1,i)*db(1)+b(3,i)*db(3)
            kl=kl+1
            s(kl)=s(kl)+b(2,i+1)*db(2)+b(3,i+1)*db(3)
            kl=kl+1
         ENDDO
         kl=kl+18-j
      ENDDO
      kl=18+1
      DO j=2,18,2
         DO k=1,3
            db(k)=d(k,2)*b(2,j)+d(k,3)*b(3,j)
            db(k)=db(k)*fac
         ENDDO
         ks=kl
         DO i=j,18,2
            s(ks)=s(ks)+b(2,i)*db(2)+b(3,i)*db(3)
            ks=ks+2
         ENDDO
         IF(j-18<0)THEN
            k=j+1
            ks=kl+1
            DO ii=k,18,2
               s(ks)=s(ks)+b(1,ii)*db(1)+b(3,ii)*db(3)
               ks=ks+2
            ENDDO
         ENDIF
         kl=kl+2*18-2*j+1
      ENDDO
      kl=1
      DO j=1,18,2
         db(1)=d(1,4)*b(4,j)*fac
         db(2)=d(2,4)*b(4,j)*fac
         db(3)=d(3,4)*b(4,j)*fac
         db(4)=d(4,1)*b(1,j)+d(4,3)*b(3,j)+d(4,4)*b(4,j)
         db(4)=db(4)*fac
         DO i=j,18,2
            s(kl)=s(kl)+b(1,i)*db(1)+b(3,i)*db(3)+b(4,i)*db(4)
            kl=kl+1
            s(kl)=s(kl)+b(2,i+1)*db(2)+b(3,i+1)*db(3)
            kl=kl+1
         ENDDO
         kl=kl+18-j
      ENDDO
      kl=18+1
      DO j=2,18,2
         db(4)=d(4,2)*b(2,j)+d(4,3)*b(3,j)
         db(4)=db(4)*fac
         DO i=j,18
            s(kl)=s(kl)+b(4,i)*db(4)
            kl=kl+1
         ENDDO
         kl=kl+18-j
      ENDDO
   ENDDO
ENDDO

RETURN
END
