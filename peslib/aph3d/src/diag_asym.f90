SUBROUTINE diag_asym(n, w0, w1, w2, rhoval, nrho, read_p)
USE Numeric_Kinds_Module
IMPLICIT NONE

LOGICAL read_p
INTEGER n, nrho, i, j
REAL(Kind=WP_Kind) rhoval(nrho), w0(n,n), w1(n,n), w2(n,n), &
        rhof1, rhof2, rhof3, &
        xtrafac, f_cor, f_asy
COMMON /xtra/ xtrafac, f_cor(0:1000), f_asy(0:1000)
!----------------------------------------------------------------------
! V(rhocent)=0 therefore w1 contains only:
!          2*mu*rho**2*(A-B)
!
!      xtrafac=-1/(2mu)*jtot*(jtot+1)*(-1)**(jtot+parity)/4.0d0
!----------------------------------------------------------------------
IF(read_p)THEN
   rhof1=xtrafac/(rhoval(1)*rhoval(1))
   rhof3=xtrafac/(rhoval(3)*rhoval(3))
   rhof2=xtrafac/(rhoval(2)*rhoval(2))
   DO j=1,n
      DO i=j,n
         w0(i,j) = w0(i,j)+w1(i,j)*rhof1
         w2(i,j) = w2(i,j)+w1(i,j)*rhof3
         w1(i,j) = w1(i,j)*rhof2
         w0(j,i) = w0(i,j)
         w2(j,i) = w2(i,j)
         w1(j,i) = w1(i,j)
      ENDDO
   ENDDO
ENDIF

RETURN
END
