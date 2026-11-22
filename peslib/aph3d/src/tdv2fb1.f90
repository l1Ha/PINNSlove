SUBROUTINE tdv2fb1(t,chi,w,jdeg,mord,nmax,peigv)
USE Numbers_Module
USE FileUnits_Module
USE Max_Module
USE Matrices_Module
USE CMatt_Module
!
!     construct transformation matrix between angles and rotational
!     functions. first index of t is angles, second index labels
!     rotational functions. chi is filled with dvr points; w with
!     corresponding weights. note that dimension of resulting matrices
!     is (jdeg-mord+1)x(jdeg-mord+1)
!
IMPLICIT NONE
LOGICAL peigv
EXTERNAL mprint
INTEGER nmax, jdeg, mord, n, l, i, j, k
REAL(Kind=WP_Kind) t(nmax,jdeg-mord+1), chi(nmax), w(nmax), pj(maxang)
n=jdeg-mord+1
l=mord
IF(n>maxang.or.n>nmax)THEN
   WRITE(Out_Unit,*)'ERROR in tdv2fb1: n>maxang.or.n>nmax',  n,maxang,nmax
   WRITE(Msg_Unit,*)'ERROR in tdv2fb1: n>maxang.or.n>nmax',  n,maxang,nmax
   STOP 'tdv2fb1'
ENDIF
!
!     construct the matrix of cos over the legendre polynomials
!
DO i=1,n
   j=i-1+l
   IF(i/=n)cmat(i,i+1)=sqrt((j+l+One)*(j-l+One)/((Two*j+One)*(Two*j+Three)))
   IF(i/=1)cmat(i,i-1)=sqrt((j+l)*(j-l)/((Two*j+One)*(Two*j-One)))
ENDDO
!
!     diagonalize cmat and put eigenvectors in t and eigenvalues in chi
!
k=0
DO i=1,n
   DO j=1,i
      k=k+1
      sd(k)=cmat(i,j)
   ENDDO
ENDDO
CALL tdiagrw(sd,chi,rdt,subd,n,maxang,maxangv)
DO i=1,n
   DO j=1,n
      t(i,j)=rdt(j,i)
   ENDDO
ENDDO
IF(peigv) CALL printm(mprint,8,"t   ",t       ,n     ,n     ,maxang)
!
WRITE(dvr11,*) 'ntheta = ',n
DO i=1,n
  IF(ABS(Chi(i))<=One)THEN
    WRITE(dvr11,*) 'theta(i) = ',acos(chi(i))/Two
  ELSE
    WRITE(Msg_Unit,*)'Error: argument of ACOS is out of range:',Chi(i) 
    WRITE(Out_Unit,*)'Error: argument of ACOS is out of range:',Chi(i)
   !STOP 'tdv2fb1'
  ENDIF
ENDDO
!
RETURN
END
