SUBROUTINE td2fbcpp(t,chi,w,jdeg,nmax,peigv)
USE Max_Module
USE Matrices_Module
USE CMatt_Module
!
!     construct transformation matrix between angles and cosine
!     functions of symmetry a(p)rime(p)rime by diagonalization of
!     cos(2*theta) matrix.
!     first index of t is angles,
!     second index labels cosine functions. chi is filled
!     with dvr points; w with corresponding weights. note
!     that dimension of resulting matrices is (jdeg)x(jdeg)
!
IMPLICIT NONE
LOGICAL peigv
EXTERNAL mprint
INTEGER nmax, jdeg, n, n2, k, j, i
REAL(Kind=WP_Kind) t(nmax,jdeg), chi(nmax), w(nmax)
n=jdeg
n2=n/2
!
!     construct the matrix of cos(2*theta) over sin((2*i+1)*theta)
!
DO i=1,n2
   IF(i/=1)cmat(i,i-1)=0.5d0
   IF(i/=n2)cmat(i,i+1)=0.5d0
ENDDO
cmat(1,1)=-0.5d0
!
!     construct the matrix of cos(2*theta) over cos((2*i+1)*theta)
!
DO i=n2+1,n
   IF(i/=(n2+1))cmat(i,i-1)=0.5d0
   IF(i/=n)cmat(i,i+1)=0.5d0
ENDDO
cmat(n2+1,n2+1)=0.5d0
!
IF(peigv) CALL printm(mprint,8,"cmat", cmat    ,jdeg  ,jdeg  ,maxang)
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
!
!     replace chi(i) by cos(acos(chi(i))/2) so that we have correct
!     chi values to use in later SUBROUTINEs
!
DO i=1,n
   chi(i)=cos(acos(chi(i))/2.d0)
ENDDO
IF(peigv) CALL printm(mprint,8,"t   ",t       ,n     ,n     ,maxang)
RETURN
END
