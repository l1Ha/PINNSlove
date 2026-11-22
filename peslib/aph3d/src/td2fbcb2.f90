SUBROUTINE td2fbcb2(t,chi,w,jdeg,nmax,peigv)
USE Numeric_Kinds_Module
USE Max_Module
USE Matrices_Module
USE CMatt_Module
!
!     construct transformation matrix between angles and cosine
!     functions of symmetry b2 by diagonalization of cos(2*theta) matrix.
!     first index of t is angles,
!     second index labels cosine functions. chi is filled
!     with dvr points; w with corresponding weights. note
!     that dimension of resulting matrices is (jdeg)x(jdeg)
!
LOGICAL peigv
INTEGER NMax, JDeg, N, I, K, J
EXTERNAL mprint
REAL(Kind=WP_Kind) t(nmax,jdeg),chi(nmax),w(nmax)
n=jdeg
!
!     construct the matrix of cos(2*theta) over appropriate functions
!
DO 3 i=1,n
   IF(i/=1)cmat(i,i-1)=0.5d0
   IF(i/=n)cmat(i,i+1)=0.5d0
 3    CONTINUE
cmat(1,1)=0.5d0
IF(peigv) CALL printm(mprint,8,"cmat",cmat    ,jdeg  ,jdeg  ,maxang)
!
!     diagonalize cmat and put eigenvectors in t and eigenvalues in chi
!
k=0
DO 12 i=1,n
   DO 18 j=1,i
      k=k+1
      sd(k)=cmat(i,j)
 18      CONTINUE
 12   CONTINUE
CALL tdiagrw(sd,chi,rdt,subd,n,maxang,maxangv)
DO 14 i=1,n
   DO 15 j=1,n
      t(i,j)=rdt(j,i)
 15      CONTINUE
 14   CONTINUE
!
!     replace chi(i) by cos(acos(chi(i))/2) so that we have correct
!     chi values to use in later SUBROUTINEs
!
DO 16 i=1,n
   chi(i)=cos(acos(chi(i))/2.d0)
 16   CONTINUE
IF(peigv) CALL printm(mprint,8,"t   ",t       ,n     ,n     ,maxang)
RETURN
END
