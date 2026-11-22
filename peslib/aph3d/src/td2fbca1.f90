SUBROUTINE td2fbca1(t,chi,w,jdeg,nmax,peigv)
USE FileUnits_Module
USE Max_Module
USE Matrices_Module
USE CMatt_Module
!
!     construct transformation matrix between angles and cosine
!     functions of symmetry a1 by diagonalization of cos(2*theta) matrix.
!     first index of t is angles,
!     second index labels cosine functions. chi is filled
!     with dvr points; w with corresponding weights. note
!     that dimension of resulting matrices is (jdeg)x(jdeg)
!
IMPLICIT NONE
EXTERNAL MPrint
LOGICAL peigv
INTEGER nmax, jdeg, n, i, k, j
REAL(Kind=WP_Kind) t(nmax,jdeg), chi(nmax), w(nmax)
n=jdeg
!
!     construct the matrix of cos(2*theta) over appropriate functions
!
DO i=2,n
   IF(i/=1)cmat(i,i-1)=0.5d0
   IF(i/=n)cmat(i,i+1)=0.5d0
ENDDO
cmat(1,2)=1./sqrt(2.d0)
cmat(2,1)=1./sqrt(2.d0)
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
WRITE(Out_Unit,40)
!
!     replace chi(i) by cos(acos(chi(i))/2) so that we have correct
!     chi values to use in later SUBROUTINEs
!
WRITE(dvr11,*) 'nchi = ',n
DO 16 i=1,n
   chi(i)=cos(acos(chi(i))/2.d0)
   WRITE(dvr11,*) 'chi(i) = ',acos(chi(i))
 16   CONTINUE
!
WRITE(Out_Unit,60)
!
IF(peigv) CALL printm(mprint,8,"t   ",t       ,n     ,n     ,maxang)
!
 40   FORMAT(/5x,'chi dvr matrix eigenvalues and angles'/)
 60   FORMAT(/5x,'chi dvr points and angles'/)
RETURN
END
