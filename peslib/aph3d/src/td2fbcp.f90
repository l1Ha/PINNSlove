SUBROUTINE td2fbcp(t,chi,w,jdeg,nmax,peigv)
USE FileUnits_Module
USE Max_Module
USE Matrices_Module
USE CMatt_Module
!
!     construct transformation matrix between angles and cosine
!     functions of symmetry a(p)rime by diagonalization of cos(2*theta)
!     matrix.
!     first index of t is angles,
!     second index labels cosine functions. chi is filled
!     with dvr points; w with corresponding weights. note
!     that dimension of resulting matrices is (jdeg)x(jdeg)
!
IMPLICIT NONE
LOGICAL peigv
EXTERNAL mprint
INTEGER nmax, jdeg, n, n2, n4, i, j, k
REAL(Kind=WP_Kind) t(nmax,jdeg),chi(nmax),w(nmax)
REAL(Kind=WP_Kind) r2inv, pi2
n=jdeg
n2=n/2
n4=n/4
!
!     construct the matrix of cos(2*theta) over cos(2*i*theta)
!
DO i=1,n2
   DO j=1,n2
      cmat(i,j)=0.0d0
   ENDDO
ENDDO
!
DO i=2,n4
   IF(i/=1)cmat(i,i-1)=0.5d0
   IF(i/=n4)cmat(i,i+1)=0.5d0
ENDDO
r2inv=1./sqrt(2.d0)
cmat(1,2)=r2inv
cmat(2,1)=r2inv
!
!     construct the matrix of cos(2*theta) over sin(2*i*theta)
!
DO i=n4+1,n2
   IF(i/=(n4+1))cmat(i,i-1)=0.5d0
   IF(i/=n2)cmat(i,i+1)=0.5d0
ENDDO
!
IF(peigv) CALL printm(mprint,8,"cmat",cmat    ,jdeg  ,jdeg  ,maxang)
!
!     diagonalize cmat and put eigenvectors in t and eigenvalues in chi
!
k=0
DO i=1,n2
   DO j=1,i
      k=k+1
      sd(k)=cmat(i,j)
   ENDDO
ENDDO
CALL tdiagrw(sd,chi,rdt,subd,n2,maxang,maxangv)
DO i=1,n2
   DO j=1,n2
      t(i,j)=r2inv*rdt(j,i)
   ENDDO
ENDDO
!
!     replace chi(i) by cos(acos(chi(i))/2) so that we have correct
!     chi values to use in later SUBROUTINEs
!
DO i=1,n2
   w(i)=cos(acos(chi(i)))
ENDDO
!
!     construct the matrix of cos(2*theta) over cos(2*i*theta)
!
DO i=1,n2
   DO j=1,n2
      cmat(i,j)=0.0d0
   ENDDO
ENDDO
!
DO i=2,n4
   cmat(i,i+n4)=0.5d0
   cmat(i+n4,i)=0.5d0
   IF( i/=(n4-1) .AND. i/=n4 )THEN
     cmat(i,i+2+n4)=-0.5d0
     cmat(i+2+n4,i)=-0.5d0
   ENDIF
ENDDO
cmat(1,1+n4)=r2inv
cmat(1+n4,1)=r2inv
!
IF(peigv) CALL printm(mprint,8,"cmat",cmat    ,jdeg  ,jdeg  ,maxang)
!
!     diagonalize cmat and put eigenvectors in t and eigenvalues in chi
!
k=0
DO i=1,n2
   DO j=1,i
      k=k+1
      sd(k)=cmat(i,j)
   ENDDO
ENDDO
CALL tdiagrw(sd,chi,rdt,subd,n2,maxang,maxangv)
DO i=1,n2
   DO j=1,n2
      t(i+n2,j+n2)=r2inv*rdt(j,i)
      t(i+n2,j)=0.0d0
      t(j,i+n2)=0.0d0
   ENDDO
ENDDO
!
WRITE(Out_Unit,40)
!
!     replace chi(i) by cos(acos(chi(i))/2) so that we have correct
!     chi values to use in later SUBROUTINEs
!
pi2=acos(-1.0d0)*0.5d0
DO i=1,n2
   w(i+n2)=cos(asin(chi(i))+pi2)
ENDDO
!
!
WRITE(dvr11,*) 'nchi   = ',n
DO i=1,n
  chi(i)=w(i)
  WRITE(dvr11,*) 'chi(i)   = ',acos(chi(i))
ENDDO
!
WRITE(Out_Unit,60)
!
IF(peigv) CALL printm(mprint,8,"t   ",t       ,n     ,n     ,maxang)
!
 40   FORMAT(/5x,'chi dvr matrix eigenvalues and angles'/)
 60   FORMAT(/5x,'chi dvr points and angles'/)
RETURN
END
