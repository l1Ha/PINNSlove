SUBROUTINE trbak3rw(nm,n,nv,a,m,z,nmax,nvmax)
USE Numeric_Kinds_Module
INTEGER i,j,k,l,m,n,ik,iz,nm,nv,nmax,nvmax
REAL(Kind=WP_Kind) a(nvmax),z(nmax,nmax)
REAL(Kind=WP_Kind) h,s
!  -------------------------------------------------------------------
!     this routine is a translation of the algol procedure trbak3,
!     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
!     this routine forms the eigenvectors of a REAL symmetric
!     matrix by back transforming those of the corresponding
!     symmetric tridiagonal matrix determined bytred3.
!     on input:
!     nm must be set to the row dimension of two-dimensional
!     array parameters as declared in the calling program
!     dimension statement;
!     n is the order of the matrix;
!     nv must be set to the dimension of the array parameter a
!     as declared in the calling program dimension statement;
!     a contains information about the orthogonal transformations
!     used in the reduction bytred3in its first
!     n*(n+1)/2 positions;
!     m is the number of eigenvectors to be back transformed;
!     z contains the eigenvectors to be back transformed
!     in its first m columns.
!     on output:
!     z contains the transformed eigenvectors
!     in its first m columns.
!     note that trbak3 preserves vector euclidean norms.
!     questions and comments should be directed to b. s. garbow,
!     applied mathematics division, argonne national laboratory
!  -------------------------------------------------------------------
IF(m == 0) GOTO 50
IF(n == 1) GOTO 50
DO 40 i = 2, n
l = i - 1
iz = (i * l) / 2
ik = iz + i
h = a(ik)
IF(h == 0.0d0) GOTO 40
DO 30 j = 1, m
s = 0.0d0
ik = iz
DO 10 k = 1, l
ik = ik + 1
s = s + a(ik) * z(k,j)
   10 CONTINUE
s = (s / h) / h
ik = iz
DO 20 k = 1, l
ik = ik + 1
z(k,j) = z(k,j) - s * a(ik)
   20 CONTINUE
   30 CONTINUE
   40 CONTINUE
   50 RETURN
!  --------------***END-trbak3***---------------------------------------
END
