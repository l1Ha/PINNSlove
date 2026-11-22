SUBROUTINE tred3rw(n,nv,a,d,e,e2,nmax,nvmax)
USE Numeric_Kinds_Module
INTEGER i,j,k,l,n,ii,iz,jk,nv,nmax,nvmax
REAL(Kind=WP_Kind) a(nvmax),d(nmax),e(nmax),e2(nmax)
REAL(Kind=WP_Kind) f,g,h,hh,scale
REAL(Kind=WP_Kind) sqrt,abs,sign
!  -------------------------------------------------------------------
! this routine is a translation of the algol procedure tred3,
! num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
! handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
! this routine reduces a REAL symmetric matrix, stored as
! a one-dimensional array, to a symmetric tridiagonal matrix
! using orthogonal similarity transformations.
! on input:
! n is the order of the matrix;
! nv must be set to the dimension of the array parameter a
! as declared in the calling program dimension statement;
! a contains the lower triangle of the REAL symmetric
! input matrix, stored row-wise as a one-dimensional
! array, in its first n*(n+1)/2 positions.
! on output:
! a contains information about the orthogonal
! transformations used in the reduction;
! d contains the diagonal elements of the tridiagonal matrix;
! e contains the subdiagonal elements of the tridiagonal
! matrix in its last n-1 positions.e(1) is set to zero;
! e2 contains the squares of the corresponding elements of e.
! e2 may coincide with e IF the squares are not needed.
! questions and comments should be directed to b. s. garbow,
! applied mathematics division, argonne national laboratory
!  -------------------------------------------------------------------
DO ii = 1, n
i = n + 1 - ii
l = i - 1
iz = (i * l) / 2
h = 0.0d0
scale = 0.0d0
IF(l < 1) GOTO 20
DO 10 k = 1, l
iz = iz + 1
d(k) = a(iz)
scale = scale + ABS(d(k))
   10 CONTINUE
IF(scale /= 0.0d0) GOTO 30
   20 e(i) = 0.0d0
e2(i) = 0.0d0
GOTO 80
   30 DO 40 k = 1, l
d(k) = d(k) / scale
h = h + d(k) * d(k)
   40 CONTINUE
e2(i) = scale * scale * h
f = d(l)
g = -sign(sqrt(h),f)
e(i) = scale * g
h = h - f * g
d(l) = f - g
a(iz) = scale * d(l)
IF(l == 1) GOTO 80
f = 0.0d0
DO 60 j = 1, l
g = 0.0d0
jk = (j * (j-1)) / 2
DO 50 k = 1, l
jk = jk + 1
IF(k > j) jk = jk + k - 2
g = g + a(jk) * d(k)
   50 CONTINUE
e(j) = g / h
f = f + e(j) * d(j)
   60 CONTINUE
hh = f / (h + h)
jk = 0
DO j = 1, l
   f = d(j)
   g = e(j) - hh * f
   e(j) = g
   DO k = 1, j
      jk = jk + 1
      a(jk) = a(jk) - f * e(k) - g * d(k)
   ENDDO
ENDDO
   80 d(i) = a(iz+1)
a(iz+1) = scale * sqrt(h)
ENDDO
RETURN
!  --------------***END-tred3***---------------------------------------
END
