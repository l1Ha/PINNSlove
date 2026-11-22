SUBROUTINE tql2rw(nm,n,d,e,z,ierr,nmax,nvmax)
USE Numeric_Kinds_Module
INTEGER i,j,k,l,m,n,ii,l1,nm,mml,ierr,nmax,nvmax
REAL(Kind=WP_Kind) d(nmax),e(nmax),z(nmax,nmax)
REAL(Kind=WP_Kind) b,c,f,g,h,p,r,s,machep
REAL(Kind=WP_Kind) sqrt,abs,sign
!  -------------------------------------------------------------------
! machep is a machine dependent parameter specifying
! the relative precision of floating point arithmetic.
! machep = 16.0d0**(-13) for long form arithmetic
! on s360
!  -------------------------------------------------------------------
!     DATA machep/z3410000000000000/
DATA machep/1.d-10/
ierr = 0
IF(n == 1) GOTO 160
DO i = 2, n
   e(i-1) = e(i)
ENDDO
f = 0.0d0
b = 0.0d0
e(n) = 0.0d0
DO 110 l = 1, n
j = 0
h = machep * (ABS(d(l)) + ABS(e(l)))
IF(b < h) b = h
!  -------------------------------------------------------------------
! look for small sub-diagonal element
!  -------------------------------------------------------------------
DO 20 m = l, n
IF(ABS(e(m)) <= b) GOTO 30
!  -------------------------------------------------------------------
! e(n) is always zero, so there is no exit
! through the bottom of the loop
!  -------------------------------------------------------------------
   20 CONTINUE
   30 IF(m == l) GOTO 100
   40 IF(j == 30) GOTO 150
j = j + 1
!  -------------------------------------------------------------------
! form shift
!  -------------------------------------------------------------------
l1 = l + 1
g = d(l)
p = (d(l1) - g) / (2.0d0 * e(l))
r = sqrt(p*p+1.0d0)
d(l) = e(l) / (p + sign(r,p))
h = g - d(l)
DO i = l1, n
   d(i) = d(i) - h
ENDDO
f = f + h
!  -------------------------------------------------------------------
! ql transformation
!  -------------------------------------------------------------------
p = d(m)
c = 1.0d0
s = 0.0d0
mml = m - l
!  -------------------------------------------------------------------
! for i=m-1 step -1 until l DO --
!  -------------------------------------------------------------------
DO 90 ii = 1, mml
i = m - ii
g = c * e(i)
h = c * p
IF(ABS(p) < ABS(e(i))) GOTO 60
c = e(i) / p
r = sqrt(c*c+1.0d0)
e(i+1) = s * p * r
s = c / r
c = 1.0d0 / r
GOTO 70
   60 c = p / e(i)
r = sqrt(c*c+1.0d0)
e(i+1) = s * e(i) * r
s = 1.0d0 / r
c = c * s
   70 p = c * d(i) - s * g
d(i+1) = h + s * (c * g + s * d(i))
!  -------------------------------------------------------------------
! form vector
DO 80 k = 1, n
!  -------------------------------------------------------------------
h = z(k,i+1)
z(k,i+1) = s * z(k,i) + c * h
z(k,i) = c * z(k,i) - s * h
   80 CONTINUE
   90 CONTINUE
e(l) = s * p
d(l) = c * p
IF(ABS(e(l)) > b) GOTO 40
  100 d(l) = d(l) + f
  110 CONTINUE
!  -------------------------------------------------------------------
! order eigenvalues and eigenvectors
!  -------------------------------------------------------------------
DO 140 ii = 2, n
i = ii - 1
k = i
p = d(i)
DO 120 j = ii, n
IF(d(j) >= p) GOTO 120
k = j
p = d(j)
  120 CONTINUE
IF(k == i) GOTO 140
d(k) = d(i)
d(i) = p
DO 130 j = 1, n
p = z(j,i)
z(j,i) = z(j,k)
z(j,k) = p
  130 CONTINUE
  140 CONTINUE
GOTO 160
!  -------------------------------------------------------------------
! set error -- no convergence to an
! eigenvalue after 30 iterations
!  -------------------------------------------------------------------
  150 ierr = l
  160 RETURN
!  --------------***END-tql2***---------------------------------------
END
