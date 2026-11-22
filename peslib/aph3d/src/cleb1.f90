REAL(Kind=WP_Kind) FUNCTION cleb1 (ja, jb, jc, ma, mb, mc)
!
USE FileUnits_Module
!  REAL(Kind=WP_Kind) FUNCTION cleb1 (ja, jb, jc, ma, mb, mc)
IMPLICIT NONE
INTEGER ja, jb, jc, ma, mb, mc, j1, j2, j3, m1, m2, m3, jx, mx
INTEGER j1max, j1min, j2old, j3old, m1old, m2old, ncoef, mid
INTEGER n2, mid1, i, ip1, im1, j, j3m3
!----------------------------------------------------------------------
!  this program was written by g. a. parker 1/3/77.  for a
!  description fo the algorithm used see klaus schulten and r. g.
!  gordon  jmp 16,1961,(1975).
!  this program is most efficient when a series of clebsch-
!  gordan coefficients are desired
!----------------------------------------------------------------------
REAL(Kind=WP_Kind) f(500), sj2mj3, sj2j31, xj2m, xj3m, xm1, y
REAL(Kind=WP_Kind) sj1, xj1, xm, x, xj1p1, sj1p1, g, coef, sum, sig, sm1
DATA j2old, j3old, m1old, m2old / -1000000, -1000000, -10000000, -1000000 /
!----------------------------------------------------------------------
!  statement functions used in recurrence relations.
!----------------------------------------------------------------------
!a(sj1)=sqrt((sj1-sj2mj3)*(sj2j31-sj1)*(sj1-sm1))
!b(xj1)=-(2.0d0*xj1+1)*(xj2m-xj3m-xj1*(xj1+1)*xm1)
!----------------------------------------------------------------------
!  check triangle inequalities.
!----------------------------------------------------------------------
IF(iabs(mc) > jc) GOTO 130
IF(iabs(mb) > jb) GOTO 130
IF(iabs(ma) > ja) GOTO 130
IF(mc /= ma+mb) GOTO 130
!----------------------------------------------------------------------
!  calculate maximum and minimum allowed values of j1
!----------------------------------------------------------------------
j1=ja
j2=jb
j3=jc
m1=ma
m2=mb
m3=-mc
jx=iabs(j2-j3)
mx=iabs(m1)
j1max=j2+j3
j1min=max0(jx,mx)
IF(j1max-j1min+1 >= 500) GOTO 140
IF(j1 > j1max .or. j1 < j1min) GOTO 130
!----------------------------------------------------------------------
!  check to see if the clebsch-gordan coefficient has already
!  determined.
!----------------------------------------------------------------------
IF(j2old /= j2) GOTO 10
IF(j3old /= j3) GOTO 10
IF(m1old /= m1) GOTO 10
IF(m2old /= m2) GOTO 10
GOTO 100
!----------------------------------------------------------------------
!  coefficients have not been previously determined.
!----------------------------------------------------------------------
10   j2old=j2
j3old=j3
m1old=m1
m2old=m2
!----------------------------------------------------------------------
!  determine constants for statement functions.
!----------------------------------------------------------------------
sj2mj3=jx*jx
sj2j31=(j1max+1)*(j1max+1)
sm1=mx*mx
xj2m=j2*(j2+1)*m1
xj3m=j3*(j3+1)*m1
xm1=-m1-m2-m2
!----------------------------------------------------------------------
!  determine where the middle coefficient is located and the
!  number of coefficients.
!----------------------------------------------------------------------
ncoef=j1max-j1min+1
mid=ncoef/2+1
IF(ncoef <= 2) GOTO 110
!----------------------------------------------------------------------
!  set initial values for the coefficietts f(1) and f(ncoef).
!----------------------------------------------------------------------
f(1)=1.d0
f(ncoef)=1.d0
!----------------------------------------------------------------------
!  use equations 16 and 17 from the reference to determine
!  f(2) and f(ncoef-1)
!----------------------------------------------------------------------
IF(j1min /= 0) GOTO 20
xm=2*m2
x=(2*j2+1)**2-1
f(2)=-f(1)*xm/sqrt(x)
GOTO 30
20   xj1=j1min
xj1p1=j1min+1
sj1p1=xj1p1*xj1p1
f(2)=f(1)*b(xj1)/(xj1*a(sj1p1))
30   g=f(2)
n2=ncoef-1
xj1=j1max
sj1=xj1*xj1
f(n2)=f(ncoef)*b(xj1)/((xj1+1.0d0)*a(sj1))
mid1=mid-1
IF(ncoef == 3) GOTO 50
!----------------------------------------------------------------------
!  start recursive evaluation of the coefficients using
!  equation 6'.
!----------------------------------------------------------------------
DO 40 i=2, mid1
    ip1=i+1
    im1=i-1
    xj1=j1min+i-1
    xj1p1=xj1+1.0d0
    sj1=xj1*xj1
    sj1p1=xj1p1*xj1p1
    g=(b(xj1)*f(i)-xj1p1*a(sj1)*f(im1))/(xj1*a(sj1p1))
    n2=ncoef-i+1
    IF(n2 == mid) GOTO 50
    f(ip1)=g
    xj1=j1max-i+1
    xj1p1=xj1+1.0d0
    ip1=n2+1
    im1=n2-1
    sj1=xj1*xj1
    sj1p1=xj1p1*xj1p1
    f(im1)=(-xj1*a(sj1p1)*f(ip1)+b(xj1)*f(n2))/(xj1p1*a(sj1))
40   CONTINUE
!----------------------------------------------------------------------
!  match at the midpoint.
!----------------------------------------------------------------------
50   IF(ABS(g) > 1.d-8) GOTO 60
i=mid
im1=mid1
xj1=j1min+mid1
xj1p1=xj1+1.0d0
sj1=xj1*xj1
sj1p1=xj1p1*xj1p1
g=(b(xj1)*f(i)-xj1p1*a(sj1)*f(im1))/(xj1*a(sj1p1))
mid1=mid
mid=mid+1
60   coef=f(mid)/g
DO 70 i=1, mid1
    f(i)=coef*f(i)
70   CONTINUE
!----------------------------------------------------------------------
!  find magnitude of the normalization coefficient.
!----------------------------------------------------------------------
sum=0.0d0
x=j1min*2-1
DO 80 i=1, ncoef
    x=x+2.0d0
    sum=sum+x*f(i)*f(i)
80   CONTINUE
coef=sqrt(1.0d0/sum)
!----------------------------------------------------------------------
!  determine sign of the normalization coefficient.
!----------------------------------------------------------------------
sig=1.0d0
j=j2-j3-m1
IF((j/2)*2 /= j) sig=-1.0d0
coef=sign(coef, sig)
!----------------------------------------------------------------------
!  normalize to obtain the clebsch-gordan coefficients.
!----------------------------------------------------------------------
y=2*j3+1
y=sqrt(y)
j3m3=j3+m3
IF((j3m3/2)*2 /= j3m3) y=-y
DO 90 i=1, ncoef
    f(i)=coef*f(i)*y
90   CONTINUE
100  i=j1-j1min+1
cleb1=f(i)
RETURN
110  IF(ncoef == 2) GOTO 120
!----------------------------------------------------------------------
!  special case for only 1 coefficient.
!----------------------------------------------------------------------
y=2*j3+1
x=2*j1+1
j=j2-m2
sig=1.0d0
IF((j/2)*2 /= j) sig=-1.0d0
f(1)=sqrt(y/x)*sig
GOTO 100
!----------------------------------------------------------------------
!  special case for only 2 coefficients.
!----------------------------------------------------------------------
120  f(1)=-1.0d0
xj1=j1min
xj1p1=j1min+1
sj1p1=xj1p1*xj1p1
f(2)=f(1)*b(xj1)/(xj1*a(sj1p1))
x=j1min*2+1
sum=x*f(1)*f(1)+(x+2.0d0)*f(2)*f(2)
coef=sqrt(1.0d0/sum)
sig=1.0d0
j=j2-j3-m1
IF((j/2)*2 /= j) sig=-1.0d0
sig=sig*f(2)
coef=sign(coef, sig)
y=2*j3+1
y=sqrt(y)
j3m3=j3+m3
IF((j3m3/2)*2 /= j3m3) y=-y
f(1)=coef*y*f(1)
f(2)=coef*y*f(2)
GOTO 100
!----------------------------------------------------------------------
!  triangle inequalities were not satisfied.
!----------------------------------------------------------------------
130  cleb1=0.0d0
RETURN
140  WRITE(Out_Unit,150) j1min, j1max
!----------------***end-cleb1***---------------------------------------
!
150  FORMAT(//, 'error in cleb. the number of j1 values exceeds',   &
    ''' dimensions.', 5x, 'j1min=', 5x, i5, 'j1max=', 5x, i5)
CONTAINS
   REAL(Kind=WP_Kind) FUNCTION A(sj1)
      USE Numeric_Kinds_Module
      REAL(Kind=WP_Kind) sj1
      a = sqrt((sj1-sj2mj3)*(sj2j31-sj1)*(sj1-sm1))
   END FUNCTION A
   REAL(Kind=WP_Kind) FUNCTION B(xj1)
      USE Numeric_Kinds_Module
      REAL(Kind=WP_Kind) xj1
      b =-(2.0d0*xj1+1)*(xj2m-xj3m-xj1*(xj1+1)*xm1)
   END FUNCTION B
END FUNCTION Cleb1
