REAL(Kind=WP_Kind) FUNCTION c6j (j1, j2, j3, l1, l2, l3)
USE FileUnits_Module
IMPLICIT NONE
!----------------------------------------------------------------------
!  this program was written by g. a. parker 1/3/77.  for a
!  description of the algorithm used see klaus schulten and r. g.
!  gordon  jmp 16,1961,(1975).
!  programme modified for the computation of a series of 6j-coefficients
!  by g. d. kohler 7/13/82.
!----------------------------------------------------------------------
INTEGER j, j1, j2, j3, l1, l2, l3, jx, lx, j1min, j1max, i, ncoef, mid
INTEGER j2old, j3old, l1old, l2old, l3old, n2, mid1, ip1, im1
REAL(Kind=WP_Kind) h(0:500)
REAL(Kind=WP_Kind) c0, c1, sj2mj3, sj2j31, sl2ml3, sl2l31, sj1
REAL(Kind=WP_Kind) y, x, xj1, sj1p1, xj1p1, g, sum, coef
DATA j2old, j3old, l1old, l2old / -1000000, -1000000, -10000000, -1000000 /
!----------------------------------------------------------------------
!  calculate maximum and minimum allowed values of j1
!----------------------------------------------------------------------
jx=iabs(j2-j3)
lx=iabs(l2-l3)
j1max=min(j2+j3,l2+l3)
j1min=max(jx,lx)
!----------------------------------------------------------------------
!  check triangle inequalities.
!----------------------------------------------------------------------
IF(j2 > l1+l3) GOTO 10
IF(l1 > j2+l3) GOTO 10
IF(l3 > l1+j2) GOTO 10
IF(j3 > l1+l2) GOTO 10
IF(l1 > j3+l2) GOTO 10
IF(l2 > j3+l1) GOTO 10
GOTO 30
10   DO i=j1min, j1max
    h(i)=0.d0
ENDDO
GOTO 120
30   IF(j2old /= j2) GOTO 40
IF(j3old /= j3) GOTO 40
IF(l1old /= l1) GOTO 40
IF(l2old /= l2) GOTO 40
IF(l3old /= l3) GOTO 40
GOTO 120
!----------------------------------------------------------------------
!  coefficients have not been previously determined.
!----------------------------------------------------------------------
40   j2old=j2
j3old=j3
l1old=l1
l2old=l2
l3old=l3
!----------------------------------------------------------------------
!  determine constants for statement functions.
!----------------------------------------------------------------------
sj2mj3=jx*jx
sj2j31=(j2+j3+1)*(j2+j3+1)
sl2ml3=(l2-l3)*(l2-l3)
sl2l31=(l2+l3+1)*(l2+l3+1)
c1=j2*(j2+1)+j3*(j3+1)+l2*(l2+1)+l3*(l3+1)-2*l1*(l1+1)
c0=(l2*(l2+1)-l3*(l3+1))*(j2*(j2+1)-j3*(j3+1))
!----------------------------------------------------------------------
!  determine where the middle coefficient is located and the
!  number of coefficients.
!----------------------------------------------------------------------
ncoef=j1max-j1min+1
mid=(j1max+j1min+1)/2
IF(ncoef <= 2) GOTO 130
!----------------------------------------------------------------------
!  set initial values for the coefficients h(j1min) and h(j1max).
!----------------------------------------------------------------------
h(j1min)=1.d0
h(j1max)=1.d0
!----------------------------------------------------------------------
!  use equations 24 and 25 from the reference to determine
!  h(j1min+1) and h(j1max-1)
!----------------------------------------------------------------------
IF(j1min /= 0) GOTO 50
y=2*(l1*(l1+1)-l2*(l2+1)-j2*(j2+1))
x=((2*j2+1)**2-1)*((2*l2+1)**2-1)
h(j1min+1)=h(j1min)*y/sqrt(x)
GOTO 60
50   xj1=j1min
xj1p1=j1min+1.d0
sj1p1=xj1p1*xj1p1
h(j1min+1)=-h(j1min)*f(j1min)/(xj1*e(sj1p1))
60   g=h(j1min+1)
n2=j1max-1
xj1=j1max
sj1=xj1*xj1
h(n2)=-h(j1max)*f(j1max)/((xj1+1.0d0)*e(sj1))
mid1=mid-1
IF(ncoef == 3) GOTO 80
!----------------------------------------------------------------------
!  start recursive evaluation of the coefficients using
!  equation 12'.
!----------------------------------------------------------------------
DO i=j1min+1, mid1
    ip1=i+1
    im1=i-1
    xj1=i
    xj1p1=xj1+1.0d0
    sj1=xj1*xj1
    sj1p1=xj1p1*xj1p1
    g=(-f(i)*h(i)-xj1p1*e(sj1)*h(im1))/(xj1*e(sj1p1))
    n2=j1max+j1min-i
    IF(n2 == mid) GOTO 80
    h(ip1)=g
    xj1=n2
    xj1p1=xj1+1.0d0
    ip1=n2+1
    im1=n2-1
    sj1=xj1*xj1
    sj1p1=xj1p1*xj1p1
    h(im1)=(-xj1*e(sj1p1)*h(ip1)-f(n2)*h(n2))/(xj1p1*e(sj1))
ENDDO
!----------------------------------------------------------------------
!  match at the midpoint.
!----------------------------------------------------------------------
80   IF(ABS(g) < 1.d-8)THEN
    i=mid
    im1=mid1
    xj1=mid
    xj1p1=xj1+1.0d0
    sj1=xj1*xj1
    sj1p1=xj1p1*xj1p1
    g=(-f(i)*h(i)-xj1p1*e(sj1)*h(im1))/(xj1*e(sj1p1))
    mid1=mid
    mid=mid+1
ENDIF
coef=h(mid)/g
DO i=j1min, mid1
    h(i)=coef*h(i)
ENDDO
!----------------------------------------------------------------------
!  find magnitude of the normalization coefficient.
!----------------------------------------------------------------------
sum=0.0d0
x=j1min*2-1
DO i=j1min,j1max
    x=x+2.0d0
    sum=sum+x*h(i)*h(i)
ENDDO
y=2*l1+1
sum=sum*y
coef=sqrt(1.0d0/sum)
!----------------------------------------------------------------------
!  determine sign of normalization coefficient
!----------------------------------------------------------------------
j=j2+j3+l2+l3
IF(j /= 2*(j/2)) coef=-coef
!----------------------------------------------------------------------
!  normalize to obtain proper absolute value
!----------------------------------------------------------------------
DO i=j1min, j1max
    h(i)=coef*h(i)
ENDDO
120  c6j=h(j1)
RETURN
130  IF(ncoef == 2) GOTO 140
!----------------------------------------------------------------------
!  special case for only 1 coefficient.
!----------------------------------------------------------------------
x=(2*j1min+1)*(2*l1+1)
h(j1min)=sqrt(1.d0/x)
j=j2+j3+l2+l3
IF((j/2)*2 /= j) h(j1min)=-h(j1min)
GOTO 120
!----------------------------------------------------------------------
!  special case for only 2 coefficients.
!----------------------------------------------------------------------
140  h(j1max)=1.d0
xj1=j1max
sj1=xj1*xj1
h(j1min)=-f(j1max)/((xj1+1)*e(sj1))
x=j1min*2+1
sum=x*h(j1min)*h(j1min)+(x+2.d0)
sum=sum*(2*l1+1)
coef=sqrt(1.d0/sum)
j=j2+j3+l2+l3
IF(j /= 2*(j/2)) coef=-coef
h(j1min)=coef*h(j1min)
h(j1max)=coef*h(j1max)
GOTO 120
!----------------***end-f6j***-----------------------------------------
CONTAINS
   REAL(Kind=WP_Kind) FUNCTION E(sj1)
      USE Numeric_Kinds_Module
      REAL(Kind=WP_Kind) sj1
      e=sqrt((sj1-sj2mj3)*(sj2j31-sj1)*(sj1-sl2ml3)*(sl2l31-sj1))
   END FUNCTION E
   REAL(Kind=WP_Kind) FUNCTION F(j1)
      USE Numeric_Kinds_Module
      INTEGER j1
      f=(2.0d0*j1+1)*((-j1*(j1+1)+c1)*j1*(j1+1)+c0)
   END FUNCTION F
END FUNCTION C6J
