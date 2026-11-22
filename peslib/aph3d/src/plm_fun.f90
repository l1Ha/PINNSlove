REAL(Kind=WP_Kind) FUNCTION plm_fun(l, mm, x)
USE Numeric_Kinds_Module
USE fileunits_Module
!  
IMPLICIT NONE
INTEGER l, mm, m, mtest, i, lm1, mu
REAL(Kind=WP_Kind) x, xm, p1, p2, p3, z, xtest, xi, rat, xl
REAL(Kind=WP_Kind) xmu, xnorm, xlm, xlp
CHARACTER(LEN=21), PARAMETER:: ProcName='plm_fun'

!  computes normalized associated legendre polynomials by recursion.
!  the phase is that of abramowitz and stegun and gradshteyn and ryzhik.
!  recurs up in l via gr 8.832.2 then up in abs(mm) via 8.735.1 and 8.73
!  for large abs(x), recurs down in abs(mm) via gr 8.733.1d and 8.733.3
!  then each plm is normalized to 1.0 on (-1.,1.).
!  for either positive or negative mm, multiplication by exp(i*mm*phi)/
!  sqrt(2.*pi) gives the spherical harmonics of condon and shortley, ros
!  messiah, and jackson.
!  to get the phase of mclean and yoshimine, replace statement 333 by 44
!  program of r. nerf modified by s. green and r. pack.
m=iabs(mm)
xl=float(l)
xm=float(m)
IF(m > l .or. l < 0) go to 100
IF(abs(x) > 1.d0) go to 100
p1=1.d0
p2=x
IF(l == 0) go to 70
IF(m == 0) go to 10
!  if x is small enough, use the recursion expected to be faster.
xtest=0.499d0*(1.d0+1.d0/xm)
mtest=l/3
IF(m > mtest .or. abs(x) > xtest) go to 40
10   DO 20 i=1, l
   xi=float(i)
   p3=((2.d0*xi+1.d0)*x*p2-xi*p1)/(xi+1.d0)
   p1=p2
   p2=p3
20   continue
IF(m == 0) go to 70
!  at end of loop p1=p(l, 0, x)
z=sqrt(1.d0-x*x)
IF(z < 1.d-6) go to 110
!  if z=0, then x=1 and plm(1.)=0 for m > 0.
!  recur up in abs(mm)
p2=(xl+1.d0)*(p2-x*p1)/z
DO 30 i=1,m
   xi=float(i)
   p3=-2.d0*x/z*p2*xi-(xl+xi)*(xl-xi+1.d0)*p1
   p1=p2
   p2=p3
30   continue
go to 70
!  recur down in abs(mm) if abs(x) is large
40   z=sqrt(1.d0-x*x)
IF(z < 1.d-6) go to 110
!  if z=0 then x=1 and plm=0 for m>0
!  calc ratio of factorials for pll
rat=1.d0
DO 50 i=1, l
   xi=float(i)
   rat=0.5d0*(xl+xi)*rat
50   continue
!  calculate pll
p1=rat*(-z)**l
IF(m == l) go to 70
p2=p1
p1=-x*p2/z
lm1=l-1
IF(m == lm1) go to 70
lm1=l-m-1
!  recur downward in abs(mm)
DO 60 i=1, lm1
   mu=l-i-1
   xmu=float(mu)
   p3=p2
   p2=p1
   p1=2.d0*(xmu+1.d0)*x*p2/z+p3
   p1=-p1/((xl-xmu)*(xl+xmu+1.d0))
60   continue
!  normalization.....
70   xnorm=(2.d0*xl+1.d0)/2.d0
IF(m == 0) go to 90
xlm=xl+1.d0
xlp=xl
DO 80 i=1, m
   xlm=xlm-1.d0
   xlp=xlp+1.d0
   xnorm=xnorm/(xlm*xlp)
80   continue
90   plm_fun=p1*sqrt(xnorm)
IF(mm < 0) plm_fun=plm_fun*(-1)**(m)
!  444 IF(m /= 0) plm=plm*(-1)**m
RETURN
100  WRITE(Out_Unit,120) l, mm, x
110  plm_fun=0.d0
RETURN
!  
120  FORMAT(//,'       error.  argument out of range for plm(', 2i6 &
      , e14.7, '  ).')
END
