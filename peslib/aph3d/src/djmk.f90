REAL(Kind=WP_Kind) FUNCTION djmk (j, mm, kk, z, a, n)
USE NUmeric_Kinds_Module
USE FileUnits_Module
USE factrl_Module
!
!  djmk -- computes a sequence of wigner rotation functions by
!  recursion on the angular momentum index, holding the
!  projection quantum numbers fixed
!
!  --- on input the calling routine supplies --
!
!  j  -- the final desired angular momentum quantum number
!  mm -- the (space-fixed) projection quantum number
!  kk -- the (body-fixed) projection quantum number
!  z  -- the angular argument, in radians
!  n  -- the maximum length allowed for the output array a (below)
!
!  --- on exit the routine returns --
!
!  a  -- an array of length l, containing the wigner d-functions
!  where a(1)=d(jstar,mm,kk), a(2)=d(jstar+1,mm,kk), ...,
!  and a(l)=d(j,mm,kk), and l=j-jstar+1, where
!  jstar = max(ABS(mm),ABS(kk)).  all of this assumes that
!  l is less than the imposed maximum length of the
!  array a (n, the input variable).  IF l would be greater
!  than n, THEN the recursion stops before the array a
!  would be overwritten, and a(n)=d(n+jstar-1,mm,kk)
!  djmk -- the function value returned is a(l), which is usually
!  d(j,mm,kk), unless insufficient space has been provided
!  in the output array a to CONTINUE the recursion to the
!  requested final j.
!
!  --- note --
!
!  this routine also expects the common block factrl to be fill
!  with the logarithms of the factorial function.  to DO this,
!  precede the first CALL to djmk with a CALL to the routine
!  factor, which will properly initialize the common block
!  a CALL to numgen should precede the CALL to factor.
!
IMPLICIT NONE
INTEGER kk, mm, ma, iabs, ka, k, m, jj, nmax, j, jp, jpk, i, jmk, n
LOGICAL mgtk, kpos, mpos
REAL(KIND=WP_Kind) a(n), fkm, z, zz, czz, szz, ajs, ams, b, aks, du
REAL(KIND=WP_Kind) aj, cz, f, d, df, amk, dd, dn, dof, ajps
!
fkm=(-1)**(kk-mm)
ma=iabs(mm)
ka=iabs(kk)
mgtk=ma > ka
mpos=mm > 0
kpos=kk > 0
!
IF(mgtk) GOTO 20
IF(kpos) GOTO 10
k=-mm
m=-kk
f=1.d0
GOTO 40
 10   m=kk
k=mm
f=fkm
GOTO 40
!
 20   IF(mpos) GOTO 30
m=-mm
k=-kk
f=fkm
GOTO 40
 30   m=mm
k=kk
f=1.d0
!
 40   jj=m
nmax=j-m+1
IF(nmax > n) print 60, nmax, j, mm, kk, n
zz=z*.5d0
czz=cos(zz)
szz=sin(zz)
cz=czz*czz-szz*szz
jp=jj+1
jpk=jp+k
jmk=jp-k
du=(-1)**(jj-k)*exp(.5d0*(ftl(jj+jp)-ftl(jpk)-ftl(jmk)))*f
IF( jpk > 1 ) du = du*czz**(jpk-1)
IF( jmk > 1 ) du = du*szz**(jmk-1)
a(1)=du
djmk=du
IF(nmax == 1) RETURN
!
d=(jp*cz-k)*sqrt(REAL(jj+jp,WP_Kind)/(jpk*jmk))*du
a(2)=d
djmk=d
IF(nmax == 2) RETURN
!
aks=k*k
ams=m*m
amk=m*k
!
DO 50 i=3, nmax
   jp=m+i-1
   jj=jp-1
   aj=jj
   ajps=jp*jp
   ajs=aj*aj
   dd=sqrt((ajps-ams)*(ajps-aks))
   b=sqrt((ajs-ams)*(ajs-aks))
   df=(cz*d-b*du/ajs)*aj*jp/dd
   dof=(cz*ajps/dd-amk*(jj+jp)/(dd*aj))*d
   dn=df+dof
   a(i)=dn
   du=d
   d=dn
 50   CONTINUE
!
djmk=a(nmax)
RETURN
!
 60   FORMAT(1x, "djmk  -1", 3x, 'overWRITEs array', 5i6)
END
