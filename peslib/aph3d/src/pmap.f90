SUBROUTINE pmap (kase, iproj, chij, theta, chi, u, v)
USE FileUnits_Module
USE Numeric_Kinds_Module
USE Numbers_Module
INTEGER IPROJ, Kase
REAL(KIND=WP_Kind) a, ca, chi, chij, crot, r, rot, theta, u, up, v, vp
REAL(KIND=WP_Kind) x, y, z
!-----------------------------------------------------------------------
!...kase=1  thetap=0 and therefore 
!              a=theta
!              rot=chi
!...kase=2  thetap=pi/2 and therefore
!              cos(a)=cos(chij)sin(theta)cos(chi)+sin(chij)sin(theta)sin(chi)
!              cos(rot)=up/sqrt(up**2+vp**2)
!              where
!                 up=z
!                 vp=-xsin(chij)+y*cos(chij)
!...we have the following:
!...u=r*sin(a)
!...v=r*cos(a)
!...where (see ncar manual for details)
!...iproj=1 stereographic      r=tan(a/2)            0<=a<pi
!...iproj=2 orthographic       r=sin(a)              0<=a<=pi
!...iproj=3 lambert equal area r=2sin(a)/sqrt(2(1+cos(a))) 0<=a<pi
!...iproj=4 gnomonic           r=tan(a)                   0<=a<pi/2
!...iproj=5 equidistance       r=a                         0<=a<=pi
!-----------------------------------------------------------------------
IF(kase==1)THEN
   a = theta
   rot = chi
ELSE
   ca = cos (chi - chij) * sin (theta)
   a = acos (ca)
   x = sin (theta) * cos (chi)
   y = sin (theta) * sin (chi)
   z = cos (theta)
   up = y
   vp = x * sin (chij) + z * cos (chij)
   crot = 0.0
   IF(up/=0.0) crot = up / sqrt (up**2 + vp**2)
   IF(crot>1.0) crot = 1.0
   IF(crot< - 1.0) crot = - 1.0
   rot = acos (crot)
ENDIF

IF(iproj==1)THEN
   r = tan (a / 2)
ELSEIF (iproj==2)THEN
   r = sin (a)
ELSEIF (iproj==3)THEN
   r = 2.0 * sin (a) / sqrt (2 * (1 + cos (a) ) )
ELSEIF (iproj==4)THEN
   IF(a>=halfpi)THEN
      WRITE(Out_Unit, * ) 'error angle a is too large: a=', a
      STOP 'pmap'
   ENDIF
   r = tan (a)
ELSEIF (iproj==5)THEN
   r = a
ELSEIF (iproj>=6)THEN
   WRITE(Out_Unit, * ) 'error iproj is too large: iproj=', iproj
   STOP 'pmap'
ENDIF

u = r * sin (rot)
v = r * cos (rot)
RETURN
ENDSUBROUTINE pmap
