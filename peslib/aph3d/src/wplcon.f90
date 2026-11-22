SUBROUTINE wplcon(z,x,y,nx,ny,nxd,cval, npmax,ncmax, xcon,ycon,ntot, icrv,npts,ncrv, ierr )
USE Numeric_Kinds_Module
USE FileUnits_Module
INTEGER	andi, ifake, i
INTEGER icrv(ncmax), npts(ncmax)
REAL(Kind=WP_Kind) z(nxd,1), x(nxd,1), y(nxd,1)
REAL(Kind=WP_Kind) xcon(1), ycon(1), cval
REAL(Kind=WP_Kind)	andx, orx, xfake
! PORT NOTE (2026-10-01): statement-function dummies renamed (ii,xx) so they no
! longer collide with the host array names i/x; result types and the VMS
! intrinsics and/or/not are unchanged.
INTEGER ii
REAL(Kind=WP_Kind) xx
andi(ii,xx) = and(ii,ifake(xx))
andx(xx,ii) = xfake(and(ifake(xx),ii))
orx (xx,ii) = xfake(or (ifake(xx),ii))
mskw = 1
mskw = not(mskw)
msk1 = 1
xmx = x(1,1)
xmn = xmx
ymx = y(1,1)
ymn = ymx
DO i = 1,nx
  xmx = max(xmx,x(i,ny))
  xmn = min(xmn,x(i,1))
ENDDO
DO j = 1, ny
  ymx = max(ymx,y(1,j))
	ymn = min(ymn,y(nx,j))
ENDDO
rngx = abs(xmx - xmn)
rngy = abs(ymx - ymn)
ncrv = 0
ntot = 0
ierr = 0
icrv(1) = 1
npts(1) = 0
nx1 = nx-1
ny1 = ny-1
k = 0
xs = -1000.
ys = -1000.
DO j = 1,ny
  DO i = 1,nx
  x(i,j) = andx(x(i,j),mskw)
  y(i,j) = andx(y(i,j),mskw)
  ENDDO
ENDDO
ia = 1
i = 1
  190 IF ((cval-z(i,1))*(cval-z(i+1,1)).gt.0.) go to 200
ic = i
jc = 1
ib = -1
in = 1
go to 725
  200 i = i+1
i1 = 0
IF (i.le.nx1) go to 190
ia = 2
i = 1
  220 IF ((cval-z(i,ny))*(cval-z(i+1,ny)).gt.0.) go to 240
ic = i
jc = ny-1
in = 3
ib = -1
go to 745
  240 i = i+1
i1 = 0
IF (i.le.nx1) go to 220
ia = 3
j = 1
  250 IF ((cval-z(1,j))*(cval-z(1,j+1)).gt.0.) go to 260
ic = 1
jc = j
ib = -1
in = 4
go to 755
  260 j = j+1
i1 = 0
IF (j.le.ny1) go to 250
ia = 4
j = 1
  270 IF ((cval-z(nx,j))*(cval-z(nx,j+1)).gt.0.) go to 280
ic = nx-1
jc = j
ib = -1
in = 2
go to 735
  280 j = j+1
i1 = 0
IF (j.le.ny1) go to 270
ia = 5
j = 2
  380 i = 1
  390 IF ((cval-z(i,j))*(cval-z(i+1,j)).gt.0.) go to 400
ic = i
jc = j
ib = -1
in = 1
go to 725
  400 i = i+1
i1 = 0
IF (i.le.nx1) go to 390
j = j+1
IF (j.le.ny1) go to 380
go to 2000
  500 c = 0.
il = 10
IF (z(i2,j2).ne.z(i1,j1)) c = (cval-z(i1,j1))/(z(i2,j2)-z(i1,j1))
x1 = c*(x(i2,j2)-x(i1,j1))+x(i1,j1)
y1 = c*(y(i2,j2)-y(i1,j1))+y(i1,j1)
  550 ip = 1
IF (z(i2,j2).ge.z(i1,j1)) go to 630
  590 IF ((cval-z(i1,j1))*(cval-z(i4,j4)).gt.0.) go to (600,700), ip
id = 1
c  = 1.
IF (z(i1,j1).ne.z(i4,j4)) c = (cval-z(i1,j1))/(z(i4,j4)-z(i1,j1))
x2 = c*(x(i4,j4)-x(i1,j1))+x(i1,j1)
y2 = c*(y(i4,j4)-y(i1,j1))+y(i1,j1)
go to 700
  600 IF ((cval-z(i4,j4))*(cval-z(i3,j3)).gt.0.) go to (640,590), ip
id = 0
c = (cval-z(i4,j4))/(z(i3,j3)-z(i4,j4))
x2 = c*(x(i3,j3)-x(i4,j4))+x(i4,j4)
y2 = c*(y(i3,j3)-y(i4,j4))+y(i4,j4)
go to 700
  630 ip = 2
  640 IF ((cval-z(i3,j3))*(cval-z(i2,j2)).gt.0.) go to (700,600), ip
  650 id = 3
c = 1.
IF (z(i2,j2).ne.z(i3,j3)) c = (cval-z(i2,j2))/(z(i3,j3)-z(i2,j2))
x2 = c*(x(i3,j3)-x(i2,j2))+x(i2,j2)
y2 = c*(y(i3,j3)-y(i2,j2))+y(i2,j2)
  700 IF (abs(x2-x1)/rngx + abs(y2-y1)/rngy .lt. 1e-7) go to 710
xb = x2
xa = x1
yb = y2
ya = y1
 705  CONTINUE
IF ( abs(xa-xs)/rngx + abs(ya-ys)/rngy .lt. 1e-4 ) go to 706
ncrv = ncrv + 1
ntot = ntot + 1
IF ( ntot .gt. npmax ) ierr = ierr + 1
IF ( ncrv .gt. ncmax ) ierr = ierr + 2
IF ( ierr .ne. 0 ) go to 2001
icrv(ncrv) = ntot
xcon(ntot) = xa
ycon(ntot) = ya
 706  CONTINUE
xs = xb
ys = yb
ntot = ntot + 1
IF ( ntot .gt. npmax ) ierr = ierr + 1
IF ( ierr .ne. 0 ) go to 2001
xcon(ntot) = xb
ycon(ntot) = yb
  707 x1 = x2
y1 = y2
  710 in = mod(in+id-1,4)+1
go to (720,730,740,750), in
  720 jc = jc+1
IF (jc.lt.ny) go to 725
x(ic,ny) = orx(x(ic,ny),msk1)
go to (200,240,260,280,400), ia
  725 CONTINUE
ma = andi(msk1,x(ic,jc))
IF (ma.eq.1) go to (200,240,260,280,400), ia
j1 = jc
i1 = ic
i2 = i1+1
j2 = j1
i3 = i2
j3 = j1+1
j4 = j3
i4 = i1
x(ic,jc) = orx(x(ic,jc),msk1)
ib = ib+1
IF (ib) 500,500,550
  730 ic = ic-1
IF (ic.gt.0) go to 735
y(1,jc) = orx(y(1,jc),msk1)
go to (200,240,260,280,400), ia
  735 CONTINUE
ma = andi(msk1,y(ic+1,jc))
IF (ma.eq.1) go to (200,240,260,280,400), ia
i1 = ic+1
j1 = jc
i2 = i1
j2 = j1+1
i3 = ic
j3 = j2
i4 = i3
j4 = jc
y(i1,j1) = orx(y(i1,j1),msk1)
ib = ib+1
IF (ib) 500,500,550
  740 jc = jc-1
IF (jc.gt.0) go to 745
x(ic,1) = orx(x(ic,1),msk1)
go to (200,240,260,280,400), ia
  745 CONTINUE
ma = andi(msk1,x(ic,jc+1))
IF (ma.eq.1) go to (200,240,260,280,400), ia
j1 = jc+1
i1 = ic+1
i2 = ic
j2 = j1
i3 = i2
j3 = jc
i4 = i1
j4 = j3
x(i2,j2) = orx(x(i2,j2),msk1)
ib = ib+1
IF (ib) 500,500,550
  750 ic = ic+1
IF (ic.lt.nx) go to 755
y(nx,jc) = orx(y(nx,jc),msk1)
go to (200,240,260,280,400), ia
  755 CONTINUE
ma = andi(msk1,y(ic,jc))
IF (ma.eq.1) go to (200,240,260,280,400), ia
i1 = ic
j1 = jc+1
i2 = i1
j2 = jc
i3 = ic+1
j3 = j2
i4 = i3
j4 = j1
y(i2,j2) = orx(y(i2,j2),msk1)
ib = ib+1
IF (ib) 500,500,550
 2000 CONTINUE
 2001 CONTINUE
IF ( ierr .eq. 0 ) go to 2099
ncrv = min0(ncrv,ncmax)
ntot = min0(ntot,npmax)
 2099 CONTINUE
IF ( ncrv .eq. 0 ) RETURN
ncm1 = ncrv - 1
IF ( ncm1 .eq. 0 ) go to 2101
DO i = 1, ncm1
  npts(i) = icrv(i+1)-icrv(i)
ENDDO
 2101 CONTINUE
npts(ncrv) = ntot + 1 - icrv(ncrv)
RETURN
END
