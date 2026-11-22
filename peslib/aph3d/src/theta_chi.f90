SUBROUTINE theta_chi(nrho,ntheta,nchi,theta_val,chi_val,func, filename, eig, ichanl, nchanl, vert, poly, vert_p, poly_p, fx, fy, f)
USE Numeric_Kinds_Module
USE FileUnits_Module
USE InputFile_Module
IMPLICIT NONE
INTEGER, PARAMETER :: maxcrv=1001
INTEGER, PARAMETER :: maxcnt=5120
INTEGER mxgrid, ntot, ncrv, ierr
INTEGER mxcont, i, j, nstart, nend, ic, i1, npc, i2
INTEGER icrv(maxcrv), npts(maxcrv)
REAL(KIND=WP_Kind) xcon(maxcnt), ycon(maxcnt), seperation, cont_start
CHARACTER(LEN=*) filename
INTEGER nrho,nchi,ntheta,itheta,ichi,inode, nvertex, npolygon, ncont
INTEGER icont, iflag, node, surf, counter, nchanl, ichanl, node_p, inode_p, dimensions, ncontours
REAL(Kind=WP_Kind) theta_val(ntheta), chi_val(nchi), func(nrho,ntheta,nchi), theta, chi, eig(nchanl)
REAL(KIND=WP_Kind) x, y, z, fx(ntheta,nchi), fy(ntheta,nchi), f(ntheta,nchi), cval(100)
REAL(KIND=WP_Kind) fmax, fmin, vert(3,ntheta*nchi), vert_p(3,nchi+1), zmax, cut_off, pot_min, zscale
INTEGER poly(ntheta*nchi), poly_p(nchi+1)
NAMELIST/potpar/ cut_off

DATA iflag/0/

OPEN(Unit=Theta_Chi_Plt_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//filename)

IF(nrho/=1)THEN
   WRITE(Out_Unit,*)'Error: nrho/=1, nrho=',nrho
   STOP 'theta_chi'
ENDIF

IF(theta_val(1)/=0.d0)THEN
   WRITE(Out_Unit,*)'Error'
   STOP 'theta_chi'
ENDIF

zmax=0
DO itheta=2,ntheta
  DO ichi=1,nchi
    IF(ABS(func(nrho,itheta,ichi))>zmax)THEN
       zmax=ABS(func(nrho,itheta,ichi))
    ENDIF
  ENDDO
ENDDO

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile)
READ(In_Unit,NML=potpar)
CLOSE(Unit=In_Unit)
WRITE(Out_Unit,NML=potpar)

!      zscale=0.4*theta_val(ntheta)/cut_off
zscale=1

IF(filename(1:1)=='v')THEN
  pot_min=cut_off
  DO itheta=2,ntheta
    DO ichi=1,nchi
      IF(func(nrho,itheta,ichi)<pot_min)THEN
         pot_min=func(nrho,itheta,ichi)
      ENDIF
    ENDDO
  ENDDO
ENDIF

x=0.0
y=0.0
z=func(1,1,1)
IF(filename(1:1)=='p')THEN
  z=zscale*((1/zmax)*((eig(ichanl+1)-eig(ichanl))/2)*z+eig(ichanl))
ELSE
  z=zscale*z
ENDIF

IF(z>(zscale*cut_off))THEN
  z=zscale*cut_off
ENDIF
WRITE(Theta_Chi_Plt_Unit,70)x,y,z

node=1
vert(1,node)=x
vert(2,node)=y
vert(3,node)=z
node_p=0
DO itheta=2,ntheta
   DO ichi=1,nchi
      theta=theta_val(itheta)
      chi=chi_val(ichi)
      x=theta*cos(chi)
      y=theta*sin(chi)
      z=func(nrho,itheta,ichi)
      IF(filename(1:1)=='p')THEN
        z=zscale*((1/zmax)*((cut_off-pot_min)/2)*z+eig(ichanl))
      ELSE
        z=z*zscale
      ENDIF
      IF(z>(zscale*cut_off))THEN
        z=zscale*cut_off
      ENDIF
      WRITE(Theta_Chi_Plt_Unit,70)x,y,z
      node=node+1
      vert(1,node)=x
      vert(2,node)=y
      vert(3,node)=z
      IF(itheta==ntheta)THEN
        node_p=node_p+1
        vert_p(1,node_p)=x
        vert_p(2,node_p)=y
        vert_p(3,node_p)=zscale*eig(ichanl)
      ENDIF
   ENDDO
ENDDO

nvertex=(ntheta-1)*nchi+1

IF (nvertex/=node)THEN
   WRITE(Out_Unit,*)'Error with node'
   STOP 'node'
ENDIF

surf=1
counter=0
inode=1
DO ichi=1,nchi
   IF(ichi/=nchi)THEN
      poly(surf+counter)=3
      poly(surf+1+counter)=inode-1
      poly(surf+2+counter)=ichi
      poly(surf+3+counter)=ichi+1
      WRITE(Theta_Chi_Plt_Unit,80)inode,ichi+1,ichi+2
   ELSE
      poly(surf+counter)=3
      poly(surf+1+counter)=inode-1
      poly(surf+2+counter)=nchi
      poly(surf+3+counter)=inode
      WRITE(Theta_Chi_Plt_Unit,80)inode,nchi+1,inode+1
   ENDIF
   counter=counter+3
   surf=surf+1
ENDDO

DO itheta=2,ntheta-1
   DO ichi=1,nchi
      inode=inode+1
      IF(ichi/=nchi)THEN
         poly(surf+counter)=4
         poly(surf+1+counter)=inode-1
         poly(surf+2+counter)=inode+nchi-1
         poly(surf+3+counter)=inode+nchi
         poly(surf+4+counter)=inode
         WRITE(Theta_Chi_Plt_Unit,80)inode,inode+nchi,inode+nchi+1,inode+1
      ELSE
         poly(surf+counter)=4
         poly(surf+1+counter)=inode-1
         poly(surf+2+counter)=inode+nchi-1
         poly(surf+3+counter)=inode
         poly(surf+4+counter)=inode-nchi
         WRITE(Theta_Chi_Plt_Unit,80)inode,inode+nchi,inode+1,inode-nchi+1
      ENDIF
      counter=counter+4
      surf=surf+1
   ENDDO
ENDDO

npolygon=nchi+(ntheta-2)*nchi

IF(npolygon/=(surf-1))THEN
    WRITE(Out_Unit,*)'Error at poly'
    STOP 'poly'
ENDIF

poly_p(1)=node_p
DO inode_p=2,node_p+1
  poly_p(inode_p)=inode_p-2
ENDDO

CALL lwo (nvertex, vert, npolygon, poly, filename)

IF(filename(1:1)=='p')THEN
  CALL plate (node_p, vert_p, 1, poly_p, filename)
ENDIF

dimensions=3
ncontours=4
seperation=0.01
cont_start=0

IF(filename(1:1)=='v')THEN
  CALL contour (nvertex, vert, npolygon, poly, filename, dimensions, ncontours, seperation, cont_start)
ENDIF

IF(iflag==0)THEN
WRITE(Out_Unit,*)'usrsub input values'
WRITE(Out_Unit,*)ntheta,nchi
fmin=func(nrho,1,1)
fmax=func(nrho,1,1)
DO ichi=1,nchi
   DO itheta=1,ntheta
      theta=theta_val(itheta)
      chi=chi_val(ichi)
      x=theta*cos(chi)
      y=theta*sin(chi)
      z=zscale*func(nrho,itheta,ichi)
      fx(itheta,ichi)=theta
      fy(itheta,ichi)=chi
      f(itheta,ichi)=z
      IF(z>fmax)fmax=z
      IF(z<fmin)fmin=z
   ENDDO
ENDDO

DO itheta=1,ntheta
   fx(itheta,nchi+1)=fx(itheta,1)
   fy(itheta,nchi+1)=8.d0*atan(1.d0)
   f(itheta,nchi+1)=f(itheta,1)
ENDDO


ncont=3
DO icont=1,ncont
  cval(icont)=fmin+icont*(fmax-fmin)/(ncont+1)
ENDDO
cval(1)=0.02
cval(2)=0.03
cval(3)=0.04
WRITE(Out_Unit,*)'ncont, fmin, fmax', ncont, fmin, fmax
WRITE(Out_Unit,*)(cval(icont),icont=1,ncont)
mxgrid=ntheta
mxcont=maxcnt
theta=theta_val(ntheta)
!ccc      DO ichi=1,nchi
!ccc         x=theta*sin(chi_val(ichi))
!ccc         y=theta*cos(chi_val(ichi))
!ccc         WRITE(Theta_Chi_ConGF_Unit,383)x,y
!ccc      ENDDO
!ccc      x=theta*sin(chi_val(1))
!ccc      y=theta*cos(chi_val(1))
!ccc      WRITE(Theta_Chi_ConGF_Unit,383)x,y
!ccc      WRITE(Theta_Chi_ConGF_Unit,*)'&'
DO icont = 1, ncont
   CALL wplcon (f, fx, fy, ntheta, nchi + 1, mxgrid, cval (icont), &
       mxcont, maxcrv, xcon, ycon, ntot, icrv, npts, ncrv, ierr)
   IF ( ncrv >= 1 )THEN
      WRITE(Theta_Chi_Con_Unit,*)'icont=',icont,'   contour value=',cval(icont)
      WRITE(Theta_Chi_Con_Unit,*)ntot,ncrv,(npts(i),icrv(i),i=1,ncrv)
      DO  i=1,ntot
         WRITE(Theta_Chi_Con_Unit,383)xcon(i),ycon(i)
      ENDDO
      nstart=1
      nend=npts(1)
      DO i=1,ncrv
         DO j=nstart,nend
            WRITE(Theta_Chi_ConGF_Unit,383)xcon(j),ycon(j)
         ENDDO
         WRITE(Theta_Chi_ConGF_Unit,*)'&'
         IF(i/=ncrv)THEN
            nstart=nstart+npts(i)
            nend=nend+npts(i+1)
         ENDIF
      ENDDO
      DO ic = 1, ncrv
         i1 = icrv(ic)
         npc = npts(ic)
         i2 = i1 + npc - 1
      ENDDO
   ENDIF
 ENDDO
iflag=100
ENDIF
npolygon=nchi+(ntheta-2)*nchi

CLOSE(Unit=Theta_Chi_Plt_Unit)

RETURN
   70 FORMAT('v',3f12.7)
   80 FORMAT('f',4I8)
383   FORMAT(1x,2e14.7)
ENDSUBROUTINE theta_chi
