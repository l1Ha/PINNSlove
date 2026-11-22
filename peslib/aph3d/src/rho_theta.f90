SUBROUTINE rho_theta (nrho, ntheta, nchi, rho_val, theta_val, func, filename, vert, poly)
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
CHARACTER(LEN=*) :: filename
INTEGER :: nrho, nchi, ntheta, irho, itheta, inode, irho_frst, nvertex, npolygon, node, surf, counter
REAL(Kind=WP_Kind) :: rho_val (nrho), theta_val (ntheta), func (nrho, ntheta, nchi), rho, theta
REAL(KIND=WP_Kind) :: x, y, z, vert (3, nrho * ntheta)

INTEGER :: poly (5 * nrho * ntheta)

OPEN(Unit=Rho_Theta_Plt_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//filename)
IF(nchi/=1)THEN
   WRITE(Out_Unit, * ) 'Error: nchi/=1, nchi=', nchi
   STOP 'rho_theta'
ENDIF

IF(rho_val (1) ==0.d0)THEN
   irho_frst = 2
ELSE
   irho_frst = 1
ENDIF

node = 0
nvertex = 0
IF(irho_frst==2)THEN
   x = 0.0
   y = 0.0
   z = func (1, 1, 1)
   nvertex = 1
   node = node+1
   vert (3, node) = x
   vert (2, node) = y
   vert (3, node) = z
   WRITE(Rho_Theta_Plt_Unit, 70) x, y, z
ENDIF

DO irho = irho_frst, nrho
DO itheta = 1, ntheta
rho = rho_val (irho)
theta = theta_val (itheta)
x = rho * cos (theta)
y = rho * sin (theta)
z = func (irho, itheta, nchi)
WRITE(Rho_Theta_Plt_Unit, 70) x, y, z
node = node+1
vert (1, node) = x
vert (2, node) = y
vert (3, node) = z
ENDDO
ENDDO

nvertex = (nrho - irho_frst + 1) * ntheta + nvertex
IF(nvertex/=node)THEN
   WRITE(Out_Unit, * ) 'Error with node'
   STOP 'node'
ENDIF

counter = 0
surf = 0

npolygon = 0
IF(irho_frst==1)THEN
   inode = 0
ELSE
   npolygon = ntheta
   inode = 1
   DO itheta = 1, ntheta
   surf = surf + 1
   IF(itheta/=ntheta)THEN
      poly (surf + counter) = 3
      poly (surf + 1 + counter) = inode-1
      poly (surf + 2 + counter) = itheta
      poly (surf + 3 + counter) = itheta + 1
      WRITE(Rho_Theta_Plt_Unit, 80) inode, itheta + 1, itheta + 2
   ELSE
      poly (surf + counter) = 3
      poly (surf + 1 + counter) = inode-1
      poly (surf + 2 + counter) = ntheta
      poly (surf + 3 + counter) = inode
      WRITE(Rho_Theta_Plt_Unit, 80) inode, ntheta + 1, inode+1
   ENDIF
   counter = counter + 3
   ENDDO
ENDIF

DO irho = irho_frst, nrho - 1
DO itheta = 1, ntheta
surf = surf + 1
inode = inode+1
IF(itheta/=ntheta)THEN
   poly (surf + counter) = 4
   poly (surf + 1 + counter) = inode-1
   poly (surf + 2 + counter) = inode+ntheta - 1
   poly (surf + 3 + counter) = inode+ntheta
   poly (surf + 4 + counter) = inode
   WRITE(Rho_Theta_Plt_Unit, 80) inode, inode+1, inode+ntheta, inode+ntheta + 1
ELSE
   poly (surf + counter) = 4
   poly (surf + 1 + counter) = inode-1
   poly (surf + 2 + counter) = inode+ntheta - 1
   poly (surf + 3 + counter) = inode
   poly (surf + 4 + counter) = inode-ntheta
   WRITE(Rho_Theta_Plt_Unit, 80) inode, inode-ntheta + 1, inode+1, inode+ntheta
ENDIF
counter = counter + 4
ENDDO
ENDDO

npolygon = npolygon + (nrho - irho_frst) * ntheta
IF(npolygon/=surf)THEN
   WRITE(Out_Unit, * ) 'Error at poly'
   STOP 'poly'
ENDIF

CALL lwo (nvertex, vert, npolygon, poly, filename)

CLOSE(unit=Rho_Theta_Plt_Unit)
RETURN
   70 FORMAT('v',3f12.7)
   80 FORMAT('f',4I8)
ENDSUBROUTINE rho_theta
