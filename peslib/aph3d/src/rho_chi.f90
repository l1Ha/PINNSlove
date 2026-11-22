SUBROUTINE rho_chi (nrho, ntheta, nchi, rho_val, chi_val, func, flname, vert, poly)
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
CHARACTER(LEN=*) :: flname
INTEGER :: nrho, nchi, ntheta, irho, ichi, inode, irho_frst, nvertex, npolygon, node, surf, counter
REAL(Kind=WP_Kind) :: rho_val (nrho), chi_val (nchi), func (nrho, ntheta, nchi), rho, chi
REAL(KIND=WP_Kind) :: x, y, z, vert (3, nrho * nchi)

INTEGER :: poly (5 * nrho * nchi)

OPEN(Unit=Rho_Chi_Plt_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//flname)
IF(ntheta/=1)THEN
   WRITE(Out_Unit, * ) 'Error: ntheta/=1, ntheta=', ntheta
   STOP 'rho_chi'
ENDIF

IF(rho_val (1) ==0.d0)THEN
   irho_frst = 2
ELSE
   irho_frst = 1
ENDIF

nvertex = 0
node = 0
IF(irho_frst==2)THEN
   x = 0.0
   y = 0.0
   z = func (1, 1, 1)
   WRITE(Rho_Chi_Plt_Unit, 70) x, y, z
   nvertex = 1
   node = node+1
   vert (3, node) = x
   vert (2, node) = y
   IF(z> (0.2) )THEN
      z = (0.2)
   ENDIF
   vert (3, node) = z
ENDIF

DO irho = irho_frst, nrho
   DO ichi = 1, nchi
      rho = rho_val (irho)
      chi = chi_val (ichi)
      x = rho * cos (chi)
      y = rho * sin (chi)
      z = func (irho, ntheta, ichi)
      WRITE(Rho_Chi_Plt_Unit, 70) x, y, z
      node = node+1
      vert (1, node) = x
      vert (2, node) = y
      IF(z> (0.2) )THEN
         z = 0.2
      ENDIF
      vert (3, node) = z
   ENDDO
ENDDO

nvertex = (nrho - irho_frst + 1) * nchi + nvertex
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
   npolygon = nchi
   inode = 1
   DO ichi = 1, nchi
   surf = surf + 1
   IF(ichi/=nchi)THEN
      poly (surf + counter) = 3
      poly (surf + 1 + counter) = inode-1
      poly (surf + 2 + counter) = ichi
      poly (surf + 3 + counter) = ichi + 1
      WRITE(Rho_Chi_Plt_Unit, 80) inode, ichi + 1, ichi + 2
   ELSE
      poly (surf + counter) = 3
      poly (surf + 1 + counter) = inode-1
      poly (surf + 2 + counter) = nchi
      poly (surf + 3 + counter) = inode
      WRITE(Rho_Chi_Plt_Unit, 80) inode, nchi + 1, inode+1
   ENDIF
   counter = counter + 3
   ENDDO
ENDIF

DO irho = irho_frst, nrho - 1
   DO ichi = 1, nchi
      inode = inode+1
      surf = surf + 1
      IF(ichi/=nchi)THEN
         poly (surf + counter) = 4
         poly (surf + 1 + counter) = inode-1
         poly (surf + 2 + counter) = inode+nchi - 1
         poly (surf + 3 + counter) = inode+nchi
         poly (surf + 4 + counter) = inode
         WRITE(Rho_Chi_Plt_Unit, 80) inode, inode+1, inode+nchi, inode+nchi + 1
      ELSE
         poly (surf + counter) = 4
         poly (surf + 1 + counter) = inode-1
         poly (surf + 2 + counter) = inode+nchi - 1
         poly (surf + 3 + counter) = inode
         poly (surf + 4 + counter) = inode-nchi
         WRITE(Rho_Chi_Plt_Unit, 80) inode, inode-nchi + 1, inode+1, inode+nchi
      ENDIF
      counter = counter + 4
   ENDDO
ENDDO

npolygon = npolygon + (nrho - irho_frst) * nchi
IF(npolygon/=surf)THEN
   WRITE(Out_Unit, * ) 'Error at poly'
   STOP 'poly'
ENDIF

CLOSE(unit=Rho_Chi_Plt_Unit)

CALL lwo (nvertex, vert, npolygon, poly, flname)
RETURN
   70 FORMAT('v',3f12.7)
   80 FORMAT('f',4I8)
ENDSUBROUTINE rho_chi
