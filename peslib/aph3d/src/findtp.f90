SUBROUTINE findtp (n, v, x, e, ntp, xtp)
USE Numeric_Kinds_Module
USE Numbers_Module
USE FileUnits_Module
!
!-------------------------------- findtp -------------------------------
!
! >>> findtp   -- finds the turning points of a potential at energy e
!
! *** calling sequence:   CALL findtp ( n, v, e, ntp, xtp )
!
! === latest revision:   | 900605 at 22:01:08 |  original typed in
!
! --- on entry the calling routine supplies --
!
!     n     -- the number of grid points on which potential is tabulated
!     v     -- array containing the tabulated potential
!     x     -- array containing the values of x on the grid
!     e     -- energy at which to locate turning points
!
! --- on exit the routine returns --
!
!     ntp   -- number of turning points found
!     xtp   -- array holding x values at turning points

INTEGER :: n, nmax
PARAMETER (nmax = 8192)
REAL(Kind=WP_Kind) ::       v (n)
REAL(Kind=WP_Kind) ::       x (n)
REAL(Kind=WP_Kind) ::       e
!
INTEGER ::        ntp
REAL(Kind=WP_Kind) ::       xtp (204)
!
REAL(Kind=WP_Kind) ::       aa (nmax), bb (nmax), cc (nmax), wt (nmax)
REAL(Kind=WP_Kind) ::       eold, enew
REAL(Kind=WP_Kind) ::       xz, xzn, cor
REAL(Kind=WP_Kind) ::       tab (3)
!
INTEGER ::        iop (2)
INTEGER ::        i, it
!
LOGICAL ::        last
!
!---------------------------- execution begins here --------------------
!
!
IF(n>=nmax)THEN
   WRITE(Msg_Unit, * ) 'Error: n>=nmax', n, nmax

ENDIF
!
iop (1)      = 5
iop (2)      = 5
!
CALL spl1d1 (n, x, v, wt (1), iop, 1, aa (1), bb (1), cc (1) )
!
!        --- loop through potential until we find a crossing ---
!
ntp    = 0
eold   = v (1) - e
DO 1 i = 2, n
   enew   = v (i) - e
   IF(eold * enew<Zero)THEN
      it     = 0
      last   = .false.
      xz     = x (i)
!
   11       CONTINUE
      it     = it + 1
      CALL spl1d2 (n, x, v, wt (1), 1, xz, tab)
      cor    = (tab (1) - e) / tab (2)
      xzn    = xz - cor
!
      IF(last) goto 12
      last = 1d0 + .5d0 * abs (cor / xzn) ==1d0
      IF(it>30) goto 12
!
      xz     = xzn
      goto 11
!
   12       CONTINUE
      ntp    = ntp + 1
      xtp (ntp) = xzn
   ENDIF
!
   eold   = enew
    1 ENDDO
!
RETURN
ENDSUBROUTINE findtp
