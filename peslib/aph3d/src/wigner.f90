SUBROUTINE wigner(nx, ny, nxdim, nydim, jrjtmax, betqf, djmkp, ic, jtot, lam, a, cnorm, parity)
!
!            P U R P O S E   O F    S U B R O U T I N E
! Calculates parity possessing Wigner d functions.

!            I N P U T    A R G U M E N T S
!  nx
!  ny
!  nxdim
!  nydmi
!  jr
!  jrjtmax
!  betqf
!  djmkp
!  ic
!  jtot
!  lam is the aph index.
!  a
!  cnorm
!  parity      Desired parity.

!            O U T P U T    A R G U M E N T S
!
!this routine is called by: AphDel_Old or AphDel_New
!this routine calls
!
!-----------------------------------------------------------------------
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

USE FileUnits_Module
USE Numeric_Kinds_Module

IMPLICIT NONE

REAL(Kind=WP_Kind) betqf, djmkp, beta, a, cnorm, djmk
INTEGER nx, ny, nxdim, nydim, jrjtmax, ic, jtot, k,lam, parity, iy, ix
EXTERNAL djmk
DIMENSION betqf(nxdim, nydim, 3), djmkp(0:jrjtmax, nxdim, nydim, 3), a(*)

!--------------------------------------------------------------------
!  now construct parity possessing Wigner small d functions
!  construct djlamk+(-1)**(jtot+lam+parity)*dj-lamk in Pack notation.
!  Bob Walker's function djmk uses Rose's notation, so we CALL
!  djklam to get Pack's djlamk
!-------------------------------------------------------------------
 WRITE(Out_Unit,*)'begining of routine wigner'
DO  iy=1,ny
   DO  ix=1,nx
      beta=betqf(ix,iy,ic)
      DO  k=0, jrjtmax
         djmkp(k,ix,iy,ic)=djmk(jtot,k,lam,beta,a,200)+(-1)**(jtot+lam+parity)*djmk(jtot,k,-lam,beta,a,200)
         djmkp(k,ix,iy,ic)=cnorm*djmkp(k,ix,iy,ic)
      ENDDO
   ENDDO
ENDDO
!djmkp=cnorm*djmkp
!
 WRITE(Out_Unit,*)'END of routine wigner'
RETURN
END
