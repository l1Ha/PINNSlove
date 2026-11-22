SUBROUTINE kbes (n,z,zjn,zjnp,zyn,zynp)
USE Numeric_Kinds_Module
!
!-------------------------------------------------------------------
! <kbes> -- computes ricatti bessel functions, reg + irreg.  + deriv
! l -- order of bessel function   (INTEGER(KIND=IW_Kind))
! x -- argument of bessel function    (REAL(Kind=WP_Kind))
! zjn -- returned value of regular bessel function
! zjnp -- returned value of deriv. of reg. bessel fn.
! zyn -- returned value of irregular bessel function
! zynp -- returned value of deriv. of irreg. bessel fn.
! x, zjn, zjnp, zjn, zjnp -- all REAL(Kind=WP_Kind)
!-------------------------------------------------------------------
USE FileUnits_Module
IMPLICIT NONE
INTEGER ithcll, ithsub, n, ito
REAL(Kind=WP_Kind) z,r, rp, zyn, zynp, zjn, zjnp
!#include <Inc.unitno>
LOGICAL little,medium,full
DATA ithcll/0/,ithsub/0/
DATA little/.false./,medium/.false./,full/.false./

CALL popt('kbes    ',little,medium,full,ithcll,ithsub)
!--------------------------------------------------------------------
! statement function ito
!--------------------------------------------------------------------
CALL rankel(n,  z, r, rp)
zyn=r
zynp=rp

CALL rankel(n, -z, r, rp)
zjn=0.5d0*(r-zyn)
!      WRITE(Msg_Unit,*)'zoe=',r,zyn
zjnp=0.5d0*(-rp-zynp)
IF(medium)THEN
    WRITE(Out_Unit,*)'n,z,zyn,zynp,zjn,zjnp,ito', n,z,zyn,zynp,zjn,zjnp,ito
ENDIF

RETURN
!----------------***end-kbes***---------------------------------------
END
