SUBROUTINE scaleit(nsize, a, njump, fac)
USE Numeric_Kinds_Module
!
!   $RCSfile: scale.f,v $ $Revision: 1.15 $
!   $Date: 89/11/16 11:03:29 $
!   $State: Stable $
!
!            P U R P O S E   O F    S U B R O U T I N E

!  This routine scales the array a by fac.
!   nsize   Lenght of the array a.
!   a       Array to be scaled.
!   njump   Scale every njump element.
!   fac     Scaling factor.


!            I N P U T    A R G U M E N T S
!    nsize
!    a
!    njump
!    fac

IMPLICIT NONE

!            L O G I C A L S


!            I N T E G E R S

INTEGER n, njump, nsize

!            R E A L S
REAL(Kind=WP_Kind) fac

!            D I M E N S I O N S

REAL(Kind=WP_Kind) a(nsize)

!            C O M M O N S

!            I N T R I N S I C    F U N C T I O N S
!            E X T E R N A L S


DO 1 n = 1, nsize, njump
   a(n) = a(n)*fac
    1 CONTINUE
RETURN
END
