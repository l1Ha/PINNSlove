SUBROUTINE savrms(phirms, ntheta, nchi, megamin, megamax)
USE Numeric_Kinds_Module
USE FileUnits_Module
!
!   $RCSfile: savrms.f,v $ $Revision: 1.15 $
!   $Date: 89/11/16 11:03:28 $
!   $State: Stable $
!

!            P U R P O S E   O F    S U B R O U T I N E


!            I N P U T    A R G U M E N T S
!    phirms
!    ntheta
!    nchi
!    megamin
!    megamax


IMPLICIT NONE
INTEGER ntheta, nchi, megamin, megamax, imega
REAL(Kind=WP_Kind) phirms(ntheta, nchi)

WRITE(Out_unit,*)'savrms: ntheta, nchi = ',ntheta,nchi

REWIND PhiRms1_Unit
REWIND PhiRms2_Unit

DO imega = megamin, megamax
   READ(PhiRms1_Unit)phirms
   WRITE(PhiRms2_Unit) phirms
ENDDO

REWIND PhiRms1_Unit
REWIND PhiRms2_Unit

RETURN
END
