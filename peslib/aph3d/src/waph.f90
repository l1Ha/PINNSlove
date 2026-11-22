SUBROUTINE waph(f2, naph, numnp)
USE Numeric_Kinds_Module
USE FileUnits_Module
!
!   $RCSfile: waph.f,v $ $Revision: 1.15 $
!   $Date: 89/11/16 11:03:35 $
!   $State: Stable $
!

!            P U R P O S E   O F    S U B R O U T I N E


!            I N P U T    A R G U M E N T S
!    f2
!    naph
!    numnp


IMPLICIT NONE

!            L O G I C A L S


!            I N T E G E R S

INTEGER naph, numnp, istart, iend, iaph, inode

!            R E A L S
REAL(Kind=WP_Kind) f2

!            D I M E N S I O N S
DIMENSION f2(numnp, naph)


!            C O M M O N S

!            I N T R I N S I C    F U N C T I O N S
!            E X T E R N A L S

DATA istart/0/
!-----------------------------------------------------------------------
! On the first CALL OPEN unit=FEM_SF_Bin_Unit.
!-----------------------------------------------------------------------
IF(istart==0)THEN
   OPEN(Unit=FEM_SF_Bin_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Sfunf2.bin', form='unformatted',  status='unknown')
ENDIF
!-----------------------------------------------------------------------
! WRITE surface functions in a FORMAT suitable for routine ovr.
!-----------------------------------------------------------------------
REWIND FEM_SF_Bin_Unit
DO 16 istart = 1, naph, 100
   iend = min(istart+99, naph)
   WRITE(FEM_SF_Bin_Unit)((f2(inode,iaph),inode=1,numnp),iaph=istart,iend)
   16 CONTINUE
REWIND FEM_SF_Bin_Unit
RETURN
END
