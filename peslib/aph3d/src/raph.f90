!REAL(Kind=WP_Kind) function ranf(dummy)
!
!   $RCSfile: ranf.f,v $ $Revision: 1.15 $
!   $Date: 89/11/16 11:03:22 $
!   $State: Stable $
!
!
!    This function uses Alliant fortran SUBROUTINEs.
!
!    This version created 880126 by Jan Linderberg.
!
!    Generates random numbers with a psuedo-normal distribution.
!
!INTEGER dummy
!REAL(Kind=WP_Kind) avrge
!REAL(KIND=WP_Kind) rand
!LOGICAL first
!DATA first/.true./,avrge/1.5/
!SAVE first
!
!IF(first)THEN
!   ranf = rand(3)
!   first = .false.
!ENDIF
!
!ranf = rand(0) + rand(0) + rand(0) - avrge
!RETURN
!END
SUBROUTINE raph(f2, naph, numnp, istart)
USE Numeric_Kinds_Module
USE FileUnits_Module
!
!   $RCSfile: raph.f,v $ $Revision: 1.15 $
!   $Date: 89/11/16 11:03:22 $
!   $State: Stable $
!

!            P U R P O S E   O F    S U B R O U T I N E


!-----------------------------------------------------------------------
! This routine reads in 100 surface functions at a time to reduce memory
! storage.
!-----------------------------------------------------------------------

!            I N P U T    A R G U M E N T S
!    f2
!    naph
!    numnp
!    istart


IMPLICIT NONE

!            L O G I C A L S


!            I N T E G E R S

INTEGER naph, inode, iaph, numnp, istart, iend

!            R E A L S
REAL(Kind=WP_Kind) f2

!            D I M E N S I O N S

DIMENSION f2(numnp, 100)

!            C O M M O N S

!            I N T R I N S I C    F U N C T I O N S
!            E X T E R N A L S


iend = min(naph+1-istart,100)
READ(FEM_SF_Bin_Unit)((f2(inode,iaph),inode=1,numnp),iaph=1,iend)
RETURN
END
