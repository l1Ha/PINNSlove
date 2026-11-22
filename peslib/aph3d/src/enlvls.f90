SUBROUTINE enlvls(r, n, w)

!         P U R P O S E   O F    S U B R O U T I N E
! Calculates the adiabatic energy levels.

!         I N P U T    A R G U M E N T S
!  r    Current propagation distance.
!  w    Interaction matrix to be diagonalized. The eigenvalues are
!    the adabatic energy levels.
!  n    Dimension of the matrix w.

!         O U T P U T    A R G U M E N T S
! NONE

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

USE FileUnits_Module
USE nstate_Module
USE doenlvls_Module
USE TotalEng_Module
USE Masses_Module
USE convrsns_Module
USE engpro_Module
IMPLICIT NONE
INTEGER n, ithcall, ithsub, nn2, i
LOGICAL little, medium, full
REAL(Kind=WP_Kind) w(n*n), temp(nstate*(nstate+1)/2), elevls(nstate), scr(nstate), vec(nstate*nstate)
REAL(Kind=WP_Kind) r, xkin
INTEGER, PARAMETER :: ione=1
!------------------------------------------------------------------
! Obtain printing options.
!------------------------------------------------------------------
DATA ithcall/0/, ithsub/0/
DATA little/.false./, medium/.false./, full/.false./
CALL popt('enlvls  ', little, medium, full, ithcall, ithsub)
!------------------------------------------------------------------
!      WRITE(Msg_Unit,*)'here*********************',lenlvls,kenergy
IF(.NOT.lenlvls.or.kenergy>1)RETURN

!trp>> Following line commented out 9/12/97 to conform to Sefano's code. c<<trp
!      xkin=usys2*Etot
IF(little)WRITE(Out_Unit,*)'xkin=',xkin

!------------------------------------------------------------------
! Store the interaction matrix in temp so that w will not be
! overwritten.
!------------------------------------------------------------------

! revised by X.L.  it looks wrong
nn2=n*(n+1)/2
DO 1 i=1,nn2
    temp(i)=w(i)
1 CONTINUE

!------------------------------------------------------------------
! Diagonalize the interaction matrix to obtain the adiabatic
! energy levels.
!------------------------------------------------------------------
CALL tdiag(temp,elevls,vec,scr,n)
IF(full)WRITE(Out_Unit,*)'temp before=',temp
!------------------------------------------------------------------
! Convert to electron volts (eV).
!trp>> modified 9/12/97 to conform to Stefano's version <<trp
!------------------------------------------------------------------
DO 2 i=1,n
!         elevls(i)=(elevls(i)+xkin)*autoev/usys2
    elevls(i)=elevls(i)*autoev/usys2
2 CONTINUE
IF(full)WRITE(Out_Unit,*)'elevls after=',elevls

!------------------------------------------------------------------
! Write out adabatic energy levels to a file for later plotting.
!------------------------------------------------------------------
WRITE(elev_unit,5)
WRITE(elev_unit,6)ione
WRITE(elev_unit,3)r
WRITE(elev_unit,7)
WRITE(elev_unit,6)n
WRITE(elev_unit,4)(elevls(i),i=1,n)
WRITE(elev_unit,8)

3 FORMAT(/,f10.5)
4 FORMAT(/,e14.7)
5 FORMAT(//,' $readx',5x,'distance in bohr')
6 FORMAT(/,i5)
7 FORMAT(' $ready',5x,'energy in ev.')
8 FORMAT(//,' $CONTINUE')

RETURN
END
