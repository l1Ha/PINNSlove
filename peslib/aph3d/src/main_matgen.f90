SUBROUTINE main_matgen
!
!

!            P U R P O S E   O F    S U B R O U T I N E
! This is the main program for the APH 3-Dimensional Reactive
! scattering propagator code.


!            I N P U T    A R G U M E N T S
!  NONE

!            O U T P U T    A R G U M E N T S
!  NONE

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
USE FileUnits_Module
USE nstate_Module
IMPLICIT NONE
CHARACTER(LEN=10) hour
CHARACTER(LEN=8)  today
CHARACTER(LEN=5) curzone
CHARACTER(LEN=200) pmatfile, ovrfile
INTEGER dtvalues(8)
INTEGER nidm
REAL(Kind=WP_Kind) rmat(nstate,nstate),  eye11(nstate,nstate), w(nstate,nstate)
REAL(Kind=WP_Kind) vecnew(nstate,nstate), tstore(nstate,nstate)
REAL(Kind=WP_Kind) vecold(nstate,nstate),    xsq(nstate), xk(nstate)
EXTERNAL dateh, timeh, title, constant, readaph_ph1, aph3d1, open_matph

WRITE(Out_Unit,*)'Called    Main_Matgen'
CLOSE(Out_Unit)
OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/Main_Matgen.txt',Form='FORMATTED')

!-----------------------------------------------------------------------
! Determine the time and todays date.
!-----------------------------------------------------------------------
CALL dateh(today)
CALL timeh(hour)

!-----------------------------------------------------------------------
!   Write the time and date to unit Out_Unit.
!-----------------------------------------------------------------------
WRITE(Out_Unit,1111)today
WRITE(Out_Unit,2222)hour
CALL Date_And_Time(today, hour, curzone, dtvalues)
CALL Print_Date_And_Time(Today, Hour, CurZone, DtValues, Out_Unit)

!-----------------------------------------------------------------------
!   Write out the title page.
!-----------------------------------------------------------------------
CALL title

!-----------------------------------------------------------------------
! Establish all conversion factors and physical constants
!-----------------------------------------------------------------------

CALL constant

!-----------------------------------------------------------------------
! Read in all input data
!-----------------------------------------------------------------------
CALL readaph_ph1 (pmatfile, ovrfile)
WRITE(Out_Unit,*)'main: pmatfile=',pmatfile
WRITE(Out_Unit,*)'main: ovrfile =',ovrfile
WRITE(Msg_Unit,*)'main: pmatfile=',pmatfile,'main: ovrfile =',ovrfile

!-----------------------------------------------------------------------
! Open files for input and output of data.
!-----------------------------------------------------------------------

CALL open_matph (pmatfile, ovrfile)

WRITE(Out_Unit,*)'main after open_matph: pmatfile=',pmatfile
WRITE(Out_Unit,*)'main after open_matph: ovrfile =',ovrfile

!-----------------------------------------------------------------------
! Propagate the set of coupled radial equations to determine transition
! probabilities and S-matrices.
!-----------------------------------------------------------------------
REWIND pmat_unit
REWIND Ovr_Unit
REWIND 22
REWIND 23
!
REWIND fph_unit
!
nidm=5

CALL aph3d1 (rmat, eye11, w, vecnew, tstore, vecold, xsq, xk, nidm)
! PORT (2026-10-02): the classic propagation consumes the per-sector records
! (Ovr_Unit/pmat_unit, including the Delves-region loop); rewind them so the
! following per-energy propagation (logprop/propph1) reads from the start.
REWIND Ovr_Unit
REWIND pmat_unit
REWIND frst_unit

1111  FORMAT(40x,'date: ',a12)
2222  FORMAT(40x,'time: ',a10)

WRITE(Out_Unit,*)'Completed Main_Matgen'
WRITE(Out_Unit,*)
CLOSE(Out_Unit)

OPEN(Unit=Out_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Output/BasisForScatt.txt',Form='FORMATTED',ACCESS='APPEND')
WRITE(Out_Unit,*)'Completed Main_Matgen'
WRITE(Out_Unit,*)

RETURN
END
