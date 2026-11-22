SUBROUTINE CalcPot (nang, rho, ThetAPH, chiaph, v)

!            P U R P O S E   O F    S U B R O U T I N E
! This routine calculates the interaction potential at all of the
! quadrature points.

!            I N P U T    A R G U M E N T S
!  nang      number of APH theta and chi angles.
!  rho       hyperradius
!  ThetAPH   an array of APH theta angles.
!  chiaph    an array of APH chi angles.

!            O U T P U T    A R G U M E N T S
!  v         an array of interaction potential values.

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

USE Numeric_Kinds_Module
USE FileUnits_Module
USE Masses_Module
IMPLICIT NONE


!            I N T E G E R S
INTEGER  iang, nang

!            R E A L S
REAL(Kind=WP_Kind) rho, tmurs, potef

!            D I M E N S I O N S
REAL(Kind=WP_Kind) ThetAPH(nang), chiaph(nang), v(nang)

!            E X T E R N A L S
EXTERNAL popt, poteff

!-----------------------------------------------------------------------
! Determine printing options.
!-----------------------------------------------------------------------
LOGICAL little, medium, full
INTEGER ithcall, ithsub
DATA ithcall /0/, ithsub /0/
DATA little /.false./, medium /.false./, full /.false./
CALL popt ('calcpot ', little, medium, full, ithcall, ithsub)

!-----------------------------------------------------------------------
! calculate interaction potential.
!-----------------------------------------------------------------------

tmurs = usys2*rho*rho

! Pre-warm in serial: the potential chain performs one-time lazy
! initialization (file i/o, RKHS kernel loading) on its first call, which is
! not thread-safe.  Evaluating the first quadrature point here - its value is
! recomputed in the loop below - guarantees that the initialization has
! completed before the threaded loop starts.  All remaining per-call state in
! the chain (PotAPH distance scratch, RKHS gamma buffers) is local, so the
! loop itself is race-free.
CALL poteff(rho, ThetAPH(1), chiaph(1), potef )

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(iang, potef)
DO iang = 1, nang
   CALL poteff(rho, ThetAPH(iang), chiaph(iang), potef )
   v(iang) = potef/tmurs
ENDDO
!$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!  Print potential IF desired.
!-----------------------------------------------------------------------
IF(little)THEN
   WRITE(Out_Unit,*)'routine calcpot'
   WRITE(Out_Unit,*)'hyperradius = ',rho
   WRITE(Out_Unit,*)'potential calulated at ',nang,' angles.'
ENDIF
IF(medium)THEN
   WRITE(Out_Unit,*)' iang, ThetAPH, chiaph,    v'
   DO iang = 1, nang
      WRITE(Out_Unit,*)iang,ThetAPH(iang),chiaph(iang),v(iang)
   ENDDO
ENDIF

RETURN
ENDSUBROUTINE CalcPot
