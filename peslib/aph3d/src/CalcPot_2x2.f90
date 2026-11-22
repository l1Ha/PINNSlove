SUBROUTINE CalcPot_2x2 (nang, rho, ThetAPH, chiaph, v1, v2, v12)

!            P U R P O S E   O F    S U B R O U T I N E
! Two-electronic-state (2x2) analogue of CalcPot.  Calculates the three
! independent elements of the 2x2 diabatic potential at all of the
! quadrature points, normalized in the APH convention V/(usys2*rho**2).

!            I N P U T    A R G U M E N T S
!  nang      number of APH theta and chi angles.
!  rho       hyperradius
!  ThetAPH   an array of APH theta angles.
!  chiaph    an array of APH chi angles.

!            O U T P U T    A R G U M E N T S
!  v1, v2, v12   arrays of the 2x2 potential matrix elements.

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

USE Numeric_Kinds_Module
USE FileUnits_Module
USE Masses_Module
IMPLICIT NONE


!            I N T E G E R S
INTEGER  iang, nang

!            R E A L S
REAL(Kind=WP_Kind) rho, tmurs, p11, p22, p12

!            D I M E N S I O N S
REAL(Kind=WP_Kind) ThetAPH(nang), chiaph(nang), v1(nang), v2(nang), v12(nang)

!            E X T E R N A L S
EXTERNAL popt, PotAPH_2x2

!-----------------------------------------------------------------------
! Determine printing options.
!-----------------------------------------------------------------------
LOGICAL little, medium, full
INTEGER ithcall, ithsub
DATA ithcall /0/, ithsub /0/
DATA little /.false./, medium /.false./, full /.false./
CALL popt ('calcpot2', little, medium, full, ithcall, ithsub)

!-----------------------------------------------------------------------
! calculate the 2x2 potential matrix on the angular grid.
!-----------------------------------------------------------------------

tmurs = usys2*rho*rho

! Pre-warm in serial: the potential chain performs one-time lazy
! initialization (file i/o, RKHS kernel loading) on its first call, which is
! not thread-safe.  Evaluating the first quadrature point here - its value is
! recomputed in the loop below - guarantees that the initialization has
! completed before the threaded loop starts.  All remaining per-call state in
! the chain (PotAPH_2x2 distance scratch, RKHS gamma buffers) is local, so the
! loop itself is race-free.
CALL PotAPH_2x2(rho, ThetAPH(1), chiaph(1), p11, p22, p12)

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(iang, p11, p22, p12)
DO iang = 1, nang
   CALL PotAPH_2x2(rho, ThetAPH(iang), chiaph(iang), p11, p22, p12)
   v1(iang)  = p11/tmurs
   v2(iang)  = p22/tmurs
   v12(iang) = p12/tmurs
ENDDO
!$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!  Print potential IF desired.
!-----------------------------------------------------------------------
IF(little)THEN
   WRITE(Out_Unit,*)'routine calcpot_2x2'
   WRITE(Out_Unit,*)'hyperradius = ',rho
   WRITE(Out_Unit,*)'2x2 potential calculated at ',nang,' angles.'
ENDIF
IF(medium)THEN
   WRITE(Out_Unit,*)' iang, ThetAPH, chiaph,    v1,    v2,    v12'
   DO iang = 1, nang
      WRITE(Out_Unit,*)iang,ThetAPH(iang),chiaph(iang),v1(iang),v2(iang),v12(iang)
   ENDDO
ENDIF

RETURN
ENDSUBROUTINE CalcPot_2x2
