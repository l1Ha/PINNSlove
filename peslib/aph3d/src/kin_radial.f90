SUBROUTINE Kin_Radial(n, toe, kinetic, kmax)
USE Numeric_Kinds_Module
USE FileUnits_Common_Module
USE FileUnits_OneDim_Module
IMPLICIT NONE
LOGICAL, PARAMETER:: idbug=.false.
INTEGER i, j, n, index, kmax
REAL(KIND=WP_Kind) toe(0:kmax), kinetic(n,n)

WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Begin(Kin_Radial)'

kinetic=0.D0
DO i=1,n
   DO j=1,n
      index=iabs(i-j)
      IF(index.lt.kmax)THEN
         kinetic(i,j)=kinetic(i,j)+toe(index)
      ENDIF
   ENDDO
ENDDO

IF(idbug)THEN
   WRITE(Dbug_Unit,*)'Radial-Kinetic Energy Matrix'
   CALL MxOut_Unit(kinetic, n, n, Dbug_Unit)
ENDIF
WRITE(Out_Unit,*)'End(Kin_Radial)'
WRITE(Out_Unit,*)

RETURN
END SUBROUTINE Kin_Radial
