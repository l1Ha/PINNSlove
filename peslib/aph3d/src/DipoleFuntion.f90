SUBROUTINE DipoleFunction(rvals, Dipole, nmax, ithelec, kthelec)
USE Numeric_Kinds_Module
USE FileUnits_Module
USE FileUnits_OneDim_Module
USE Numbers_Module, ONLY: Zero
IMPLICIT NONE
SAVE
LOGICAL IRead
INTEGER n, ithelec, kthelec, nmax, i
REAL(KIND=WP_Kind) r, rvals(0:nmax), Dipole(0:nmax), Rcalc(118), Dcalc(118,7)
REAL(KIND=WP_Kind), ALLOCATABLE :: Dinterp(:,:)
DATA IRead/.False./


IF(IRead)THEN
   ALLOCATE(Dinterp(0:nmax,7))
   OPEN(Unit=In_Unit,File=TRIM(InputDir)//"TransitionDipoleMoments.txt")
   IRead=.False.
   DO n=1,11
      READ(In_unit,*)
   ENDDO
   DO n=1,118
      READ(In_unit,*)Rcalc(n),(Dcalc(n,i),i=1,7)
      !WRITE(*,'(F6.3,7E14.6)')Rcalc(n),(Dcalc(n,i),i=1,7)
   ENDDO
   CLOSE(In_Unit)
   DO i=1,7
      CALL AKIMA(3,118,rcalc,dcalc(:,i),nmax+1,rvals, Dinterp(:,i))
   ENDDO
   OPEN(Unit=Dipole_Unit,File=TRIM(OutDIR)//"GraphicsOut/Dipole.csv")
   DO n=0,nmax
      WRITE(Dipole_Unit,'(8(e14.7,","))')rvals(n),(Dinterp(n,i),i=1,7)
   ENDDO
   CLOSE(Dipole_Unit)
ELSE
   ! PORT FIX (2026-10-01): Dinterp is ALLOCATABLE and unallocated on this
   ! early-return path; assigning to it caused a SIGSEGV.  Only Dipole needs
   ! clearing here.
   !Dinterp=Zero
   Dipole=Zero
   RETURN
ENDIF

!Dipole moment function for a specific electronic state
IF(ithelec.eq.kthelec)THEN
   DO n=0,nmax
      r=rvals(n)
      Dipole(n)=Dinterp(n,1)
   ENDDO
ELSE
!Transition dipole moment function for a transition
!from the jth electronic state to the kth electronic state
   DO n=0,nmax
      r=rvals(n)
      Dipole(n)=Dinterp(n,2)
   ENDDO
ENDIF

RETURN
END SUBROUTINE DipoleFunction