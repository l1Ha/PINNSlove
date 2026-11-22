SUBROUTINE COMPL(BETA, BETAPI, NU, LAMBDA, S, MAXNU, MINNU, PMAX, LGEB)
USE Numeric_Kinds_Module
SAVE
!
!   $RCSfile: compl.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:05 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE -
!
!         COMPUTES MAX AND MIN(NU(I)), I=1, ..., P, LAMBDA, AND BETAPI
!         FOR ALL NOT TOO SMALL, NU-VALUES.
!
!
!     OUTPUT PARAMETERS -
!
!         LGEB  = TRUE IF A CONVERGED EIGENVALUE IS >= B (SUCH
!                 EIGENVALUES ARE THROWN AWAY), AND FALSE OTHERWISE.
!                 IF LGEB = T, IT WILL CAUSE THE NEXT SHIFT TO BECOME B.
!         MAXNU = MAX(NU(1), ..., NU(P))
!         MINNU = MIN(NU(1), ..., NU(P))
!
!         IF ABS(NU(I)) <= TOLZNU THE VALUE IS NOT USED.
!
!     PLEASE SEE THE PROGRAMMERS GUIDE FOR INFORMATION ABOUT
!     PARAMETERS NOT EXPLAINED ABOVE, AND FOR MORE DETAILS ABOUT
!     THE FUNCTION OF THE ROUTINE.
!
!
! **********************************************************************
!
INTEGER   PMAX
REAL(Kind=WP_Kind)  ALTMU, B, BETA(1), BETAPI(1), FACTOR, LAMBDA(1), SRELPR, MAXNU, MINNU, MU, NEXTMU, NU(1), NUI, OLDMU, S(PMAX, 1)
REAL(Kind=WP_Kind) TOLBPI, TOLLDL, TOLPDM, TOLS1I, TOLZBT, TOLZNU, A, BB, CUT
INTEGER   CNEG, CPOS, I, ITERNO, N, OLCPOS, P, REST, RNEW, ROLD,  TCONV, NUMEIG, MAXW, MAXIW
LOGICAL   USEMX, ZERBET, REACHB, USSMXR, USEDB, T, F, LGEB
COMMON   /STLMCT/ N, ITERNO, TCONV, CNEG, CPOS, OLCPOS, RNEW, ROLD, REST, P, USEMX, ZERBET
COMMON   /STLMIN/ A, B, NUMEIG, MAXW, MAXIW
COMMON   /STLMMU/ MU, OLDMU, NEXTMU, ALTMU
COMMON   /STLMTL/ SRELPR, TOLBPI, TOLS1I, FACTOR, TOLPDM, TOLZBT, TOLZNU, TOLLDL(3)
COMMON   /STLMTF/ T, F
COMMON   /STLMTS/ REACHB, USSMXR, USEDB
!
!     ********************************************************
!     LOOK FOR A NU-VALUE WHERE 1/NU WILL BE SAFE TO COMPUTE.
!     INITIALIZE MINNU AND MAXNU.
!     ********************************************************
DO 10 I = 1, P
  IF(.NOT. ABS(NU(I)) > TOLZNU) GOTO 10
    MINNU = NU(I)
    MAXNU = MINNU
    GOTO 20
10      CONTINUE
!
MINNU = -1.0D0
MAXNU = -1.0D0
!
20    BB = BETA(P)
!     ****************************************************************
!     THE ONLY CASE WHEN CUT /= B, IS WHEN THE LDL(T) FACTORIZATION
!     FAILED FOR MU = B.
!     ****************************************************************
CUT = B
LGEB = F
IF(REACHB) CUT = MU
DO 40 I = 1, P
  NUI = NU(I)
  IF(.NOT. ABS(NUI) > TOLZNU) GOTO 30
!         **********************************************************
!         FOR ALL NU-VALUES NOT NEAR ZERO, COMPUTE THE CORRESPONDING
!         LAMBDA AND CONVERGENCE VALUES.
!         UPDATE MINNU AND MAXNU.
!         **********************************************************
    LAMBDA(I) = MU + 1.0D0 / NUI
    BETAPI(I) = BB * S(P, I) / NUI
    IF(.NOT. (LAMBDA(I) > CUT  .AND.  ABS(BETAPI(I)) <= TOLBPI)  ) GOTO 25
!           **********************************************************
!           THROW AWAY CONVERGED EIGENVALUES GREATER THAN CUT, AND SET
!           LGEB.
!           **********************************************************
      BETAPI(I) = -1.0D20
      LGEB      = T
!
25        MINNU = MIN(MINNU, NUI)
    MAXNU = MAX(MAXNU, NUI)
    GOTO 40
!
!       **************************************************
!       PUT IN VALUES SO THAT WRINFO WILL WORK.
!       BETAPI=-1E20 MEANS THAT IT WILL NEVER BE ACCEPTED.
!       **************************************************
30      LAMBDA(I) = -1.0D20
  BETAPI(I) = -1.0D20
!
40      CONTINUE
!
RETURN
END
