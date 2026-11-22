SUBROUTINE SETUP2 ( PRINT, HBARSQ )
USE DVRMasses_Module
USE Numeric_Kinds_Module
USE Numbers_Module
INTEGER I, J, K
LOGICAL       PRINT
REAL(Kind=WP_Kind) Hbarsq, d, r0, a, b, c, e, q, x, de, beta, delta, vbase

COMMON / LEPS /   DE(3), BETA(3), R0(3), DELTA(3), VBASE
!
!       VBASE had been used only in defining QCHAN in RXN3D.
!       It has been necessary to allow VBASE to change the zero
!       of energy in SLICE functions.
!
!     VBASE = -MAX ( DE(1), DE(2), DE(3) )
!

MTOT = MASS(1) + MASS(2) + MASS(3)
REDMAS = SQRT ( MASS(1)*MASS(2)*MASS(3)/MTOT )
TWOMU = 2D0 * REDMAS * HBARSQ
DO 1 I = 1, 3
J=MOD(I,3)+1
K=MOD(J,3)+1
SCALE(I) = ( MTOT*MASS(J)*MASS(K)/(MASS(I)*(MASS(J)+MASS(K))**2))** .25D0
D = (MASS(I) + MASS(J))*(MASS(I) + MASS(K))
COSN(I) = -SQRT(MASS(J)*MASS(K)/D)
SINE(I) =  SQRT(MTOT*MASS(I)/D)
TANG(I) = SINE(I)/COSN(I)
ROT(I) = ATAN2 ( SINE(I), COSN(I) )
IF( ROT(I) < 0D0 ) ROT(I) = ROT(I) + PI + PI
!
! 5/85 MODIFICATION to allow ZETA to be input to specify polar region.
!
IF(ZETA(I)<=0.0D0)THEN
ETA(I) = PI - ROT(I)
ZETA(I) = .5D0*PI - ETA(I)
RSTAR(I) = RMATCH(I)*SINE(I)
DSTAR(I) = -RMATCH(I)*COSN(I)
ELSE
ETA(I) = .5D0*PI - ZETA(I)
RSTAR(I) = RMATCH(I)*COS(ZETA(I))
DSTAR(I) = RMATCH(I)*SIN(ZETA(I))
ENDIF
! 5/85 END MODIFICATION
REQ(I) = R0(I)*SCALE(I)
A=MASS(I)*MASS(I)
B=MASS(J)*MASS(J)
C=MASS(K)*MASS(K)
D=MASS(J)-MASS(K)
E=MASS(J)+MASS(K)
RMASS(I)=MASS(J)/E
DMATCH(I)=SCALE(I)*E
Q=A+D*D+2D0*MASS(I)*(B+C)/E
Q=SQRT(Q)
X=MTOT*D/(Q*E)
XSTAR(I)=ACOS(X)
    1 CONTINUE
IF( .NOT. PRINT ) RETURN
PRINT 1000, (I, MASS(I), SCALE(I),     SINE(I), COSN(I), ROT(I), ETA(I), ZETA(I), XSTAR(I), REQ(I), &
                 RMATCH(I), RSTAR(I), DSTAR(I),I=1,3)
PRINT 1001, (I, R0(I), DE(I), BETA(I), DELTA(I), I = 1, 3 )
 1000 FORMAT( 1X, I1, F5.1, F7.3, 2X, 2F8.4, 3F8.3, 5X, 5F9.4)
 1001 FORMAT( 1X, I1, 2X, 4F8.4)
RETURN
ENDSUBROUTINE Setup2
