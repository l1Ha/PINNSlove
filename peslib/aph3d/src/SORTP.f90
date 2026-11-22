SUBROUTINE SORTP(POINTR, LAMBDA, LAST)
USE Numeric_Kinds_Module
!
!   $RCSfile: sortp.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:43 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE -
!
!         THIS IS A USUAL SELECTION SORT WITH ONE DIFFERENCE. IT
!         SORTS THE POINTR VECTOR WITH RESPECT TO LAMBDA (I.E. WE
!         CHANGE POINTR AND NOT LAMBDA).
!
!
!     INPUT PARAMETERS -
!
!         LAMBDA       = LAMBDA EIGENVALUES.
!         POINTR(LAST) = NUMBER OF ELEMENTS TO BE SORTED.
!
!
!     OUTPUT PARAMETERS -
!
!         POINTR IS SORTED SO THAT  (IF P =POINTR(LAST))
!         LAMBDA(IABS(POINTR(1))), ..., LAMBDA(IABS(POINTR(P))) ARE IN
!         ASCENDING ORDER.
!
! **********************************************************************
!
REAL(Kind=WP_Kind)     LAMBDA(1), TLAM
INTEGER   I, IP1, J, K, LAST, P, PM1, POINT, POINTR(1)
!
P = POINTR(LAST)
PM1 = P - 1
IF(PM1 <= 0) GOTO 9999
!
DO 20 I = 1, PM1
  K = I
  POINT = IABS(POINTR(K))
  TLAM = LAMBDA(POINT)
  IP1 = I + 1
!
  DO 10 J = IP1, P
    POINT = IABS(POINTR(J))
!         ********
!         COMPARE.
!         ********
    IF(.NOT. LAMBDA(POINT) < TLAM) GOTO 10
     K = J
      POINT = IABS(POINTR(K))
      TLAM = LAMBDA(POINT)
10        CONTINUE
!
!       ***************************
!       INTERCHANGE POINTER VALUES.
!       ***************************
  IF(.NOT. I /= K) GOTO 20
    POINT = POINTR(I)
    POINTR(I) = POINTR(K)
    POINTR(K) = POINT
!
20      CONTINUE
!
9999  RETURN
END
