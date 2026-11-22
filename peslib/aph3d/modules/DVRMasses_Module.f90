MODULE DVRMasses_Module
   USE Numeric_Kinds_Module
   USE Masses_Module
   IMPLICIT NONE
   SAVE
   REAL(Kind=WP_Kind) scale(3)   !
   REAL(Kind=WP_Kind) sine(3)    !
   REAL(Kind=WP_Kind) cosn(3)    !
   REAL(Kind=WP_Kind) tang(3)    !
   REAL(Kind=WP_Kind) rot(3)     !
   REAL(Kind=WP_Kind) eta(3)     !
   REAL(Kind=WP_Kind) zeta(3)    !
   REAL(Kind=WP_Kind) rmatch(3)  !
   REAL(Kind=WP_Kind) rstar(3)   !
   REAL(Kind=WP_Kind) dstar(3)   !
   REAL(Kind=WP_Kind) req(3)     !
   REAL(Kind=WP_Kind) redmas     !
   REAL(Kind=WP_Kind) twomu      !
   REAL(Kind=WP_Kind) qchan(3)   !
   REAL(Kind=WP_Kind) dmatch(3)  !
   REAL(Kind=WP_Kind) rmass(3)   !
   REAL(Kind=WP_Kind) xstar(3)   !
END MODULE DVRMasses_Module
