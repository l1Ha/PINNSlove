MODULE CriticalDist_Module
   USE Numeric_Kinds_Module
   IMPLICIT NONE
   SAVE
   
   REAL(KIND=WP_Kind) re         ! Location of Potential Minimum
   REAL(KIND=WP_Kind) v0re       ! V(re)
   REAL(KIND=WP_Kind) v1re       ! dV(r)/dr evaluated at r=re
   REAL(KIND=WP_Kind) v2re       ! d2V(r)/dr2 evaluated at r=re
   REAL(KIND=WP_KIND) rmin(0:1000)! re for each electronic state
   REAL(KIND=WP_KIND) vmin(0:1000)! V(re) for each electronic state
   
   REAL(KIND=WP_Kind) rc         ! Location of Barrier Maximum
   REAL(KIND=WP_Kind) v0rc       ! V(rc)
   REAL(KIND=WP_Kind) v1rc       ! dV(r)/dr evaluated at r=rc
   REAL(KIND=WP_Kind) v2rc       ! d2V(r)/dr2 evaluated at r=rc
   REAL(KIND=WP_KIND) vmax(0:1000)! V(rc) for each electronic state
   
   REAL(KIND=WP_Kind) rb         ! Inner inflection point
   REAL(KIND=WP_Kind) v0rb       ! V(rb)
   REAL(KIND=WP_Kind) v1rb       ! dV(r)/dr evaluated at r=rb
   REAL(KIND=WP_Kind) v2rb       ! d2V(r)/dr2 evaluated at r=rb
   
   REAL(KIND=WP_Kind) rd         ! Outer inflection point
   REAL(KIND=WP_Kind) v0rd       ! V(rd)
   REAL(KIND=WP_Kind) v1rd       ! dV(r)/dr evaluated at r=rd
   REAL(KIND=WP_Kind) v2rd       ! d2V(r)/dr2 evaluated at r=rd
   
   REAL(KIND=WP_Kind) sigma0     ! Inner point where V(r) crosses it asymptotic value
   REAL(KIND=WP_Kind) v0sigma0   ! V(sigma0)
   REAL(KIND=WP_Kind) v1sigma0   ! dV(r)/dr evaluated at r=sigma0
   REAL(KIND=WP_Kind) v2sigma0   ! d2V(r)/dr2 evaluated at sigma0
   
   REAL(KIND=WP_Kind) sigma1     ! Outer point where V(r) crosses it asymptotic value
   REAL(KIND=WP_Kind) v0sigma1   ! V(sigma1)
   REAL(KIND=WP_Kind) v1sigma1   ! dV(r)/dr evaluated at r=sigma1
   REAL(KIND=WP_Kind) v2sigma1   ! d2V(r)/dr2 evaluated at r=sigma1
   
   REAL(KIND=WP_Kind) ri         ! Smallest for which the potential is indistinguishable from its asympotic value
   REAL(KIND=WP_Kind) vasy       ! Asymptotic value of the potential
   REAL(KIND=WP_KIND) vinf(0:1000)! vasy for each electronic state
   
   
   REAL(KIND=WP_Kind) k, omega, alpha, de, be, beta
   
   
END MODULE CriticalDist_Module