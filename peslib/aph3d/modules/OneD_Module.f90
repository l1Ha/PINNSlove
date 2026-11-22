MODULE OneD_Module
   USE Numeric_Kinds_Module
   USE DVR_Module
   IMPLICIT NONE
   SAVE
   INTEGER(KIND=IW_Kind) nev(nthetmax)
   INTEGER(KIND=IW_Kind) nptchi(nmax,nthetmax)         ! keeps track tells us that basis point nctem corresponds to dvr point ichi at t
   INTEGER(KIND=IW_Kind) nc2tot(nthetmax)              ! number of dvr points kept at each theta
   REAL(Kind=WP_Kind) cc1dtot(nmax,nmax,nthetmax)   ! contains 1-d eigenvectors
   REAL(Kind=WP_Kind) s1dtot(nmax,nmax)             ! not used
   REAL(Kind=WP_Kind) h1d(nmax,nmax,nthetmax)       ! 1-d hamiltonian matrices
   REAL(Kind=WP_Kind) e1d(nmax,nthetmax)            ! 1-d eigenvalues
   REAL(Kind=WP_Kind) e2d(nt2max)                   ! 2-d eigenvalues
   REAL(Kind=WP_Kind) cc1d(nmax,nmax,nthetmax)      ! cc1d--explained in chctrans
END MODULE OneD_Module
