MODULE System_Module
   USE Numeric_Kinds_Module
   USE DVRMasses_Module
   USE Boundary_Module
   IMPLICIT NONE
   SAVE
   CHARACTER(LEN=3) :: pointg = 'c2v'       ! point group of PES
   LOGICAL :: debug   = .False.                ! LOGICAL variable that print out lots of stuff IF you set them to true.
   LOGICAL :: medium  = .False.                ! LOGICAL variable. .true. prints a few matrix elements.
   LOGICAL :: matelem =.True.               ! IF true, program calcs matrix elements.
   LOGICAL :: h3sys = .False.                 ! h3sys=.false. except for the system H3.
   LOGICAL :: logspace = .True.              ! IF true, rho spacing is a constant fraction of rho.
   LOGICAL :: peigv = .False.                 ! LOGICAL variable that print out lots of stuff IF you set them to true.
   LOGICAL :: fbrtr = .False.                 ! IF true, dvr to fbr transformation is performed.
   INTEGER :: nthesec(10) = [20,   50 ,   60,   0, 0, 0, 0, 0, 0, 0]          ! order of legendre polynomial basis in theta
   INTEGER :: nchisec(10) = [40,   100,  120,   0, 0, 0, 0, 0, 0, 0]           ! number of dvr chi points for this representation
   INTEGER :: ltheta = 0                ! degree of legendre polynomial basis (????)
   INTEGER :: nsfunc = 136                ! number of surface functions saved.
   INTEGER :: nrow = 5                  ! Not used
   INTEGER :: ncol = 5                  ! Not used
   INTEGER :: ngood1 = 25                ! ngood1<ngood2 are the numbers of surface function energies included in energy averages used i
   INTEGER :: ngood2 = 25                ! ngood1<ngood2 are the numbers of surface function energies included in energy averages used i
   INTEGER :: Num_DVR_Sectors = 3       ! number of sectors with different ntheta and nchi.
   INTEGER :: initsec = 1               ! ntheta, nchi sector of initial rho.
   INTEGER :: nrhomx  = 5               ! number of rho values wanted.
   INTEGER :: nfbr = 200                  ! the number of surface functions written to one file.
   REAL(Kind=WP_Kind) :: rendsec(10) = [6.0d0,8.0d0,12.0d0,0.0d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0] ! ends of sectors where ntheta and nchi change.
   REAL(Kind=WP_Kind) :: potmsec(10) = [3.9d0,3.9d0, 3.9d0,0.0d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0]! maximum potential energy for which points are kept
   REAL(Kind=WP_Kind) :: pevmsec(10) = [3.5d0,2.6d0, 2.6d0,0.0d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0]! max energy for which 1-d eigenvectors are kept
   REAL(Kind=WP_Kind) :: rhofbr1 = 11.80d0    ! beginning of rho range where transf. to fbr is done.
   REAL(Kind=WP_Kind) :: rhofbr2 = 11.99d0    ! END of rho range where transformation to fbr is done.
   REAL(Kind=WP_Kind) :: fmpotma = 1.500d0     ! ratio of potmax to pev2dmax IF nrhos1>=2
   REAL(Kind=WP_Kind) :: deltar = 0.01d0     ! rho spacing.  Absolute or relative depending on logspace.
   REAL(Kind=WP_Kind) :: Rho1 = 5.d0       ! first rho
   NAMELIST/system/nthesec,nchisec,rendsec,peigv,debug, medium,potmsec,pevmsec,ltheta,             &
                   rho1,matelem,h3sys,fbrtr,rhofbr1,rhofbr2, nsfunc,nrow,ncol,ngood1,ngood2,pointg,     &
                   Num_DVR_Sectors,initsec,nrhomx,fmpotma,deltar,logspace, nfbr
   END MODULE System_Module