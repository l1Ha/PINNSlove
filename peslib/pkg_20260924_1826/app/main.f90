!     File created at Fri Jun  5 21:58:55 PDT 2020
!     Original source code: abc.f

      program abc
         use mod_driver,only:driver 
         use mod_ParSet,only:parset 
         use mod_potenl,only:writePES 
         use mod_pes_multi,only:init_pes
         use mod_common
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     CCP6 hyperspherical coordinate reactive scattering program ABC.
!     This version dated 31 March 2000.
!     -----------------------------------------------------------------
!
!     common blocks
!
!
!     input parameters
!
   integer :: NPES = 1
   character(len=20) :: PESMODE = 'INTERNAL_1STATE'
   logical :: write_pes = .false.
    namelist /input/ mass,jtot,ipar,jpar,jmax,kmax,rmax,mtr,emax,enrg,dnrg,nnrg, &
                     nout,jout,nmax,NPES,PESMODE,pes_threshold_factor,         &
                     pes_threshold_w,pes_C6,write_pes,mro_override,            &
                     mvi_override,nvi_override,rmin_override,smax_override,    &
                     ireact,rdat,smin_asym_override,iscale,ecol_min,ecol_max
    read(5,nml=input)
    call init_pes(NPES,PESMODE)
! 
!     universal constants
!     (energies in eV, lengths in bohr, masses in amu)
!     
      
      do ia = 1,3
         imass = int(mass(ia))
         if (imass .eq.  1) then
            mass(ia) = 1.007825d0
         else if (imass .eq.  2) then
            mass(ia) = 2.014000d0
         else if (imass .eq. 4) then 
            mass(ia) = 4.002600d0 
          else if (imass .eq. 19) then
             mass(ia) = 18.99840d0
          else if (imass .eq. 39) then
             mass(ia) = 38.9637064864d0
          else if (imass .eq. 40) then
             mass(ia) = 39.963998166d0
          else if (imass .eq. 41) then
             mass(ia) = 40.9618252579d0
          else if (imass .eq. 85) then
             mass(ia) = 84.9117897379d0
          else if (imass .eq. 87) then
             mass(ia) = 86.9091805310d0
          else if (imass .eq. 35) then
             mass(ia) = 34.96885d0
          else if (imass .eq. 37) then
            mass(ia) = 36.96590d0
         else
            stop 'abc 1'
         endif
      enddo
      hbarsq = 0.014927625d0
!
!     scale factors for MSJ coordinates
!
      mtot = mass(1)+mass(2)+mass(3)
      mred = sqrt(mass(1)*mass(2)*mass(3)/mtot)
      scales(1) = sqrt((mass(1)/mred)*(1.d0-mass(1)/mtot))
      scales(2) = sqrt((mass(2)/mred)*(1.d0-mass(2)/mtot))
      scales(3) = sqrt((mass(3)/mred)*(1.d0-mass(3)/mtot))
      rmlmda = 2.0d0*mred/hbarsq
!
!     angular momentum projections
!
      if ((-1)**jtot .eq. ipar) then
         kmin = 0
      else
         kmin = 1
      endif
      kmax = min(jtot,jmax,kmax)
      if (mass(2) .ne. mass(3)) jpar = 0
      nout = max(nout,0)
      if (jpar .eq. -1) then
         jout = max(jout,1)
      else
         jout = max(jout,0)
      endif
!
!     input summary
!
      write (6,61) mass,jtot,ipar,jpar,ireact,jmax,kmin,kmax,rmax,mtr,emax, &
                   enrg,dnrg,nnrg,nout,jout,nmax

  61  format(/1x,'INPUT:'/1x,70('-')/1x, &
      '  mass = ',3f16.3/1x, &
      '  jtot = ',i32/1x, &
      '  ipar = ',i32/1x, &
      '  jpar = ',i32/1x, &
      ' ireact = ',i32/1x, &
      '  jmax = ',i32/1x, &
      '  kmin = ',i32/1x, &
      '  kmax = ',i32/1x, &
      '  rmax = ',f32.3/1x, &
      '   mtr = ',i32/1x, &
      '  emax = ',f32.3/1x, &
      '  enrg = ',f32.3/1x, &
      '  dnrg = ',f32.3/1x, &
      '  nnrg = ',i32/1x, &
      '  nout = ',i32/1x, &
      '  jout = ',i32/1x, &
      '  nmax = ',i32/1x,70('-'))
!
   PRINT*,'calculation'
!     
      if (write_pes) call writePES()
      if (kmin .gt. kmax) stop 'abc 2'
      call parset
      call driver (nmax)
!     stop
      end
