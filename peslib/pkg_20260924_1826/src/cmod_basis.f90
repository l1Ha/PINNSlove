!     File created at Fri Jun  5 21:58:55 PDT 2020
!     Original source code: basis.f
module mod_basis
      private
      public :: basis
      contains
      subroutine basis (cvi,cro,eint,ilev,jlev,klev,nlev,nmax)
         use mod_common
         use mod_pes_multi, only: PESMODE
         use mod_potenl, only: potenl
         use mod_symevp, only: symevp 
         use mod_basort ,only: basort 
         use mod_qvib ,only: qvib 
         use mod_bclean ,only: bclean
         use mod_scratchfile, only: getrec 
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine constructs an asymptotic multiple-arrangement
!     rovibrational basis set containing all channels with j.le.jmax
!     and eint.le.emax.
!     -----------------------------------------------------------------
!
!     common blocks (migrated to mod_common)
!
!     input arrays
!
      dimension cvi(nvi,nmax),cro(3,nmax),eint(nmax)
      dimension ilev(nmax),jlev(nmax),klev(nmax),nlev(nmax)
!
!     local arrays
!
!     dimension h(nvi,nvi),t(nvi),e(nvi)
      allocatable h(:,:),t(:),e(:)

!     dimension s(0:2*nvi),v(0:2*nvi)
      allocatable s(:),v(:)

!     dimension wvi(mvi),xvi(mvi)
      allocatable wvi(:),xvi(:)

      logical :: apply_floor
      real(8), parameter :: KRB2_SPURIOUS_MARGIN_EV = 0.1d0
      real(8) :: vmin_ref, e_floor, vraw, sa_ref, smin_eff
      integer :: nvphys

      ! runtime sanity checks to avoid absurd allocations from corrupted common
      if (nvi .lt. 1 .or. nvi .gt. 10000) then
         write(*,*) 'FATAL BASIS: invalid nvi =', nvi
         stop 'basis 2'
      endif
      if (mvi .lt. 1 .or. mvi .gt. 10000) then
         write(*,*) 'FATAL BASIS: invalid mvi =', mvi
         stop 'basis 3'
      endif
      allocate (h(nvi,nvi),t(nvi),e(nvi))
      allocate (s(0:2*nvi),v(0:2*nvi))
      allocate (wvi(mvi),xvi(mvi))

!
!     vibrational quadrature rule
!
      if (smin_asym_override .gt. 0.d0) then
         smin = smin_asym_override
      else if (trim(PESMODE) .eq. 'KRb2' .or. trim(PESMODE) .eq. 'KRB2') then
         smin = 1.0d0
      else
         smin = 0.d0
      endif
      smin_eff = smin
      call qvib (smin,smax,mvi,wvi,xvi)
!
!     kinetic energy
!
      pi = acos(-1.d0)
      tscale = pi/smax
      do jvi = 1,nvi
         t(jvi) = (jvi*tscale)**2
      enddo
      do jvi = 0,2*nvi
         s(jvi) = 0.d0
      enddo
      do kvi = 1,mvi
         sa = xvi(kvi)
         weight = wvi(kvi)/smax
         onsa2 = weight/sa**2
         arg = tscale*sa
         do jvi = 0,2*nvi
            cosine = cos(jvi*arg)
            s(jvi) = s(jvi)+cosine*onsa2
         enddo
      enddo
!
!     loop over arrangements
!
      n = 0
      ered = rmlmda*emax
      na = 3-iabs(jpar)
       do ia = 1,na
         vmin_ref = 1.d300
!
!        reference potential
!
          do jvi = 0,2*nvi
             v(jvi) = 0.d0
         enddo
          do kvi = 1,mvi
             sa = xvi(kvi)
             sa_ref = max(sa,smin_eff)
             weight = wvi(kvi)/smax
             call potenl (100.d0,sa_ref,0.d0,va,ia)
             vraw = va
             vmin_ref = min(vmin_ref,vraw)
             va = weight*va
             arg = tscale*sa
             do jvi = 0,2*nvi
               cosine = cos(jvi*arg)
               v(jvi) = v(jvi)+cosine*va
            enddo
         enddo
!
!        loop over rotational quantum numbers
!
         if (ia .eq. 1) then
            jmin = (1-jpar)/2
            jinc = 1+iabs(jpar)
         else
            jmin = 0
            jinc = 1
         endif
          do j = jmin,jmax,jinc
             kmaxt = min(j,jtot)
             kmaxj = min(j,kmax)
             apply_floor = (trim(PESMODE) .eq. 'KRb2' .or. trim(PESMODE) .eq. 'KRB2')
             if (apply_floor) then
                e_floor = vmin_ref - KRB2_SPURIOUS_MARGIN_EV*rmlmda
             else
                e_floor = -1.d300
             endif
!
!           vibrational eigenvalue problem
!
            ctot = jtot*(jtot+1)
            cent = j*(j+1)
            do jvi = 1,nvi
               do ivi = 1,jvi
                  h(ivi,jvi) = cent*(s(jvi-ivi)-s(jvi+ivi)) &
                             +      (v(jvi-ivi)-v(jvi+ivi))
                  h(jvi,ivi) = h(ivi,jvi)
               enddo
               h(jvi,jvi) = h(jvi,jvi)+t(jvi)
            enddo
            call symevp (h,nvi,nvi,e,ierr)
            if (ierr .ne. 0) stop 'basis 1'
!
!           basis functions with eint .le. ered
!
             nvphys = -1
             do kvi = 1,nvi
                if (e(kvi) .ge. e_floor .and. e(kvi) .le. ered) then
                   nvphys = nvphys+1
                   do k = kmin,kmaxj
                      n = n+1
                      if (n .le. nmax) then
!
!                       vibrational expansion coefficients
!
                        do jvi = 1,nvi
                           cvi(jvi,n) = h(jvi,kvi)
                        enddo
!
!                       Coriolis coupling matrix elements
!
                        if (k .eq. kmin) then
                           cro(1,n) = 0.d0
                        else if (k .eq. 1) then
                           cro(1,n) = -sqrt(2.d0*ctot*cent)
                        else
                           ckm = k*(k-1)
                           cro(1,n) = -sqrt((ctot-ckm)*(cent-ckm))
                        endif
                        cro(2,n) = ctot+cent-2*k*k
                        if (k .eq. kmaxt) then
                           cro(3,n) = 0.d0
                        else if (k .eq. 0) then
                           cro(3,n) = -sqrt(2.d0*ctot*cent)
                        else
                           ckp = k*(k+1)
                           cro(3,n) = -sqrt((ctot-ckp)*(cent-ckp))
                        endif
!
!                       quantum numbers
!
                         eint(n) = e(kvi)
                         ilev(n) = ia
                         jlev(n) = j
                         klev(n) = k
                         nlev(n) = nvphys
                      endif
                   enddo
                endif
            enddo
         enddo
      enddo
      if (n .gt. nmax) then
         write (6,99) n,nmax
  99     format(/1x,'Error: n = ',i6,' > nmax = ',i6/1x, &
                    'Run halted in subroutine basis')
         stop
      endif
!
!     channel ordering
!
      call basort (cvi,nvi,cro,eint,ilev,jlev,klev,nlev,n)
!
!     phase convention
!
      call bclean (cvi,nvi,n,smin,smax,mvi)
!
!     basis summary
!
      open(66,file='./output/energy_diatomic.dat')
      write (6,61)
      write (66,61)
  61  format(/1x,'BASIS:'/1x,70('-')/1x, &
      ' Channel       a       v       j       k       Energy (eV) '/1x, &
      70('-'))
      do i = 1,n
         temp = eint(i)/rmlmda
         write (6,62) i,ilev(i),nlev(i),jlev(i),klev(i),temp
         write (66,62) i,ilev(i),nlev(i),jlev(i),klev(i),temp
  62     format(1x,5i8,1x,es24.16)
      enddo
      write (6,63)
  63  format(1x,70('-'))
      return
      end
end module mod_basis 
