!     File created at Fri Jun  5 21:58:58 PDT 2020
!     Original source code: solve.f
module mod_solve 
      private 
      public :: solve 
      contains 
      subroutine solve (cvi,cro,ilev,jlev,klev,nlev,ered,nnrg)
            use mod_common, only: mro,mvi,nvi,n,nchan,rmin,rmax,smax,mtr,rmlmda,rdat
            use mod_scratchfile, only: getrec,putrec
            use mod_logder, only: logder 
            use mod_sector, only: sector 
            use mod_sbasis, only: sbasis
            use mod_overlp, only: overlp
            use mod_metric, only: metric
            use mod_pes_multi, only: NPES
       implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine solves the hyperspherical coordinate
!     coupled-channel equations between rho = rmin and rmax.
!     -----------------------------------------------------------------
!
!     common blocks
!
      ! COMMON blocks migrated to mod_common
!
!     input arrays
!
      dimension cvi(nvi,n),cro(3,n)
      dimension ilev(n),jlev(n),klev(n),nlev(n)
      dimension ered(nnrg)
!
!     local arrays
!
!     dimension cvp(nvi,n)
      allocatable cvp(:,:)

!     dimension eint(n),e(n)
      allocatable eint(:),e(:)

!     dimension h(n,n),s(n,n),c(n,n)
      allocatable h(:,:),s(:,:),c(:,:)

      allocate (cvp(nvi,n))
      allocate (eint(n),e(nchan))
      allocate (h(nchan,nchan),s(nchan,nchan),c(nchan,nchan))
      if (NPES .gt. 1 .and. rdat .gt. 0.d0) then
         stop 'exact DAT overlap chain not yet implemented'
      endif

!
!     loop over sectors
!
      npp = 0
      drho = (rmax-rmin)/mtr
      rhop = rmin-0.5d0*drho
      write(6,'(A,1x,I0,1x,A,1x,I0,1x,A,1x,I0)') 'SOLVE_REAL: start mtr', mtr, 'nchan', nchan, 'nnrg', nnrg
      do ktr = 1,mtr
!
!        surface eigenvalue problem
!
         rho = rhop+drho
         if (mod(ktr-1,10).eq.0) write(6,'(A,1x,I0,1x,A,1x,I0,1x,A,1x,ES12.4)') 'SOLVE_REAL: sector', ktr, 'of', mtr, 'rho', rho
            call sector (rho,cvi,cro,ilev,jlev,klev,nlev,h,s,c,e,np)
!
!        sector-to-sector transformation matrix
!
         if (npp .eq. 0) then
            call putrec (c,nchan,0)
         else
            call metric (rhop,cvp,npp,rho,cvi,np,ilev,jlev,klev,h,s,c)
          endif

!
!        log derivative propagation
!
          call logder (drho,ered,nnrg,h,s,c,e,nchan,np,npp)
!
!        setup for the next sector
!
         do j = 1,n
            do k = 1,nvi
               cvp(k,j) = cvi(k,j)
            enddo
         enddo
         npp = np
         rhop = rho
      enddo
!
!     transformation to final coupled-channel basis set at rmax
!
      npp = nchan+1
      rho = rmax
      call sbasis (rho,cvi,e,ilev,jlev,klev,nlev,1)
      call overlp (rhop,cvp,rho,cvi,ilev,jlev,klev,h)
!     call overlp (rho,cvi,rhop,cvp,ilev,jlev,klev,s)
!     do j = 1,n
!     do i = 1,n
!     h(i,j) = 0.5d0*(h(i,j)+s(j,i))
!     enddo
!     enddo
      call getrec (s,nchan,0)
      call dgemm  ('t','n',nchan,np,nchan,1.d0,h,nchan,s,nchan,0.d0,c,nchan)
      call logder (drho,ered,nnrg,h,s,c,e,nchan,np,npp)
      return
      end
end module mod_solve 
