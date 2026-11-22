!     File created at Fri Jun  5 21:58:57 PDT 2020
!     Original source code: metric.f
module mod_metric
      private
      public :: metric
      contains

      subroutine metric (rhop,cvp,npp,rho,cvi,np,ilev,jlev,klev,h,s,c)
      use mod_common
      use mod_scratchfile, only: putrec,getrec 
      use mod_overlp, only: overlp 
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine constructs the sector-to-sector transformation
!     matrix t between sectors centred at hyperradii rhop and rho.
!     Compare with subroutine sector.
!     -----------------------------------------------------------------
!
!     common blocks
!
      ! COMMON /arrays/ replaced by `mod_common`
!
!     input arrays
!
       dimension cvp(nvi,n),cvi(nvi,n)
       dimension ilev(n),jlev(n),klev(n)
       dimension h(nchan,nchan),s(nchan,nchan),c(nchan,nchan)
!
!     primitive overlap integrals
!
      call overlp (rhop,cvp,rho,cvi,ilev,jlev,klev,h)
!     call overlp (rho,cvi,rhop,cvp,ilev,jlev,klev,s)
!     seps = 0.d0
!     do j = 1,n
!     do i = 1,n
!     seps = max(seps,abs(h(i,j)-s(j,i)))
!     h(i,j) = 0.5d0*(h(i,j)+s(j,i))
!     enddo
!     enddo
!
!     sector-to-sector transformation matrix
!
      call dgemm ('n','n',nchan,np,nchan,1.d0,h,nchan,c,nchan,0.d0,s,nchan)
      call getrec (h,nchan,0)
      call putrec (c,nchan,0)
      call dgemm ('t','n',npp,np,nchan,1.d0,h,nchan,s,nchan,0.d0,c,nchan)
!
!     metric summary
!
!     write (6,61) rhop,rho,seps
!     61 format(1x,'METRIC:',
!     +' rhop = ',f8.3,
!     +' rho = ',f8.3,
!     +' seps = ',1p,e8.1)
      return
      end


      end module mod_metric 
      
