!     File created at Fri Jun  5 21:58:57 PDT 2020
!     Original source code: sector.f
module mod_sector
      private
      public :: sector

      contains

      subroutine sector (rho,cvi,cro,ilev,jlev,klev,nlev,h,s,c,e,np)
            use mod_sbasis, only: sbasis
            use mod_couple, only: couple
            use mod_hcsevp, only: hcsevp
            use mod_common
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine sets up and solves the hyperspherical
!     coordinate surface eigenvalue problem at hyperradius rho.
!     -----------------------------------------------------------------
!
!     common blocks (migrated to mod_common)
!
!     input arrays
!
      dimension cvi(nvi,n),cro(3,n)
      dimension ilev(n),jlev(n),klev(n),nlev(n)
      dimension h(nchan,nchan),s(nchan,nchan),c(nchan,nchan),e(nchan)

!     local arrays
!
!     dimension eint(n),sigma(nchan)
      allocatable eint(:),sigma(:)

      allocate (eint(n),sigma(nchan))

!
!     surface basis set
!
      write(6,'(A,1x,ES12.4,1x,A,1x,I0)') 'SECTOR: sbasis rho', rho, 'nlev', nlev(1)
      call sbasis (rho,cvi,eint,ilev,jlev,klev,nlev,0)
!
!     hamiltonian and overlap matrices
!
      write(6,'(A,1x,ES12.4)') 'SECTOR: couple rho', rho
      call couple (rho,cvi,cro,eint,ilev,jlev,klev,h,s)
!
!     surface eigenvalue problem
!
      write(6,'(A,1x,I0)') 'SECTOR: hcsevp nchan', nchan
      call hcsevp (h,s,c,nchan,e,sigma,np)
!
!     sector summary
!
      write (6,61) e(1)/rmlmda,e(np)/rmlmda,np

  61  format(1x,'SECTOR:', &
      '  e(1) = ',f8.3, &
      '  e(np) = ',f8.3, &
      '    np = ',i8)
      return
      end


      end module mod_sector 