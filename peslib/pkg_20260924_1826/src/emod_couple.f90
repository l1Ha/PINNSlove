!     File created at Fri Jun  5 21:58:55 PDT 2020
!     Original source code: couple.f
module mod_couple 
      private
      public :: couple
      contains 

      subroutine couple (rho,cvi,cro,eint,ilev,jlev,klev,h,s)
         use mod_common
         use mod_direct, only: direct
         use mod_exchng, only: exchng
         use mod_pes_multi, only: NPES, active_npes
      
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine constructs the Hamiltonian and overlap matrices
!     h and s that are needed to solve the hyperspherical coordinate
!     surface eigenvalue problem at the given value of rho.
!     -----------------------------------------------------------------
!
!     common blocks (migrated to mod_common)
!
!     input arrays
!
      dimension cvi(nvi,n),cro(3,n)
      dimension eint(n),ilev(n),jlev(n),klev(n)
      dimension h(nchan,nchan),s(nchan,nchan)
!
!     local arrays
!
!     dimension v(mro,mvi)
      allocatable v(:,:)

      integer ibase, i, j, ia, ib, na, pp, ic, npes_max, npes_act
      character(len=8) :: dbg_couple
      double precision :: h23max, s23max

       allocate (v(mro,mvi))
       v = 0.d0

       npes_max = max(1,NPES)
       npes_act = active_npes(rho)

!
      do j = 1,nchan
         do i = 1,nchan
            h(i,j) = 0.d0
            s(i,j) = 0.d0
         enddo
         s(j,j) = 1.d0
      enddo
      na = 3-iabs(jpar)
      do ia = 1,na
         call direct (rho,cvi,cro,ilev,jlev,klev,h,v,ia)
         do ib = 1,na
            if (ib .ne. ia) then
               call exchng (rho,cvi,cro,ilev,jlev,klev,h,s,v,ia,ib)
            endif
         enddo
      enddo
       if (na .eq. 2) then
          ia = 2
          ib = 3
          call exchng (rho,cvi,cro,ilev,jlev,klev,h,s,v,ia,ib)
       endif
       if (npes_act .lt. npes_max) then
          do ibase = 1,n
             do pp = npes_act+1,npes_max
                ic = (ibase-1)*npes_max + pp
                do j = 1,nchan
                   h(ic,j) = 0.d0
                   h(j,ic) = 0.d0
                   s(ic,j) = 0.d0
                   s(j,ic) = 0.d0
                enddo
             enddo
          enddo
       endif
       call get_environment_variable('ABC_DEBUG_COUPLE', dbg_couple)
       if (len_trim(dbg_couple) .gt. 0) then
          h23max = 0.d0
          s23max = 0.d0
          do j = 1,nchan
             ibase = (j-1)/npes_max + 1
             if (ilev(ibase) .ne. 3) cycle
             do i = 1,nchan
                ibase = (i-1)/npes_max + 1
                if (ilev(ibase) .ne. 2) cycle
                h23max = max(h23max, abs(h(i,j)))
                s23max = max(s23max, abs(s(i,j)))
             enddo
          enddo
          write (6,'(1x,A,F8.3,1x,A,1P,E12.4,1x,A,E12.4)') &
             'COUPLE_DEBUG: rho = ', rho, 'max|H23| = ', h23max, &
             'max|S23| = ', s23max
       endif
       do j = 1,nchan
          do i = 1,nchan
             ibase = (i-1)/npes_max + 1
             h(i,j) = h(i,j) + eint(ibase)*s(i,j)
          enddo
       enddo
      heps = 0.d0
      seps = 0.d0
      do j = 2,nchan
         do i = 1,j-1
            heps = max(heps,abs(h(i,j)-h(j,i)))
            h(i,j) = 0.5d0*(h(i,j)+h(j,i))
            h(j,i) = h(i,j)
            seps = max(seps,abs(s(i,j)-s(j,i)))
            s(i,j) = 0.5d0*(s(i,j)+s(j,i))
            s(j,i) = s(i,j)
         enddo
      enddo
!
!     couple summary
!
      write (6,61) rho,heps/rmlmda,seps
  61  format(/1x,'COUPLE:', &
      '   rho = ',f8.3, &
      '   heps = ',f8.3, &
      '  seps = ',1p,e8.1)
      return
      end

   end module mod_couple 
