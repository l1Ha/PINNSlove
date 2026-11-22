!     File created at Fri Jun  5 21:58:57 PDT 2020
!     Original source code: potenl.f
module mod_Potenl
   use mod_pes_multi, only: NPES, PESMODE, potenl_ref
   private
   public :: potenl,writePES 
   contains

      subroutine writePES()
      use mod_mkmp2,only:potface1
      use mod_common
         integer,parameter:: rk=kind(1.0d0) 
         integer:: i,j,k
         integer:: nr2,nr3,na23 
         real(rk):: r2min,r2max,r3min,r3max 
         real(rk):: r2,r3,a23,am3(3),dx 
         real(rk),allocatable:: v3(:,:,:),r2s(:),r3s(:),a23s(:) 
         print*,'mass',mass 
         r2min = 0.5_rk
         r2max = 5.0_rk
         r3min = 0.5_rk
         r3max = 15.0_rk
         am3(1) = mass(1)
            am3(2) = mass(2)
            am3(3) = mass(3) 
         
         nr2 = 200
         nr3 = 300

         na23 = int(3.151592653589793_rk/0.05_rk) 
         allocate (r2s(nr2),r3s(nr3),a23s(na23))
         allocate (v3(nr3,nr2,na23)) 
         do k=1,na23 
         do j=1,nr2 
         do i=1,nr3 
            r3s(i) = r3min+i*(r3max-r3min)/(nr3+1.0_rk)
            r2s(j) = r2min+j*(r2max-r2min)/(nr2+1.0_rk)
            a23s(k) = 0.0_rk+(k-1)*0.05_rk
            ! print*,'r2,r3,a23',r3s(i),r2s(j),a23s(k)

            call potenl(r3s(i),r2s(j),cos(a23s(k)),v3(i,j,k),1)
            ! print*,r3s(i),r2s(j),cos(a23s(k)),v3(i,j,k)
            ! if (sum(mass) < 3.5) then
            !  call potface1(1,1,am3,r3s(i),r2s(j),a23s(k),v3(i,j,k))
            ! else if (sum(mass) > 3.5 .and. sum(mass) < 4.5) then 
            !   call potenl(r3s(i),r2s(j),cos(a23s(k)),v3(i,j,k),1)
            ! endif
            ! call potface1(1,1,am3,r3s(i),r2s(j),a23s(k),v3(i,j,k))

      enddo 
  
      enddo 
      enddo 

      
         open(unit=66,file='./output/potxy.dat')
         write(66,'(*(g0,1x))') 0.0_rk,(r3s(i),i=1,nr3)
         do j = 1,nr2
            write(66,'(*(g0,1x))') r2s(j),(minval(v3(i,j,:)),i=1,nr3) 
         enddo
         close(66)  

         open(unit=66,file='./output/potxz.dat')
         write(66,'(*(g0,1x))') 0.0_rk,(r3s(i),i=1,nr3)
         do k=1,na23 
         write(66,'(*(g0,1x))') a23s(k),(minval(v3(i,:,k)),i=1,nr3)
         enddo 
         close(66) 
   

         open(unit=66,file='./output/potyz.dat')
         write(66,'(*(g0,1x))') 0.0_rk,(r2s(j),j=1,nr2)
         do k=1,na23
         write(66,'(*(g0,1x))') a23s(k),(minval(v3(:,j,k)),j=1,nr2)
         enddo 
         close(66) 
         
         deallocate(v3,r2s,r3s,a23s) 
      end subroutine writePES 
            
      
      
      subroutine pes (r12,r23,r31,vev)
            use mod_mkmp2,only:pot0 
            ! use mod_pot_HeStar_H2,only:pot_HeStar_H2_rrr
            implicit double precision (a-h,o-z)
      !
      !     ----------------------------------------------
      !     User defined PES should be written here in eV.
      !     ----------------------------------------------
      !  
            real(8):: x(3),v(1)
            x(1) = r12
            x(2) = r23
            x(3) = r31
             call pot0(1,v,x) 
            ! v(1) = pot_HeStar_H2_rrr(x(1),x(2),x(3))            

            vev = v(1)*27.211385d0 ! Your surface goes here.
            return
      end 

            
      subroutine potenl (ra,sa,cosa,va,ia)
       !use amod_pot_HeStar_H2,only: pot_HeStar_H2
      use mod_mkmp2,only:potface1
      use mod_common
     
      ! use apot,only: potsub1

      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine returns the potential energy surface va as a
!     function of mass-scaled Jacobi coordinates in arrangement ia.
!     -----------------------------------------------------------------
!
      dimension r(3)
      dimension am3(3)
      real(8) :: vref

      if (NPES .gt. 1 .or. trim(PESMODE) .eq. 'KRb2' .or. trim(PESMODE) .eq. 'KRB2' .or. trim(PESMODE) .eq. 'H3') then
          call potenl_ref(ra,sa,cosa,vref,ia)
          va = rmlmda * vref
          return
      endif
! 
      ib = ia+1
      if (ib .gt. 3) ib = 1
      ic = 6-ia-ib
      rap = ra/scales(ia)
      cap = 2.0d0*cosa*rap
      sap = scales(ia)*sa
      sbp = mass(ib)/(mass(ib)+mass(ic))*sap
      scp = mass(ic)/(mass(ib)+mass(ic))*sap
      r(ia) = sap
      r(ib) = sqrt(rap**2-cap*sbp+sbp**2)
      r(ic) = sqrt(rap**2+cap*scp+scp**2)
      call potsub (r,vev)
      va = rmlmda*vev

      ! print*,'PES HHH'
      ! am3(1) = 1.0d0 
      ! am3(2) = 1.0d0
      ! am3(3) = 1.0d0
      ! call potface1(1,1,am3,ra,sa,acos(cosa),va)
      ! va = va*27.211385d0 
      ! va = pot_HeStar_H2(ra,sa,acos(cosa))
      return
      end
      !     File created at Fri Jun  5 21:58:57 PDT 2020
!     Original source code: potsub.f

      subroutine potsub (r,vev)
      use mod_common
      ! use mod_pot_HeStar_H2,only:pot_HeStar_H2_rrr
         implicit double precision (a-h,o-z)
      !
      !     -----------------------------------------------------------------
      !     This subroutine chooses which potential to use on the basis
      !     of the atomic masses in common /masses/, and also ensures
      !     that the potential is called with the bond lengths in the
      !     correct order.
      !     -----------------------------------------------------------------
      !
            dimension r(3),s(3),m(3)
      !
            imax = 1
            imin = 1
            do i = 1,3
               m(i) = nint(mass(i))
               if (m(i) .gt. m(imax)) then
                  imax = i
               else if (m(i) .lt. m(imin)) then
                  imin = i
               endif
            enddo
            if (imax .eq. imin) then
               imin = 1
               imid = 2
               imax = 3
            else
               imid = 6-imax-imin
            endif
            s(1) = r(imax)
            s(2) = r(imid)
            s(3) = r(imin)
            mmax = m(imax)
            mmid = m(imid)
            mmin = m(imin)
          if (mmin .lt. 1) stop 'potsub 1'
          if (mmid .gt. 2) stop 'potsub 2'
          if (mmax .le. 2) then
            !  call hh2pot (s,vev)
            call pes(s(1),s(2),s(3),vev)
          else if (mmax .eq. 4)  then
            ! vev = pot_HeStar_H2_rrr(s(1),s(2),s(3))            
            !  call hh3pot (s,vev)
          else if (mmax .eq. 19) then
            !  call fh2pot (s,vev)
          else if (mmax .eq. 35) then
            !  call clh2pt (s,vev)
          else if (mmax .eq. 37) then
            !  call clh2pt (s,vev)
          else
             stop 'potsub 3'
          endif
            ! call pes(s(1),s(2),s(3),vev)
            return
            end
      

end module mod_Potenl 
