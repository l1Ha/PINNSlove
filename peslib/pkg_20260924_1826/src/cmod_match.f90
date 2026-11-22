module mod_match
      private
      public:: match
      contains
      subroutine match (cvr,cvi,cro,eint,ilev,jlev,klev,nlev,ered,nnrg)
          use mod_common, only: rmlmda, nout, jout, jtot, ipar, jpar, jmax, kmin, kmax, nchan, nvi, n, ireact
          use mod_pes_multi, only: NPES, asymptotic_npes
          use Math_Functions,only:bessel 
          use mod_scratchfile,only:getrec 
          use mod_frames,only:frames
          use mod_smatrx,only:smatrx
          use mod_asymptotic_real, only: asymptotic_real_diagnostics
          use mod_match_real_diag, only: build_k_matrix_real

      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine uses the asymptotic log derivative matrix
!     at each energy to calculate and output the scattering matrix.
!     -----------------------------------------------------------------
!
!     common blocks
!
   ! COMMON blocks and local `scales` replaced by `mod_common`
!
!     input arrays
!
      dimension cvr(nvi,n),cvi(nvi,n),cro(3,n)
      dimension ilev(n),jlev(n),klev(n),nlev(n)
      dimension eint(n),ered(nnrg)
!     
      !add label for energy loop
      character(5):: label_e 
!     local arrays
!
!     dimension x(n,n),y(n,n),z(n,n)
      double precision, allocatable :: x(:,:), y(:,:), z(:,:)
      double precision, allocatable :: yfull(:,:), zfull(:,:)

!     dimension llev(n),crp(3,nnrg)
      double precision, allocatable :: llev(:), crp(:,:), eintc(:)
      double precision, allocatable :: phase(:), qdiag(:), qeval(:)
      double precision, allocatable :: kprev(:,:), knext(:,:), dkde(:,:)
      integer, allocatable :: ilevc(:), jlevc(:), klevc(:), nlevc(:)
      integer :: pp, qq
      integer :: nchan_full, npes_full, nchan_work, npes_work
      integer :: old_nchan, old_npes, icfull, jcfull
      logical :: reduced_match
      integer :: nopen_diag, idiag
      double precision :: den

      nchan_full = nchan
      npes_full = max(1,NPES)
      npes_work = asymptotic_npes()
      nchan_work = n*npes_work
      reduced_match = nchan_work .lt. nchan_full

      allocate (x(nchan_work,nchan_work),y(nchan_work,nchan_work),z(nchan_work,nchan_work))
      allocate (llev(nchan_work),crp(3,nnrg))
      allocate (eintc(nchan_work),ilevc(nchan_work),jlevc(nchan_work),klevc(nchan_work),nlevc(nchan_work))
      allocate (phase(nchan_work), qdiag(nchan_work), qeval(nchan_work))
      allocate (kprev(nchan_work,nchan_work), knext(nchan_work,nchan_work), dkde(nchan_work,nchan_work))
      if (reduced_match) then
         allocate (yfull(nchan_full,nchan_full),zfull(nchan_full,nchan_full))
      endif

!
!     loop over energies
!     
      open(667,file='./output/energy.dat')
      open(670,file='./output/Eigenphase-real.dat')
      open(671,file='./output/Lifetime-real.dat')
      open(672,file='./output/Qdiag-real.dat')
      do 3 nrg = 1,nnrg
!         
         write(label_e,'(I5)') nrg
!        BF to SF transformation of Y
!
         if (reduced_match) then
            call getrec (yfull,nchan_full,nrg)
            do jc = 1,n
               jcfull = (jc-1)*npes_full + 1
               do ic = 1,n
                  icfull = (ic-1)*npes_full + 1
                  y(ic,jc) = yfull(icfull,jcfull)
               enddo
            enddo
            old_nchan = nchan
            old_npes = NPES
            nchan = nchan_work
            NPES = npes_work
         else
            call getrec (y,nchan,nrg)
         endif
          call frames (y,z,cro,jlev,klev,llev,+1)

!        build composite-channel maps (replicate base channel quantum numbers/thresholds)
         do ic = 1,nchan
            ibase = (ic-1)/max(1,NPES) + 1
            eintc(ic) = eint(ibase)
            ilevc(ic) = ilev(ibase)
            jlevc(ic) = jlev(ibase)
            klevc(ic) = klev(ibase)
            nlevc(ic) = nlev(ibase)
         enddo
!
!        asymptotic matching
!
          enrg = ered(nrg)
          call bessel (x,y,z,cvr,cvi,eint,ilev,jlev,klev,llev,nlev,enrg)
          phase = 0.d0
          qdiag = 0.d0
          qeval = 0.d0
          dkde = 0.d0
          if (nnrg .gt. 1) then
             if (nrg .eq. 1) then
                call build_k_matrix_real(2, knext, cvr, cvi, cro, eint, ilev, &
                                         jlev, klev, nlev, ered)
                den = (ered(2) - ered(1)) / rmlmda
                if (abs(den) .gt. 0.d0) dkde = (knext - x) / den
             else if (nrg .eq. nnrg) then
                call build_k_matrix_real(nnrg-1, kprev, cvr, cvi, cro, eint, &
                                         ilev, jlev, klev, nlev, ered)
                den = (ered(nnrg) - ered(nnrg-1)) / rmlmda
                if (abs(den) .gt. 0.d0) dkde = (x - kprev) / den
             else
                call build_k_matrix_real(nrg-1, kprev, cvr, cvi, cro, eint, &
                                         ilev, jlev, klev, nlev, ered)
                call build_k_matrix_real(nrg+1, knext, cvr, cvi, cro, eint, &
                                         ilev, jlev, klev, nlev, ered)
                den = (ered(nrg+1) - ered(nrg-1)) / rmlmda
                if (abs(den) .gt. 0.d0) dkde = (knext - kprev) / den
             endif
          endif
          call asymptotic_real_diagnostics(x, dkde, eintc, nchan, enrg, &
                                           nnrg .gt. 1, phase, qdiag, qeval, &
                                           nopen_diag, idiag)
          if (idiag .ne. 0) stop 'match 0'
          write(670,'(*(g0,1x))') enrg/rmlmda, nopen_diag, &
                                  (phase(i), i=1,nopen_diag)
          if (nnrg .gt. 1) then
             write(671,'(*(g0,1x))') enrg/rmlmda, nopen_diag, &
                                     (qeval(i), i=1,nopen_diag)
             write(672,'(*(g0,1x))') enrg/rmlmda, nopen_diag, &
                                     (qdiag(i), i=1,nopen_diag)
          endif
          call smatrx (x,y,z,eintc,llev,nchan,seps,ierr,enrg)
          if (ierr .ne. 0) stop 'match 1'
         enrg = enrg/rmlmda
!
!        SF to BF transformation of S
!
         call frames (x,z,cro,jlev,klev,llev,-1)
         call frames (y,z,cro,jlev,klev,llev,-1)
!
!        phase convention for parity-adapted BF S-matrix
!
          if (ipar .eq. -1) then
             do j = 1,nchan
                do i = 1,nchan
                   x(i,j) = -x(i,j)
                   y(i,j) = -y(i,j)
                enddo
             enddo
          endif
          if (max(1,NPES) .gt. 1) then
             do ibase = 1,n
                if (ilev(ibase) .eq. ireact .and. nlev(ibase) .eq. 0 .and. jlev(ibase) .eq. 0 .and. klev(ibase) .eq. 0) then
                   open(668,file='./output/S-block'//trim(adjustl(label_e))//'.dat')
                   do pp = 1,max(1,NPES)
                      ic = (ibase-1)*max(1,NPES) + pp
                      write(668,'(*(g0,1x))') (x(ic,(ibase-1)*max(1,NPES)+qq), y(ic,(ibase-1)*max(1,NPES)+qq), qq=1,max(1,NPES))
                   enddo
                   close(668)
                   exit
                endif
             enddo
          endif
!        cumulative reaction probabilities
!
         do ia = 1,3
            crp(ia,nrg) = 0.d0
         enddo
         ! accumulate reaction probabilities by mapping composite channels to base channels
         do ibase = 1,n
            if (ered(nrg).gt.eint(ibase) .and. ilev(ibase).eq.ireact) then
               do jbase = 1,n
                  if (ered(nrg).gt.eint(jbase) .and. ilev(jbase).ne.ireact) then
                      tp_sum = 0.d0
                      do pp = 1, max(1,NPES)
                         ic = (ibase-1)*max(1,NPES) + pp
                         do qq = 1, max(1,NPES)
                            jc = (jbase-1)*max(1,NPES) + qq
                            tp_sum = tp_sum + x(ic,jc)**2 + y(ic,jc)**2
                         enddo
                      enddo
                      crp(ilev(jbase),nrg) = crp(ilev(jbase),nrg) + tp_sum
                  endif
               enddo
            endif
         enddo
!
!        OUTPUT:
!        This outputs S-matrix elements and transition probabilities
!        from all open A+BC(v,j) reactant channels with v.le.nout
!        and j.le.jout to all open A+BC, B+CA, C+AB product channels:
!        
         open(66,file='./output/S-matrix'//trim(adjustl(label_e))//'.dat')
         open(666,file='./output/S-matrixsum'//trim(adjustl(label_e))//'.dat')
         open(669,file='./output/S-open'//trim(adjustl(label_e))//'.dat')
         
         do ibase = 1,n
            if (ered(nrg) .le. eint(ibase)) cycle
            ii = ilev(ibase)
            ji = jlev(ibase)
            ki = klev(ibase)
            ni = nlev(ibase)
            do jbase = 1,n
               if (ered(nrg) .le. eint(jbase)) cycle
               iif = ilev(jbase)
               jf = jlev(jbase)
               kf = klev(jbase)
               nf = nlev(jbase)
               tp = 0.d0
               sr = 0.d0
               si = 0.d0
               do pp = 1, max(1,NPES)
                  ic = (ibase-1)*max(1,NPES) + pp
                  do qq = 1, max(1,NPES)
                     jc = (jbase-1)*max(1,NPES) + qq
                     sr = sr + x(ic,jc)
                     si = si + y(ic,jc)
                     tp = tp + x(ic,jc)**2 + y(ic,jc)**2
                  enddo
               enddo
               write (669,'(*(g0,1x))') ii,ni,ji,ki,iif,nf,jf,kf,sr,si,tp
            enddo
         enddo

         unit = 0.d0
         cols = 0.d0
         write (6,61) enrg,seps
         write (667,'(*(g0,1x))') nrg,enrg,seps
         do 2 ibase = 1,n
            if (ered(nrg) .le. eint(ibase)) go to 2
            ii = ilev(ibase)
            ji = jlev(ibase)
            ki = klev(ibase)
            ni = nlev(ibase)
            if (ii .ne. ireact) go to 2
            if (ni .gt. nout) go to 2
            if (ji .gt. jout) go to 2
            jp = -1
            sumj = 0.d0
            cols = cols+1.d0
            do 1 jbase = 1,n
               if (ered(nrg) .le. eint(jbase)) go to 1
               iif = ilev(jbase)
               jf = jlev(jbase)
               kf = klev(jbase)
               nf = nlev(jbase)
               ! aggregate over PES surfaces
                 tp = 0.d0
                 sr = 0.d0
                 si = 0.d0
                 do pp = 1, max(1,NPES)
                    ic = (ibase-1)*max(1,NPES) + pp
                    do qq = 1, max(1,NPES)
                       jc = (jbase-1)*max(1,NPES) + qq
                       sr = sr + x(ic,jc)
                       si = si + y(ic,jc)
                       tp = tp + x(ic,jc)**2 + y(ic,jc)**2
                    enddo
                 enddo
               unit = unit+tp
               if (jf .lt. jp) then
                  write (6,63) ii,ni,ji,ki,ip,np,sumj
                  write (666,'(*(g0,1x))') ii,ni,ji,ki,ip,np,sumj
                  sumj = 0.d0
               endif
               ip = iif
               np = nf
               jp = jf
               sumj = sumj+tp
                write (6,62) ii,ni,ji,ki,iif,nf,jf,kf,sr,si,tp
                write (66,'(*(g0,1x))') ii,ni,ji,ki,iif,nf,jf,kf,sr,si,tp
   1        continue
            write (6,63) ii,ni,ji,ki,ip,np,sumj
            write (666,'(*(g0,1x))') ii,ni,ji,ki,ip,np,sumj
   2     continue
         if (cols .gt. 0.d0) then
            unit = unit/cols
            write (6,64) unit
            ! write (66,'(*(g0,1x))') unit
         endif
         close(66) 
         close(666)
         close(669)
   
   3  continue
      close(667)
      close(670)
      close(671)
      close(672)
!      
!     And this outputs cumulative reaction probabilities:
!     
      open(66,file='./output/ReactionProbabilities.dat')
      if (jpar .eq. 0) then
         write (6,65)
         do nrg = 1,nnrg
            enrg = ered(nrg)/rmlmda
            write (6,66) enrg,crp(2,nrg),crp(3,nrg)
            write (66,'(*(g0,1x))') enrg,crp(2,nrg),crp(3,nrg)
         enddo
         write (6,67)
      else
         write (6,68)
         
         do nrg = 1,nnrg
            enrg = ered(nrg)/rmlmda
            write (6,69) enrg,crp(2,nrg)
            write (66,'(*(g0,1x))') enrg,crp(2,nrg)
         enddo
         write (6,67)
      endif
      close(66)
      if (reduced_match) then
         nchan = nchan_full
         NPES = npes_full
      endif
      return

  61  format(/1x,'MATCH:'/1x,70('-')/1x, &
       ' E(eV) = ',es20.12,24x,'Seps = ',1p,e8.1,0p/1x,70('-')/1x, &
       '  a    v    j    k    a''   v''   j''   k''     Re(S)', &
       '     Im(S)     |S|^2 '/1x,70('-'))
  62  format(1x,i3,7i5,1x,3f10.5)
  63  format(1x,i3,5i5,5x,'  all ',f30.5)
  64  format(1x,70('-')/1x,29x,'Unitarity',f31.5/1x,70('-'))
  65  format(/1x,'Cumulative reaction probabilities:'/1x,70('-')/1x, &
       '       E(eV)          A+BC --> B+CA          A+BC --> C+AB ' &
       /1x,70('-'))
  66  format(1x,f13.5,f19.5,f23.5)
  67  format(1x,70('-'))
  68  format(/1x,'Cumulative reaction probabilities:'/1x,70('-')/1x, &
      '       E(eV)                                     A+B2 --> B+AB ' &
       /1x,70('-'))
  69  format(1x,f13.5,f46.5)
      end
   end module mod_match 
