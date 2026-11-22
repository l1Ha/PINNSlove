!     File created at Fri Jun  5 21:58:56 PDT 2020
!     Original source code: driver.f
module mod_driver
      private
      public:: driver
      contains
      subroutine driver (nmax)
            use mod_common
            use mod_match,only:match 
            use mod_match_complex, only: match_complex
            use mod_solve,only:solve 
            use mod_solve_complex, only: solve_complex
            use mod_basis,only:basis 
            use mod_scratchfile,only:setrec,endrec 
            use mod_scratchfile_complex, only: setrec_complex, endrec_complex
            use mod_pes_multi, only: NPES
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     Main driving subroutine for hyperspherical coordinate calculation
!     -----------------------------------------------------------------
!
!     common blocks (migrated to mod_common)
!
!     global arrays
!
!     parameter (nmax = 1500)

!     dimension ered(nmax)
      allocatable ered(:)

!     dimension cvi(nvi,nmax),cro(3,nmax),eint(nmax)
      allocatable cvi(:,:),cro(:,:),eint(:)

!     dimension ilev(nmax),jlev(nmax),klev(nmax),nlev(nmax)
      allocatable ilev(:),jlev(:),klev(:),nlev(:)

!     dimension cvr(nvi,nmax)
      allocatable cvr(:,:)

      allocate (ered(nmax))
      allocate (cvi(nvi,nmax),cro(3,nmax),eint(nmax))
      allocate (ilev(nmax),jlev(nmax),klev(nmax),nlev(nmax))
      allocate (cvr(nvi,nmax))

!
!     rovibrational basis set
!
      call basis (cvi,cro,eint,ilev,jlev,klev,nlev,nmax)
!
!     scattering energies
!
      eminr = 1.d+30
      do j = 1,n
         if (ilev(j) .eq. ireact) then
            eminr = min(eminr,eint(j))
         endif
      enddo
      eminr = eminr/rmlmda
      emin = eminr
      write(6,'(A,ES20.12)') ' DRIVER: reactant threshold emin (eV) = ', emin
      if (iscale .eq. 1 .and. enrg .gt. emin) then
         emin = enrg
         write(6,'(A,ES20.12)') ' DRIVER: iscale=1 incident threshold override emin (eV) = ', emin
      endif
      nrg = 0
      if (iscale .eq. 1) then
!
!        --- Log-spaced collision energy mode ---
!        ecol_min, ecol_max are collision energies (eV, > 0)
!        nnrg points are log-spaced from ecol_min to ecol_max
!        total energy = emin + ecol
!
         write(6,'(A)') ' DRIVER: iscale=1  log-spaced collision energy mode'
         write(6,'(A,ES12.4,A,ES12.4,A,I6)') &
              '         ecol_min=',ecol_min,' ecol_max=',ecol_max, &
              ' nnrg=',nnrg
         if (ecol_min .le. 0.d0 .or. ecol_max .le. 0.d0) then
            write(6,*) 'FATAL: ecol_min and ecol_max must be > 0 for iscale=1'
            stop 'driver: ecol'
         endif
         alog_min = log(ecol_min)
         alog_max = log(ecol_max)
         do inrg = 1, nnrg
            if (nnrg .gt. 1) then
               frac = dble(inrg-1)/dble(nnrg-1)
            else
               frac = 0.d0
            endif
            ecol = exp(alog_min + frac*(alog_max - alog_min))
            etot = emin + ecol
            nrg = nrg + 1
            if (nrg .le. nmax) then
               ered(nrg) = rmlmda*etot
            endif
         enddo
      else
!
!        --- Linear energy mode (original ABC) ---
!
         do inrg = 1,nnrg
             if (enrg .gt. emin) then
               nrg = nrg+1
               if (nrg .le. nmax) then
                  ered(nrg) = rmlmda*enrg
               endif
             else
                write (6,61) enrg
   61          format(1x,'*** No open reactant channels at ' &
                ,f7.4,' eV. Energy skipped.')
             endif
            enrg = enrg+dnrg
         enddo
      endif
      nnrg = nrg
      if (nnrg .gt. nmax) then
         write (6,99) nnrg,nmax
   99     format(/1x,'Error: nnrg = ',i6,' > nmax = ',i6/1x, &
                    'Run halted in subroutine driver')
         stop
      endif
      if (nnrg .eq. 0) return
!
!     coupled channel propagation
!
      nchan = n * max(1, NPES)
      if (NPES .gt. 1) then
         call setrec_complex(nchan, nnrg)
         call solve_complex(cvr, cro, ilev, jlev, klev, nlev, ered, nnrg)
         call match_complex(cvr, cvi, cro, eint, ilev, jlev, klev, nlev, ered, nnrg)
         call endrec_complex
      else
         call setrec (nchan,nnrg)
         call solve (cvr,cro,ilev,jlev,klev,nlev,ered,nnrg)
         call match (cvr,cvi,cro,eint,ilev,jlev,klev,nlev,ered,nnrg)
         call endrec
      endif
      return
      end
end module mod_driver 
