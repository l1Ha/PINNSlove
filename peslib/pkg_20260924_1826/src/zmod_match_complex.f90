module mod_match_complex
  private
  public :: match_complex

contains

  subroutine match_complex(cvr, cvi, cro, eint, ilev, jlev, klev, nlev, ered, nener)
    use mod_common
    use mod_bessel_complex, only: bessel_complex
    use mod_frames_complex, only: frames_complex
    use mod_smatrx_complex, only: smatrx_complex
    use mod_scratchfile_complex, only: getrec_complex
    implicit double precision (a-h,o-z)

    dimension cvr(nvi,n), cvi(nvi,n), cro(3,n)
    dimension ilev(n), jlev(n), klev(n), nlev(n), eint(n), ered(*)

    complex*16, allocatable :: s(:,:), q(:,:), z(:,:)
    complex*16, allocatable :: kprev2(:,:), kprev1(:,:), kcurr(:,:), knext1(:,:), knext2(:,:), dk(:,:)
    double precision, allocatable :: llev(:), cmat(:,:), phase(:)
    double precision, allocatable :: qeval(:)
    integer, intent(in) :: nener
    integer :: i, j, nrg, ierr, nopen
    double precision :: seps, enow, den

    allocate(s(nchan,nchan), q(nchan,nchan), z(nchan,nchan))
    allocate(kprev2(nchan,nchan), kprev1(nchan,nchan), kcurr(nchan,nchan))
    allocate(knext1(nchan,nchan), knext2(nchan,nchan), dk(nchan,nchan))
    allocate(cmat(nchan,nchan), llev(nchan), phase(nchan), qeval(nchan))
    open(66, file='./output/S-matrix-complex.dat')
    open(67, file='./output/Eigenphase-complex.dat')
    open(68, file='./output/Lifetime-complex.dat')
    open(69, file='./output/Qdiag-complex.dat')
    open(72, file='./output/ReactionProbabilities-complex.dat')
    open(70, file='./output/TrQ-complex.dat')
    open(71, file='./output/Observables-complex.dat')
    write(70,'(A)') '# E(eV)  E_coll(eV)  T_coll(K)  nopen  TrQ(fs)  TrQ(us)'
    write(71,'(A)') '# E_coll(mK)  Lifetime_max(fs)  P_reac  k_el(cm^3/s)'

    write(6,'(A,1x,I0,1x,A,1x,I0)') 'MATCH_COMPLEX: start nener', nener, 'nchan', nchan
    flush(6)

    if (nener .ge. 1) then
      call build_k_matrix(1, kprev2, llev, cmat, z, cvr, cvi, cro, eint, &
                          ilev, jlev, klev, nlev, ered)
    endif

    if (nener .eq. 1) then
      write(6,'(A,1x,I0,1x,A,1x,I0)') 'MATCH_COMPLEX: energy', 1, 'of', nener
      flush(6)
      dk = (0.d0,0.d0)
      call write_energy_block(1, kprev2, dk, .false., s, q, qeval, cmat, &
                              llev, phase, nopen, eint, cro, ilev, jlev, klev, &
                              ered, seps)
    else if (nener .ge. 2) then
      call build_k_matrix(2, kprev1, llev, cmat, z, cvr, cvi, cro, eint, &
                          ilev, jlev, klev, nlev, ered)
      if (nener .ge. 3) call build_k_matrix(3, kcurr, llev, cmat, z, cvr, cvi, cro, &
                          eint, ilev, jlev, klev, nlev, ered)
      if (nener .ge. 4) call build_k_matrix(4, knext1, llev, cmat, z, cvr, cvi, cro, &
                          eint, ilev, jlev, klev, nlev, ered)
      if (nener .ge. 5) call build_k_matrix(5, knext2, llev, cmat, z, cvr, cvi, cro, &
                          eint, ilev, jlev, klev, nlev, ered)

      ! Energy 1: forward difference (K2 - K1) / (E2 - E1)
      write(6,'(A,1x,I0,1x,A,1x,I0)') 'MATCH_COMPLEX: energy', 1, 'of', nener
      flush(6)
      den = (ered(2) - ered(1)) / rmlmda
      if (abs(den) .gt. 0.d0) then
        dk = (kprev1 - kprev2) / den
      else
        dk = (0.d0,0.d0)
      endif
      call write_energy_block(1, kprev2, dk, .true., s, q, qeval, cmat, &
                              llev, phase, nopen, eint, cro, ilev, jlev, klev, &
                              ered, seps)

      ! Energy 2: 3-point centered (K3 - K1) / (E3 - E1), or backward if nener==2
      write(6,'(A,1x,I0,1x,A,1x,I0)') 'MATCH_COMPLEX: energy', 2, 'of', nener
      flush(6)
      if (nener .eq. 2) then
        den = (ered(2) - ered(1)) / rmlmda
        if (abs(den) .gt. 0.d0) then
          dk = (kprev1 - kprev2) / den
        else
          dk = (0.d0,0.d0)
        endif
        call write_energy_block(2, kprev1, dk, .true., s, q, qeval, cmat, &
                                llev, phase, nopen, eint, cro, ilev, jlev, klev, &
                                ered, seps)
      else
        den = (ered(3) - ered(1)) / rmlmda
        if (abs(den) .gt. 0.d0) then
          dk = (kcurr - kprev2) / den
        else
          dk = (0.d0,0.d0)
        endif
        call write_energy_block(2, kprev1, dk, .true., s, q, qeval, cmat, &
                                llev, phase, nopen, eint, cro, ilev, jlev, klev, &
                                ered, seps)
      endif

      ! Main loop: energies 3 to nener-2 (5-point stencil interior)
      if (nener .gt. 4) then
        ! At entry: kprev2=K(1), kprev1=K(2), kcurr=K(3), knext1=K(4), knext2=K(5)
        do nrg = 3, nener - 2
          write(6,'(A,1x,I0,1x,A,1x,I0)') 'MATCH_COMPLEX: energy', nrg, 'of', nener
          flush(6)
          den = (ered(nrg + 1) - ered(nrg - 1)) / rmlmda
          if (abs(den) .gt. 0.d0) then
            dk = (kprev2 - 8*kprev1 + 8*knext1 - knext2) / (6.d0 * den)
          else
            dk = (0.d0,0.d0)
          endif
          call write_energy_block(nrg, kcurr, dk, .true., s, q, qeval, cmat, &
                                  llev, phase, nopen, eint, cro, ilev, jlev, klev, &
                                  ered, seps)
          ! Shift buffer forward: K(nrg-2..nrg+2) -> K(nrg-1..nrg+3)
          kprev2 = kprev1
          kprev1 = kcurr
          kcurr = knext1
          knext1 = knext2
          if (nrg + 3 .le. nener) then
            call build_k_matrix(nrg + 3, knext2, llev, cmat, z, cvr, cvi, cro, &
                                eint, ilev, jlev, klev, nlev, ered)
          endif
        enddo
      endif

      ! Remaining energies: nener-1 and nener (3-point stencil or backward)
      if (nener .ge. 4) then
        ! After the main loop (or for nener=4 with no main loop):
        ! kprev2 = K(nener-3), kprev1 = K(nener-2), kcurr = K(nener-1), 
        ! knext1 = K(nener), knext2 = undefined or beyond nener
        do nrg = max(3, nener - 1), nener
          write(6,'(A,1x,I0,1x,A,1x,I0)') 'MATCH_COMPLEX: energy', nrg, 'of', nener
          flush(6)
          if (nrg .eq. nener) then
            ! Last energy: backward difference (K_n - K_{n-1}) / (E_n - E_{n-1})
            if (nener .ge. nrg) then
              if (nrg .eq. nener .and. nrg .ge. 4) then
                ! Use kprev1 for K(n-1) and kcurr for K(n)
                den = (ered(nrg) - ered(nrg - 1)) / rmlmda
                if (abs(den) .gt. 0.d0) then
                  dk = (kcurr - kprev1) / den
                else
                  dk = (0.d0,0.d0)
                endif
                call write_energy_block(nrg, kcurr, dk, .true., s, q, qeval, cmat, &
                                        llev, phase, nopen, eint, cro, ilev, jlev, klev, &
                                        ered, seps)
              else if (nener .eq. 3) then
                den = (ered(3) - ered(2)) / rmlmda
                if (abs(den) .gt. 0.d0) then
                  dk = (kcurr - kprev1) / den
                else
                  dk = (0.d0,0.d0)
                endif
                call write_energy_block(3, kcurr, dk, .true., s, q, qeval, cmat, &
                                        llev, phase, nopen, eint, cro, ilev, jlev, klev, &
                                        ered, seps)
              endif
            endif
          else
            ! nrg = nener-1: 3-point centered (K_{n} - K_{n-2}) / (E_n - E_{n-2})
            den = (ered(nrg + 1) - ered(nrg - 1)) / rmlmda
            if (abs(den) .gt. 0.d0) then
              dk = (knext1 - kprev1) / den
            else
              dk = (0.d0,0.d0)
            endif
            call write_energy_block(nrg, kcurr, dk, .true., s, q, qeval, cmat, &
                                    llev, phase, nopen, eint, cro, ilev, jlev, klev, &
                                    ered, seps)
            ! Shift for next iteration
            kprev2 = kprev1
            kprev1 = kcurr
            kcurr = knext1
          endif
        enddo
      else if (nener .eq. 3) then
        ! nener=3: energy 3 needs backward difference
        write(6,'(A,1x,I0,1x,A,1x,I0)') 'MATCH_COMPLEX: energy', 3, 'of', 3
        flush(6)
        den = (ered(3) - ered(2)) / rmlmda
        if (abs(den) .gt. 0.d0) then
          dk = (kcurr - kprev1) / den
        else
          dk = (0.d0,0.d0)
        endif
        call write_energy_block(3, kcurr, dk, .true., s, q, qeval, cmat, &
                                llev, phase, nopen, eint, cro, ilev, jlev, klev, &
                                ered, seps)
      endif
    endif

    close(66)
    close(67)
    close(68)
    close(69)
    close(70)
    close(72)
    write(6,'(A)') 'MATCH_COMPLEX: done'
    flush(6)
    
    deallocate(s, q, z)
    deallocate(kprev2, kprev1, kcurr, knext1, knext2, dk)
    deallocate(cmat, llev, phase, qeval)
  end subroutine match_complex

  subroutine build_k_matrix(nrg, kmat, llev, cmat, zmat, cvr, cvi, cro,  &
                            eint, ilev, jlev, klev, nlev, ered)
    use mod_common
    use mod_bessel_complex, only: bessel_complex
    use mod_frames_complex, only: frames_complex
    use mod_scratchfile_complex, only: getrec_complex
    implicit double precision (a-h,o-z)

    integer, intent(in) :: nrg
    complex*16, intent(out) :: kmat(nchan,nchan), zmat(nchan,nchan)
    double precision, intent(out) :: llev(nchan), cmat(nchan,nchan)
    double precision, intent(in) :: cvr(nvi,n), cvi(nvi,n), cro(3,n), eint(n), ered(*)
    integer, intent(in) :: ilev(n), jlev(n), klev(n), nlev(n)

    complex*16 :: y(nchan,nchan)

    call getrec_complex(y, nchan, nrg)
    call frames_complex(y, cmat, cro, jlev, klev, llev, +1)
    call bessel_complex(kmat, y, zmat, cvr, cvi, eint, ilev, jlev, klev, &
                        llev, nlev, ered(nrg))
  end subroutine build_k_matrix

  subroutine write_energy_block(nrg, kmat, dkde, has_deriv, smat, qmat,  &
                                qeval, cmat, llev, phase, nopen, eint,   &
                                cro, ilev, jlev, klev, ered, seps)
    use mod_common
    use mod_frames_complex, only: frames_complex
    use mod_smatrx_complex, only: smatrx_complex
    use mod_pes_multi, only: NPES, get_elec_offset
    implicit double precision (a-h,o-z)

    integer, intent(in) :: nrg
    complex*16, intent(in) :: kmat(nchan,nchan), dkde(nchan,nchan)
    logical, intent(in) :: has_deriv
    complex*16, intent(inout) :: smat(nchan,nchan), qmat(nchan,nchan)
    double precision, intent(inout) :: cmat(nchan,nchan), llev(nchan)
    double precision, intent(inout) :: phase(nchan), qeval(nchan)
    integer, intent(inout) :: nopen
    double precision, intent(in) :: eint(n), cro(3,n), ered(*)
    integer, intent(in) :: ilev(n), jlev(n), klev(n)
    double precision, intent(out) :: seps
    integer :: i, j, ierr, ic, ibase, ipp
    double precision :: enow
    double precision, allocatable :: eintc(:)

    enow = ered(nrg)
    
    allocate(eintc(nchan))
    do ic = 1, nchan
      ibase = (ic-1)/max(1,NPES) + 1
      ipp = ic - (ibase-1)*max(1,NPES)
      eintc(ic) = eint(ibase) + get_elec_offset(ilev(ibase), ipp)
    enddo

    call smatrx_complex(kmat, dkde, has_deriv, smat, qmat, qeval, eintc, &
                        llev, nchan, seps, ierr, enow, phase, nopen)
    deallocate(eintc)
    if (ierr .ne. 0) stop 'match_complex 1'

    ! --- Cumulative Reaction Probabilities (CRP) by product arrangement ---
    ! Same logic as cmod_match.f90:179-201 for NPES=2 complex S-matrix
    block
      integer :: ibase_c, jbase_c, ia_c, ic_c, jc_c, pp_c, qq_c
      double precision :: tp_sum_c, crp_c(3)
      do ia_c = 1, 3
        crp_c(ia_c) = 0.d0
      enddo
      do ibase_c = 1, n
        if (enow .gt. eint(ibase_c) .and. ilev(ibase_c) .eq. ireact) then
          do jbase_c = 1, n
            if (enow .gt. eint(jbase_c) .and. ilev(jbase_c) .ne. ireact) then
              tp_sum_c = 0.d0
              do pp_c = 1, max(1,NPES)
                ic_c = (ibase_c - 1)*max(1,NPES) + pp_c
                do qq_c = 1, max(1,NPES)
                  jc_c = (jbase_c - 1)*max(1,NPES) + qq_c
                  tp_sum_c = tp_sum_c + abs(smat(ic_c, jc_c))**2
                enddo
              enddo
              crp_c(ilev(jbase_c)) = crp_c(ilev(jbase_c)) + tp_sum_c
            endif
          enddo
        endif
      enddo
      if (jpar .eq. 0) then
        write(72,'(*(g0,1x))') enow/rmlmda, crp_c(2), crp_c(3)
      else
        write(72,'(*(g0,1x))') enow/rmlmda, crp_c(2)
      endif
      flush(72)
    end block

    write(67,'(*(g0,1x))') enow/rmlmda, nopen, (phase(i), i=1,nopen)
    if (has_deriv) then
      write(68,'(*(g0,1x))') enow/rmlmda, nopen, (qeval(i), i=1,nopen)
      write(69,'(*(g0,1x))') enow/rmlmda, nopen, (dble(qmat(i,i)), i=1,nopen)
      ! --- Observables output ---
      block
        double precision :: trq_fs, trq_us, e_ev, e_coll_ev, t_coll_k, e_coll_mk
        double precision :: kB_eV, emin_ev, max_tau
        integer :: ii, j_min, ic_in, jj
        complex*16 :: s11
        double precision :: p_reac, k_el, k_bohr, mu_au, m1, m2, m3, mtot_au, mred_au
        double precision :: pi, factor_k_el
        
        pi = acos(-1.d0)
        kB_eV = 8.617333262145d-5   ! eV/K
        trq_fs = 0.d0
        max_tau = -1.d30
        do ii = 1, nopen
          trq_fs = trq_fs + qeval(ii)
          if (qeval(ii) .gt. max_tau) max_tau = qeval(ii)
        enddo
        if (nopen .le. 0) max_tau = 0.d0
        trq_us = trq_fs * 1.d-9     ! fs -> μs
        e_ev = enow / rmlmda
        
        ! Find incident channel (lowest energy in ireact)
        emin_ev = 1.d30
        j_min = 1
        do ii = 1, n
          if (ilev(ii) .eq. ireact) then
            if (eint(ii)/rmlmda .lt. emin_ev) then
               emin_ev = eint(ii)/rmlmda
               j_min = ii
            endif
          endif
        enddo
        e_coll_ev = e_ev - emin_ev
        if (e_coll_ev .lt. 0.d0) e_coll_ev = 0.d0
        t_coll_k = e_coll_ev / kB_eV
        e_coll_mk = t_coll_k * 1000.d0
        write(70,'(6(ES20.12,1x),I8)') e_ev, e_coll_ev, t_coll_k, &
             dble(nopen), trq_fs, trq_us
             
        ! Calculate P_reac and k_el
        ! P_reac = 1 - |S_ii|^2: total flux leaving the incident channel
        ! (includes both reactive and inelastic scattering)
        ic_in = (j_min-1)*max(1,NPES) + 1  ! Ground electronic state of incident channel
        s11 = smat(ic_in, ic_in)
        p_reac = 1.d0 - abs(s11)**2
        if (p_reac .lt. 0.d0) p_reac = 0.d0
        
        ! k_el calculation
        m1 = 86.9091805310d0 * 1822.8885d0
        m2 = 39.963998166d0 * 1822.8885d0
        m3 = m1
        mtot_au = m1 + m2 + m3
        mred_au = sqrt(m1 * m2 * m3 / mtot_au)
        k_bohr = sqrt(2.d0 * mred_au * (e_coll_ev / 27.211386d0))
        if (k_bohr .gt. 0.d0) then
            k_el = (pi / (mred_au * k_bohr)) * abs(1.d0 - s11)**2
            k_el = k_el * 6.12627e-9 ! AU to cm^3/s
        else
            k_el = 0.d0
        endif
        
        write(71,'(4(ES20.12,1x))') e_coll_mk, max_tau, p_reac, k_el
      end block
    endif

    call frames_complex(smat, cmat, cro, jlev, klev, llev, -1)
    write(66,'(A,1x,ES20.12,1x,A,1x,ES12.4)') 'E(eV)=', enow/rmlmda, 'Seps=', seps
    write(66,'(A)') 'S-matrix complex block follows'
    do i = 1, nchan
      write(66,'(*(ES20.10,1x,ES20.10,1x))') (dble(smat(i,j)), dimag(smat(i,j)), j=1,nchan)
    enddo
  end subroutine write_energy_block

end module mod_match_complex
