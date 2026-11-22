module mod_bessel_complex
  use mod_complex_linalg, only: solve_linear_complex
  use Math_Functions, only: rbesjy, rbessk
  implicit none
  private
  public :: bessel_complex

contains

  subroutine bessel_complex(kmat, ymat, zmat, cvr, cvi, eint, ilev, jlev, klev, llev, nlev, ered)
    use mod_common
    use mod_qvib, only: qvib
    use mod_pvib, only: pvib
    use mod_pes_multi, only: NPES, PESMODE, get_elec_offset
    implicit double precision (a-h,o-z)

    complex*16, intent(out) :: kmat(nchan,nchan), zmat(nchan,nchan)
    complex*16, intent(in) :: ymat(nchan,nchan)
    double precision, intent(in) :: cvr(nvi,n), cvi(nvi,n), eint(n), llev(nchan), ered
    integer, intent(in) :: ilev(n), jlev(n), klev(n), nlev(n)

    double precision, allocatable :: wvi(:), xvi(:), pvi(:,:), fr(:), fvi(:), gvi(:)
    double precision, allocatable :: fa(:), fb(:), fc(:), fd(:), a(:,:), b(:,:), erow(:), ecol(:), srow(:), scol(:)
    complex*16, allocatable :: x(:,:), z(:,:), zr(:,:), xr(:,:)
    integer :: i, j, k, nj, kvi, jbase, ibase, kbase, nactive, ii, jj, iidx, ipp, jpp, kpp
    integer :: ij, ji, ki, kj, kk, idx
    double precision :: pi, piby2, rtrho, tworho, smin, tmin, tmax, weight, ta, cta, sta, ra, sa
    double precision :: ell, psq, p, arg, cj, dj, ej, cy, dy, ey, ck, dk, ek, rtp, atr, btr, ctr, dtr
    logical :: same_family
    integer :: ierr
    integer, allocatable :: aidx(:)
    logical, allocatable :: isact(:)
    double precision :: maxz, ztol, ysum, maxy, ytol, tiny

    allocate(wvi(mvi), xvi(mvi), pvi(nvi,2), fr(n), fvi(n), gvi(n), fa(nchan), fb(nchan), fc(nchan), fd(nchan))
    allocate(a(nchan,0:nvi-1), b(nchan,0:nvi-1), erow(nchan), ecol(nchan), srow(nchan), scol(nchan))
    allocate(x(nchan,nchan), z(nchan,nchan))

    a = 0.d0
    b = 0.d0
    x = (0.d0, 0.d0)
    z = (0.d0, 0.d0)

    pi = acos(-1.d0)
    piby2 = 0.5d0*pi
    rtrho = sqrt(rmax)
    tworho = 2.d0*rmax
    if (smin_asym_override .gt. 0.d0) then
      smin = smin_asym_override
    else if (trim(PESMODE) .eq. 'KRb2' .or. trim(PESMODE) .eq. 'KRB2') then
      smin = 1.d0
    else
      smin = 0.d0
    endif
    tmin = 0.d0
    tmax = asin(min(1.d0, smax/rmax))
    call qvib(tmin, tmax, mvi, wvi, xvi)

    do kvi = 1, mvi
      weight = rtrho*wvi(kvi)
      ta = xvi(kvi)
      cta = cos(ta)
      sta = sin(ta)
      ra = rmax*cta
      sa = rmax*sta
      if (sa .lt. smax) then
        call pvib(tmin, ta, tmax, cvr, nvi, n, fr, 0)
        do i = 1, n
          fr(i) = weight*fr(i)
        enddo
        call pvib(smin, sa, smax, cvi, nvi, n, fvi, 0)
        call pvib(smin, sa, smax, cvi, nvi, n, gvi, 1)
        do i = 1, nchan
          ibase = (i-1)/max(1,NPES) + 1
          ipp = i - (ibase-1)*max(1,NPES)
          ell = llev(ibase)
          psq = ered - (eint(ibase) + get_elec_offset(ilev(ibase), ipp))
          if (psq .gt. 0.d0) then
            p = sqrt(psq)
            arg = p*ra
            call rbesjy(ell, arg, cj, dj, ej, cy, dy, ey)
            if (kvi .eq. 1) then
              ecol(i) = ej
              scol(i) = exp(ej)
              erow(i) = ey
              srow(i) = exp(-ey)
            else
              cj = exp(ej-ecol(i))*cj
              dj = exp(ej-ecol(i))*dj
              cy = exp(ey-erow(i))*cy
              dy = exp(ey-erow(i))*dy
            endif
            rtp = sqrt(p)
            atr = cj/rtp
            btr = cy/rtp
            ctr = dj*rtp
            dtr = dy*rtp
          else
            ! Cap psq so that p * rmax <= 500.0 to prevent overflow in exp(ek - erow(i))
            if (psq .lt. -(500.d0/rmax)**2) then
              psq = -(500.d0/rmax)**2
            endif
            p = sqrt(-psq)
            arg = p*ra
            call rbessk(ell, arg, ck, dk, ek)
            if (kvi .eq. 1) then
              ecol(i) = 0.d0
              scol(i) = 0.d0
              erow(i) = ek
              srow(i) = 0.d0
            else
              ck = exp(ek-erow(i))*ck
              dk = exp(ek-erow(i))*dk
            endif
            rtp = sqrt(p)
            atr = 0.d0
            btr = ck/rtp
            ctr = 0.d0
            dtr = dk*rtp
          endif
          fa(i) = atr*fvi(ibase)
          fb(i) = btr*fvi(ibase)
          fc(i) = cta*ctr*fvi(ibase)+sta*atr*gvi(ibase)
          fd(i) = cta*dtr*fvi(ibase)+sta*btr*gvi(ibase)
          fc(i) = fa(i)/tworho + fc(i)
          fd(i) = fb(i)/tworho + fd(i)
        enddo

        do j = 1, nchan
          jbase = (j-1)/max(1,NPES) + 1
          jpp = j - (jbase-1)*max(1,NPES)
          ij = ilev(jbase)
          jj = jlev(jbase)
          kj = klev(jbase)
          nj = nlev(jbase)
          do i = 1, nchan
            ibase = (i-1)/max(1,NPES) + 1
            ipp = i - (ibase-1)*max(1,NPES)
            ii = ilev(ibase)
            ji = jlev(ibase)
            ki = klev(ibase)
            
            ! Enforce that electronic states don't mix artificially in the matching
            same_family = (ji.eq.jj .and. ki.eq.kj .and. ipp.eq.jpp)
            if (same_family) then
              same_family = (ii.eq.ij)
              if (.not.same_family .and. jpar.eq.0) then
                if (ii+ij .eq. 5) same_family = .true.
              endif
            endif
            if (same_family) then
              a(i,nj) = a(i,nj) + fr(ibase)*fa(j)
              b(i,nj) = b(i,nj) + fr(ibase)*fb(j)
              x(i,j) = x(i,j) - dcmplx(fr(ibase)*fc(j), 0.d0)
              z(i,j) = z(i,j) - dcmplx(fr(ibase)*fd(j), 0.d0)
            endif
          enddo
        enddo
      endif
    enddo

    do j = 1, nchan
      jbase = (j-1)/max(1,NPES) + 1
      jpp = j - (jbase-1)*max(1,NPES)
      ij = ilev(jbase)
      jj = jlev(jbase)
      kj = klev(jbase)
      nj = nlev(jbase)
      do k = 1, nchan
        kbase = (k-1)/max(1,NPES) + 1
        kpp = k - (kbase-1)*max(1,NPES)
        ii = ilev(kbase)
        ji = jlev(kbase)
        ki = klev(kbase)
        
        same_family = (ji.eq.jj .and. ki.eq.kj .and. kpp.eq.jpp)
        if (same_family) then
          same_family = (ii.eq.ij)
          if (.not.same_family .and. jpar.eq.0) then
            if (ii+ij .eq. 5) same_family = .true.
          endif
        endif
        if (same_family) then
          do i = 1, nchan
            x(i,j) = x(i,j) + ymat(i,k) * dcmplx(a(k,nj), 0.d0)
            z(i,j) = z(i,j) + ymat(i,k) * dcmplx(b(k,nj), 0.d0)
          enddo
        endif
      enddo
    enddo

    allocate(aidx(nchan), isact(nchan))
    isact = .false.
    maxz = 0.d0
    maxy = 0.d0
    do i = 1, nchan
      maxz = max(maxz, abs(z(i,i)))
      do j = 1, nchan
        maxy = max(maxy, abs(ymat(i,j)))
      enddo
    enddo
    ztol = max(1.d-300, maxz*1.d-12)
    ytol = max(1.d-300, maxy*1.d-12)
    nactive = 0
    do i = 1, nchan
      ysum = 0.d0
      do j = 1, nchan
        ysum = ysum + abs(ymat(i,j))
      enddo
      if (abs(z(i,i)) .gt. ztol .or. ysum .gt. ytol) then
        nactive = nactive + 1
        aidx(nactive) = i
        isact(i) = .true.
      endif
    enddo
    if (nactive .eq. nchan) then
      call solve_linear_complex(z, x, nchan, nchan, ierr)
    else
      allocate(zr(nactive,nactive), xr(nactive,nchan))
      do ii = 1, nactive
        iidx = aidx(ii)
        do j = 1, nchan
          xr(ii,j) = x(iidx,j)
        enddo
      enddo
      do ii = 1, nactive
        do jj = 1, nactive
          zr(ii,jj) = z(aidx(ii), aidx(jj))
        enddo
      enddo
      call solve_linear_complex(zr, xr, nactive, nchan, ierr)
      if (ierr .eq. 0) then
        do ii = 1, nactive
          iidx = aidx(ii)
          do j = 1, nchan
            x(iidx,j) = xr(ii,j)
          enddo
        enddo
        do i = 1, nchan
          if (.not. isact(i)) then
            do j = 1, nchan
              x(i,j) = (0.d0,0.d0)
            enddo
          endif
        enddo
      endif
      deallocate(zr, xr)
    endif
    if (ierr .ne. 0) then
      tiny = max(1.d-300, maxz*1.d-12)
      do i = 1, nchan
        if (abs(z(i,i)) .lt. tiny) z(i,i) = dcmplx(tiny, 0.d0)
      enddo
      call solve_linear_complex(z, x, nchan, nchan, ierr)
      if (ierr .ne. 0) stop 'bessel_complex 1'
    endif
    deallocate(aidx, isact)

    do j = 1, nchan
      do i = 1, nchan
        x(i,j) = dcmplx(srow(i),0.d0) * x(i,j) * dcmplx(scol(j),0.d0)
      enddo
    enddo

    kmat = x
    zmat = z
  end subroutine bessel_complex

end module mod_bessel_complex
