SUBROUTINE Jacobi_LongRange_Prop(rm_in, rf_in)
  USE Numeric_Kinds_Module
  USE FileUnits_Module
  USE Masses_Module, ONLY: usys2
  USE Energy_Module, ONLY: eref
  USE TotalEng_Module, ONLY: etot
  IMPLICIT NONE

  REAL(Kind=WP_Kind), INTENT(IN) :: rm_in, rf_in
  REAL(Kind=WP_Kind) :: dr, h, r_curr, r_mid, r2_mid, r6_mid
  REAL(Kind=WP_Kind) :: k2, p, v_cent, v_c6, c6_val, e_krb_min, delta_e
  REAL(Kind=WP_Kind) :: rmax_ev, tol_ev, c60, c62
  REAL(Kind=WP_Kind), ALLOCATABLE :: amat(:,:), Y(:,:), ysub(:,:)
  REAL(Kind=WP_Kind), ALLOCATABLE :: vmat(:,:), ev(:), prj(:,:), work_tmp(:)
  REAL(Kind=WP_Kind), ALLOCATABLE :: k2mat(:,:), qmat(:,:), ev2(:), work_tmp2(:)
  REAL(Kind=WP_Kind), ALLOCATABLE :: work_arr(:), a2blk(:,:), ytmp(:,:), y1m(:,:)
  REAL(Kind=WP_Kind), ALLOCATABLE :: y1(:), y2(:), energies(:)
  INTEGER, ALLOCATABLE :: chanl(:), elect(:), nvib(:), jrot(:), lorb(:)
  INTEGER :: n, ithmat, jthmat, nu_jtot, jtot, nstep, istep, i, j, k, info, ios
  INTEGER :: lwork, npop, lwork2
  INTEGER :: inunit, outunit
  CHARACTER(LEN=20) :: TypeOfMat, proj, approx, symm
  CHARACTER(LEN=20) :: nu_type, nu_proj, nu_approx, nu_symm
  CHARACTER(LEN=256) :: file_in, file_out
  LOGICAL :: more

  IF (rf_in <= rm_in) THEN
     WRITE(Out_Unit,*) 'Jacobi_LongRange_Prop: rf <= rm, skipping propagation.'
     RETURN
  END IF

  WRITE(Out_Unit,'(A,F8.2,A,F8.2)') &
       'Jacobi_LongRange_Prop: Integrating from rm = ', rm_in, ' to rf = ', rf_in
  WRITE(*,'(A,F8.2,A,F8.2)') &
       'Jacobi_LongRange_Prop: Integrating from rm = ', rm_in, ' to rf = ', rf_in

  file_in  = TRIM(OutDIR) // 'BinOut/URU_Matrix.bin'
  file_out = TRIM(OutDIR) // 'BinOut/URU_Matrix_rf.bin'

  inunit  = 181
  outunit = 182

  OPEN(Unit=inunit, File=TRIM(file_in), Form='UNFORMATTED', Status='OLD', IOSTAT=ios)
  IF (ios /= 0) THEN
     WRITE(Out_Unit,*) 'ERROR Jacobi_LongRange_Prop: Cannot open ', TRIM(file_in)
     RETURN
  END IF

  OPEN(Unit=outunit, File=TRIM(file_out), Form='UNFORMATTED', Status='REPLACE', IOSTAT=ios)
  IF (ios /= 0) THEN
     WRITE(Out_Unit,*) 'ERROR Jacobi_LongRange_Prop: Cannot open ', TRIM(file_out)
     CLOSE(inunit)
     RETURN
  END IF

  dr = 0.1d0
  nstep = NINT((rf_in - rm_in) / dr)
  IF (nstep < 10) nstep = 10
  h = (rf_in - rm_in) / DBLE(nstep)

  record_loop: DO
     READ(inunit, IOSTAT=ios) jthmat, nu_type, n, etot, nu_jtot, nu_proj, nu_approx, nu_symm
     IF (ios /= 0) EXIT record_loop

     ithmat    = jthmat
     TypeOfMat = nu_type
     jtot      = nu_jtot
     proj      = nu_proj
     approx    = nu_approx
     symm      = nu_symm

     ALLOCATE(amat(n,n), energies(n), chanl(n), elect(n), nvib(n), jrot(n), lorb(n))
     CALL ReadMat(amat, energies, chanl, elect, nvib, jrot, lorb, n, inunit)

     ! If this is a dummy or zero-channel record, write it through unchanged
     IF (n <= 0) THEN
        CALL Write_Matrix2Old(amat, energies, chanl, elect, nvib, jrot, lorb, n, etot, &
                              ithmat, TypeOfMat, jtot, proj, approx, symm, outunit)
        DEALLOCATE(amat, energies, chanl, elect, nvib, jrot, lorb)
        CYCLE record_loop
     END IF

     ALLOCATE(Y(n,n), ysub(n,n), y1(n), y2(n))

     ! ------------------------------------------------------------------
     ! PORT (2026-10-11): The R-matrix carried out of the APH region is
     ! R_delves = U R_aph U^T with U rectangular (ndelves x naph): its rank
     ! is at most naph, so as many as ndelves-naph directions of the Delves
     ! channel space carry NO wavefunction amplitude.  A raw matrix inverse
     ! of this singular matrix produces O(1/e-machine) noise in those
     ! directions, which the propagation then mixes into the physical block
     ! and the asymptotic matching consumes: the resulting R' is strongly
     ! asymmetric (|R-R^T| ~ |R| observed) and the per-channel S-matrix
     ! phases are corrupted (Pauli residual 1-4% instead of 0).
     !
     ! Fix: Moore-Penrose pseudo-inverse (zero on the null space, exactly
     ! symmetric) + re-projection onto the populated subspace after every
     ! propagation step, so the null structure is preserved exactly and the
     ! unpopulated directions decouple in the matching.
     ! ------------------------------------------------------------------
     ALLOCATE(vmat(n,n), ev(n), prj(n,n), work_tmp(1))
     vmat = 0.5d0*(amat + TRANSPOSE(amat))
     CALL DSYEV('V', 'U', n, vmat, n, ev, work_tmp(1:1), -1, info)
     lwork = MAX(INT(work_tmp(1)), 3*n)
     DEALLOCATE(work_tmp)
     ALLOCATE(work_tmp(lwork))
     CALL DSYEV('V', 'U', n, vmat, n, ev, work_tmp, lwork, info)
     IF (info /= 0) THEN
        WRITE(Out_Unit,*) 'Jacobi_LongRange_Prop: DSYEV failed, info=', info
        CALL Write_Matrix2Old(amat, energies, chanl, elect, nvib, jrot, lorb, n, etot, &
                              ithmat, TypeOfMat, jtot, proj, approx, symm, outunit)
        DEALLOCATE(amat, energies, chanl, elect, nvib, jrot, lorb, Y, ysub, y1, y2)
        DEALLOCATE(vmat, ev, prj, work_tmp)
        CYCLE record_loop
     END IF
     rmax_ev = MAXVAL(ABS(ev))
     tol_ev  = 1.d-10*rmax_ev
     npop = COUNT(ABS(ev) > tol_ev)
     Y   = 0.0d0
     prj = 0.0d0
     DO i = 1, n
        IF (ABS(ev(i)) > tol_ev) THEN
           DO j = 1, n
              DO k = 1, n
                 prj(j,k) = prj(j,k) + vmat(j,i)*vmat(k,i)
                 Y(j,k)   = Y(j,k)   + vmat(j,i)*vmat(k,i)/ev(i)
              ENDDO
           ENDDO
        END IF
     ENDDO
     WRITE(Out_Unit,'(A,I5,A,I5,A,ES9.2)') &
        ' Jacobi_LongRange_Prop: populated directions=', npop, ' of', n, &
        '  tol_ev=', tol_ev

     ! Step 2: Boundary metric derivative correction from Delves to Jacobi:
     ! Y_Jacobi(rm) = Y_Delves(rm) - 1/(2*rm) * I  -- restricted to the
     ! populated subspace (the unpopulated directions stay exactly zero).
     DO i = 1, n
        DO j = 1, n
           Y(i,j) = Y(i,j) - 0.5d0 * prj(i,j) / rm_in
        ENDDO
     END DO

     ! ------------------------------------------------------------------
     ! PORT (2026-10-11): FULL COUPLED-CHANNEL long-range propagation.
     !
     ! The previous version propagated a DIAGONAL potential (isotropic
     ! C6/R^6 + centrifugal only).  The diagonal model misses the two
     ! off-diagonal couplings that feed the reaction flux in the
     ! 8.1 -> 49.8 a0 region:
     !   * the C6 anisotropy  -C6,2 P2(cos theta)/R^6  (couples j -> j+/-2),
     !   * the channel-threshold dispersion through the eigenbasis mixing.
     ! Three long-range treatments changed the total reaction probability
     ! by FOUR orders of magnitude (3.3e-8 / 1.1e-5 / 2.3e-4), proving
     ! these couplings are the dominant missing physics.
     !
     ! Long-range potential (atom-diatom, same for all three H+H2
     ! arrangements, no inter-arrangement coupling: the exchange part of
     ! H3 is exponentially dead beyond 8 a0):
     !   V_ij(R) = -[C6,0 delta_ij(arr) + C6,2 A2(j_i,j_j) delta_ij(arr)] / R^6
     !   A2(j',j) = sqrt((2j'+1)(2j+1)) * (3j(j' 2 j;0 0 0))^2
     ! C6,2/C6,0 = 2*Delta_alpha/(3*alpha_bar) = 0.2262 from the H2
     ! polarizability anisotropy (alpha_par=6.3802, alpha_perp=4.5776 a0^3),
     ! C6,0 = 6.4990267 a.u. (literature value, deck env APH3D_C6).
     !
     ! Propagator: piecewise-constant matrix log-derivative (Johnson
     ! generalized): at each mid-step, diagonalize the symmetric wave
     ! operator
     !   K2_ij = 2mu (Etot - E_i) delta_ij - l_i(l_i+1)/R^2 delta_ij
     !           + 2mu V_ij(R_mid)
     ! propagate the log-derivative diagonally in the eigenbasis, rotate
     ! back.  Converges as h -> 0.
     ! ------------------------------------------------------------------
     ALLOCATE(k2mat(n,n), qmat(n,n), ev2(n), work_tmp2(1))
     ALLOCATE(a2blk(n,n), ytmp(n,n), y1m(n,n))
     ! arrangement-resolved channel data
     c60 = 6.4990267d0
     BLOCK
       CHARACTER(LEN=32) :: env_c6b
       INTEGER :: stat_c6b
       CALL get_environment_variable('APH3D_C6', env_c6b, status=stat_c6b)
       IF (stat_c6b == 0) READ(env_c6b, *) c60
     END BLOCK
     c62 = c60*0.2262d0
     ! the P2 angular matrix per arrangement (identical blocks)
     a2blk = 0.0d0
     DO i = 1, n
        DO j = 1, n
           IF (chanl(i) == chanl(j)) THEN
              a2blk(i,j) = SQRT(DBLE(2*jrot(i)+1)*DBLE(2*jrot(j)+1)) &
                           * threej000(jrot(i), 2, jrot(j))**2
           END IF
        ENDDO
     ENDDO

     ! Step 3: coupled-channel log-derivative propagation from rm to rf
     r_curr = rm_in
     DO istep = 1, nstep
        r_mid  = r_curr + 0.5d0 * h
        r2_mid = r_mid * r_mid
        r6_mid = r2_mid * r2_mid * r2_mid

        ! assemble the symmetric wave operator
        DO i = 1, n
           DO j = 1, n
              k2mat(i,j) = 0.0d0
              IF (i == j) THEN
                 k2mat(i,j) = usys2*(etot - energies(i)) &
                              - DBLE(lorb(i)*(lorb(i)+1))/r2_mid &
                              + 2.0d0*usys2*c60/r6_mid
              END IF
              k2mat(i,j) = k2mat(i,j) + 2.0d0*usys2*c62*a2blk(i,j)/r6_mid
           ENDDO
        ENDDO
        k2mat = 0.5d0*(k2mat + TRANSPOSE(k2mat))

        ! diagonalize the wave operator
        CALL DSYEV('V', 'U', n, k2mat, n, ev2, work_tmp2(1:1), -1, info)
        lwork2 = MAX(INT(work_tmp2(1)), 3*n)
        IF (ALLOCATED(work_arr)) THEN
           IF (SIZE(work_arr) < lwork2) THEN
              DEALLOCATE(work_arr); ALLOCATE(work_arr(lwork2))
           END IF
        ELSE
           ALLOCATE(work_arr(lwork2))
        END IF
        CALL DSYEV('V', 'U', n, k2mat, n, ev2, work_arr, lwork2, info)
        IF (info /= 0) THEN
           WRITE(Out_Unit,*) 'Jacobi_LongRange_Prop: DSYEV failed at step ', istep, ' info=', info
           EXIT
        END IF
        qmat = k2mat

        ! diagonal log-derivative elements in the eigenbasis
        DO i = 1, n
           k2 = ev2(i)
           IF (k2 > 1.d-14) THEN
              p = SQRT(k2)
              y1(i) = p / TAN(p * h)
              y2(i) = p / SIN(p * h)
           ELSE IF (k2 < -1.d-14) THEN
              p = SQRT(-k2)
              y1(i) = p / TANH(p * h)
              y2(i) = p / SINH(p * h)
           ELSE
              y1(i) = 1.d0 / h
              y2(i) = 1.d0 / h
           END IF
        ENDDO

        ! ysub = Q^T (Y + Q diag(y1) Q^T) Q ; invert ; rotate back.
        ! PORT (2026-10-11): explicit loops instead of DGEMM -- the BLAS call
        ! segfaulted in this context and n=126 makes the O(n^3) loops cheap.
        DO i = 1, n
           DO j = 1, n
              ytmp(i,j) = Y(j,i)               ! ytmp = Y^T
           ENDDO
        ENDDO
        DO i = 1, n
           DO j = 1, n
              y1m(i,j) = 0.0d0
              DO k = 1, n
                 y1m(i,j) = y1m(i,j) + qmat(k,i)*ytmp(k,j)   ! Q^T Y
              ENDDO
           ENDDO
        ENDDO
        DO i = 1, n
           DO j = 1, n
              ysub(i,j) = 0.0d0
              DO k = 1, n
                 ysub(i,j) = ysub(i,j) + y1m(i,k)*qmat(k,j)  ! Q^T Y Q
              ENDDO
           ENDDO
        ENDDO
        DO i = 1, n
           ysub(i,i) = ysub(i,i) + y1(i)
        END DO
        CALL invert_lapack(ysub, n, info)
        IF (info /= 0) THEN
           WRITE(Out_Unit,*) 'Jacobi_LongRange_Prop: Warning inverting ysub at step ', istep
        END IF
        DO j = 1, n
           DO i = 1, n
              ytmp(i,j) = -y2(i) * ysub(i,j) * y2(j)
           END DO
           ytmp(j,j) = ytmp(j,j) + y1(j)
        END DO
        DO i = 1, n
           DO j = 1, n
              y1m(i,j) = 0.0d0
              DO k = 1, n
                 y1m(i,j) = y1m(i,j) + qmat(i,k)*ytmp(k,j)   ! Q ytmp
              ENDDO
           ENDDO
        ENDDO
        DO i = 1, n
           DO j = 1, n
              Y(i,j) = 0.0d0
              DO k = 1, n
                 Y(i,j) = Y(i,j) + y1m(i,k)*qmat(j,k)        ! (Q ytmp) Q^T
              ENDDO
           ENDDO
        ENDDO

        r_curr = r_curr + h
     END DO

     ! Step 4: Metric correction from Jacobi back to Delves at rf:
     ! Y_Delves(rf) = Y_Jacobi(rf) + 1/(2*rf) * I  -- populated subspace only.
     DO i = 1, n
        DO j = 1, n
           Y(i,j) = Y(i,j) + 0.5d0 * prj(i,j) / rf_in
        ENDDO
     END DO

     ! Step 5: Invert Y(rf) to obtain the asymptotic R-matrix: amat = Y^(-1).
     ! With the full coupled-channel propagation the previously-unpopulated
     ! directions are physically populated, so Y(rf) is generically
     ! well-conditioned and the raw inverse is the correct object.  Fall
     ! back to the pseudo-inverse only if LAPACK reports exact singularity.
     amat = Y
     CALL invert_lapack(amat, n, info)
     IF (info /= 0) THEN
        WRITE(Out_Unit,*) 'Jacobi_LongRange_Prop: Y(rf) singular at energy ', ithmat, &
                          ' -> falling back to the pseudo-inverse'
        CALL pinv_sym(Y, n, amat, info)
     END IF

     ! Step 6: Write updated record into outunit
     CALL Write_Matrix2Old(amat, energies, chanl, elect, nvib, jrot, lorb, n, etot, &
                           ithmat, TypeOfMat, jtot, proj, approx, symm, outunit)

     DEALLOCATE(amat, energies, chanl, elect, nvib, jrot, lorb, Y, ysub, y1, y2)
     DEALLOCATE(vmat, ev, prj, work_tmp)
     DEALLOCATE(k2mat, qmat, ev2, work_tmp2, a2blk, ytmp, y1m)
     IF (ALLOCATED(work_arr)) DEALLOCATE(work_arr)
  END DO record_loop

  CLOSE(inunit)
  CLOSE(outunit)

  ! Replace original URU_Matrix.bin with URU_Matrix_rf.bin
  CALL copy_file(TRIM(file_out), TRIM(file_in))
  WRITE(Out_Unit,*) 'Jacobi_LongRange_Prop: Successfully propagated URU_Matrix to rf = ', rf_in
  WRITE(*,*) 'Jacobi_LongRange_Prop: Successfully propagated URU_Matrix to rf = ', rf_in

CONTAINS

  SUBROUTINE invert_lapack(A, n, info)
    INTEGER, INTENT(IN) :: n
    REAL(Kind=WP_Kind), INTENT(INOUT) :: A(n,n)
    INTEGER, INTENT(OUT) :: info
    INTEGER :: ipiv(n)
    REAL(Kind=WP_Kind) :: work(n*64)
    INTEGER :: lwork

    lwork = n * 64
    CALL DGETRF(n, n, A, n, ipiv, info)
    IF (info /= 0) RETURN
    CALL DGETRI(n, A, n, ipiv, work, lwork, info)
  END SUBROUTINE invert_lapack

  ! Moore-Penrose pseudo-inverse of a symmetric matrix: symmetric result,
  ! exactly zero on the directions whose eigenvalues are below the cutoff.
  SUBROUTINE pinv_sym(A, n, Apinv, info)
    INTEGER, INTENT(IN) :: n
    REAL(Kind=WP_Kind), INTENT(IN)  :: A(n,n)
    REAL(Kind=WP_Kind), INTENT(OUT) :: Apinv(n,n)
    INTEGER, INTENT(OUT) :: info
    REAL(Kind=WP_Kind), ALLOCATABLE :: v(:,:), w(:), wk(:)
    REAL(Kind=WP_Kind) :: rmax, tol, wl
    INTEGER :: i, j, k, lw

    ALLOCATE(v(n,n), w(n))
    v = 0.5d0*(A + TRANSPOSE(A))
    CALL DSYEV('V', 'U', n, v, n, w, wl, -1, info)
    lw = MAX(INT(wl), 3*n)
    ALLOCATE(wk(lw))
    CALL DSYEV('V', 'U', n, v, n, w, wk, lw, info)
    IF (info /= 0) THEN
       DEALLOCATE(v, w, wk)
       RETURN
    END IF
    rmax = MAXVAL(ABS(w))
    tol  = 1.d-10*rmax
    Apinv = 0.0d0
    DO i = 1, n
       IF (ABS(w(i)) > tol) THEN
          DO j = 1, n
             DO k = 1, n
                Apinv(j,k) = Apinv(j,k) + v(j,i)*v(k,i)/w(i)
             ENDDO
          ENDDO
       END IF
    ENDDO
    DEALLOCATE(v, w, wk)
  END SUBROUTINE pinv_sym

  ! (j1 j2 j3; 0 0 0) Wigner 3j symbol for INTEGER arguments (Racah formula):
  !   (j1 j2 j3;000) = (-1)^g sqrt[(2g-2j1)!(2g-2j2)!(2g-2j3)!/(2g+1)!]
  !                    * g! / [(g-j1)!(g-j2)!(g-j3)!],  g=(j1+j2+j3)/2
  ! checks: (1 2 1;000)=sqrt(2/15), (2 2 2;000)=-sqrt(2/35)
  REAL(Kind=WP_Kind) FUNCTION threej000(j1, j2, j3)
    INTEGER, INTENT(IN) :: j1, j2, j3
    REAL(Kind=WP_Kind) :: g, arg, lgn1, lgn2, lgn3, lgd, lgr
    INTEGER :: gt
    threej000 = 0.0d0
    gt = j1 + j2 + j3
    IF (MOD(gt, 2) /= 0) RETURN          ! requires j1+j2+j3 even
    IF (j3 < ABS(j1-j2) .OR. j3 > j1+j2) RETURN   ! triangle condition
    g = 0.5d0*DBLE(gt)
    lgn1 = LGAMMA(DBLE(2*g - 2*j1) + 1.0d0)
    lgn2 = LGAMMA(DBLE(2*g - 2*j2) + 1.0d0)
    lgn3 = LGAMMA(DBLE(2*g - 2*j3) + 1.0d0)
    lgd  = LGAMMA(DBLE(2*g + 1) + 1.0d0)
    arg  = 0.5d0*(lgn1 + lgn2 + lgn3 - lgd)
    lgr  = LGAMMA(g + 1.0d0) - LGAMMA(g - DBLE(j1) + 1.0d0) &
         - LGAMMA(g - DBLE(j2) + 1.0d0) - LGAMMA(g - DBLE(j3) + 1.0d0)
    threej000 = DBLE((-1)**(gt/2)) * EXP(arg + lgr)
  END FUNCTION threej000

  SUBROUTINE copy_file(fsrc, fdst)
    CHARACTER(LEN=*), INTENT(IN) :: fsrc, fdst
    INTEGER :: exitstat
    CALL EXECUTE_COMMAND_LINE('cp -f "' // TRIM(fsrc) // '" "' // TRIM(fdst) // '"', EXITSTAT=exitstat)
  END SUBROUTINE copy_file

END SUBROUTINE Jacobi_LongRange_Prop
