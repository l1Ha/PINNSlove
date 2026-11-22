SUBROUTINE potmatn (n, rmid, w, g0, g1, g2, dmu, psi, gamma,  &
                    rhoval, xksq, energy, w0, w1, w2, wtemp,  &
                    diag, cor, asym, ndim)
!
USE FileUnits_MODULE
USE dip_MODULE
IMPLICIT NONE
!
!#include <Inc.dip>
!#include <Inc.unitno>
!
!---------------------------------------------------------------------
! written by g. a. parker
! calculates the interaction matrix
!---------------------------------------------------------------------
LOGICAL little, medium, full
INTEGER n, ithcll, ithsub, ndim, i
REAL(Kind=WP_Kind)  rmid, w(n,n), g0(n,n),g1(n,n),g2(n,n),  &
        dmu(n,n), psi(n), gamma(n), xksq(n), energy(n),     &
        w0(n,n), w1(n,n), w2(n,n),                          &
        wtemp(n,n), diag(n,n), cor(n,n), asym(n,n), rhoval(*)
SAVE

DATA ithcll/0/,ithsub/0/
DATA little/.false./,medium/.false./,full/.false./
CALL popt('potmatn ', little, medium, full, ithcll, ithsub)

CALL aphmatn(rmid, n, w, w0, w1, w2, xksq, energy, rhoval, wtemp, ndim, diag, cor, asym)

IF(lphoto.AND.rmid<=enddip)THEN
    CALL phot_mat(rmid, n, g0, g1, g2, dmu, psi, gamma, rhoval)
ENDIF
!

!         medium=.true.
IF(medium)THEN
    WRITE(Out_Unit,170) rmid
    CALL MxOutD(w, n, n, 1)
ENDIF

CALL enlvls(rmid, n, w)

! PORT (2026-10-02): gated trace of the interaction-matrix values (used to
! verify that the 2N data path consumes the storage identically to the
! single-state run): APH3D_PTRACE=1, first APH3D_PTRACE_N calls
CALL Ptrace(n, w, rmid)

! PORT (2026-10-02): gated dump of the final coupled-channel matrix
! W = (aphmatn potential) - xktot*I  (APH3D_WDUMP / _LO / _HI), used to
! validate the 2x2 inner propagator against the production chain.
CALL WdumpN(n, w, rmid)

RETURN
!----------------***end-potmatn***--------------------------------------
!
170 FORMAT(/,"the potential at r=",f10.5)

CONTAINS
SUBROUTINE Ptrace(nn, ww, rr)
! PORT (2026-10-02): trace the first interaction-matrix values (packed
! diagonal entries) at every call.
USE FileUnits_MODULE
IMPLICIT NONE
INTEGER, INTENT(IN) :: nn
REAL(Kind=WP_Kind), INTENT(IN) :: ww(nn,nn), rr
LOGICAL :: pt_on, pt_first
INTEGER :: pt_lim, pt_cnt, pt_len, pt_stat, ii, ip
REAL(Kind=WP_Kind) :: pt_sum
CHARACTER(LEN=64) :: pt_env
DATA pt_first/.true./
SAVE pt_on, pt_lim, pt_cnt
IF(pt_first)THEN
   pt_first = .false.; pt_on = .false.; pt_lim = 30; pt_cnt = 0
   CALL get_environment_variable('APH3D_PTRACE', pt_env, pt_len, pt_stat)
   IF(pt_stat==0 .AND. pt_len>0)THEN
      IF(pt_env(1:1)=='1' .OR. pt_env(1:1)=='T' .OR. pt_env(1:1)=='t') pt_on=.true.
   ENDIF
   CALL get_environment_variable('APH3D_PTRACE_N', pt_env, pt_len, pt_stat)
   IF(pt_stat==0 .AND. pt_len>0) READ(pt_env(1:pt_len),*) pt_lim
ENDIF
IF(.NOT.pt_on) RETURN
pt_cnt = pt_cnt + 1
IF(pt_cnt > pt_lim) RETURN
! packed parameters are stored row-wise in the FIRST column of w(.,.)
pt_sum = 0.d0
DO ii = 1, nn
   ip = ii*(ii+1)/2
   pt_sum = pt_sum + ww(ip,1)
ENDDO
OPEN(Unit=89, File='ptrace.txt', Status='UNKNOWN', Position='APPEND')
WRITE(89,'(I6,2ES22.14)') nn, rr, pt_sum
CLOSE(89)
END SUBROUTINE Ptrace
SUBROUTINE WdumpN(nn, ww, rr)
USE FileUnits_MODULE
IMPLICIT NONE
INTEGER, INTENT(IN) :: nn
REAL(Kind=WP_Kind), INTENT(IN) :: ww(nn,nn), rr
REAL(Kind=WP_Kind) :: xktot
LOGICAL :: wd_on, wd_first
REAL(Kind=WP_Kind) :: wd_lo, wd_hi
CHARACTER(LEN=64) :: wd_env
INTEGER :: wd_len, wd_stat, iw, jw, ijp
REAL(Kind=WP_Kind) :: wunp(nn,nn)
COMMON /utotal/ xktot
DATA wd_first/.true./
SAVE wd_on, wd_lo, wd_hi
IF(wd_first)THEN
   wd_first = .false.
   wd_on = .false.; wd_lo = -1.d10; wd_hi = 1.d10
   CALL get_environment_variable('APH3D_WDUMP', wd_env, wd_len, wd_stat)
   IF(wd_stat==0 .AND. wd_len>0)THEN
      IF(wd_env(1:1)=='1' .OR. wd_env(1:1)=='T' .OR. wd_env(1:1)=='t') wd_on=.true.
   ENDIF
   CALL get_environment_variable('APH3D_WDUMP_LO', wd_env, wd_len, wd_stat)
   IF(wd_stat==0 .AND. wd_len>0) READ(wd_env(1:wd_len),*) wd_lo
   CALL get_environment_variable('APH3D_WDUMP_HI', wd_env, wd_len, wd_stat)
   IF(wd_stat==0 .AND. wd_len>0) READ(wd_env(1:wd_len),*) wd_hi
ENDIF
IF(.NOT.wd_on) RETURN
IF(rr < wd_lo .OR. rr > wd_hi) RETURN
OPEN(Unit=96, File='wdump.txt', Status='UNKNOWN', Position='APPEND')
WRITE(96,'(A,ES22.14,A,I6)') '# W(r) at r=', rr, '  n=', nn
! the interaction matrix is stored packed (upper triangle, row-wise);
! unpack it here and subtract xktot from the diagonal (cf. pot_shift)
DO iw = 1, nn
   DO jw = 1, nn
      IF(iw <= jw) THEN
         ijp = iw + jw*(jw-1)/2
      ELSE
         ijp = jw + iw*(iw-1)/2
      ENDIF
      wunp(iw,jw) = ww(ijp,1) - MERGE(xktot, 0.d0, iw==jw)
   ENDDO
   WRITE(96,'(500ES22.14)') (wunp(iw,jw), jw = 1, nn)
ENDDO
CLOSE(96)
END SUBROUTINE WdumpN
END
