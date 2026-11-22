SUBROUTINE SFunBas_2x2(rho, nbasis, jtot, parity, mega, nang, marran, chanl, nangch, philast, phipr, s, v, &
                       ThetAPH, chiaph, weight, vect, hamil, Ovrlp, hnew, xvect, loverlap, ithrho, neigmin, &
                       eigen, megamin, iunit, sthaph, cthaph, nbasiss, jrot, megamax, chitau,             &
                       v12_scale, v22_same)
!-----------------------------------------------------------------------
! Two-electronic-state (2x2) surface functions for the KRb2 APH3D chain.
!
! Mirrors SFunBas (same argument list and internal flow) but assembles the
! 2N x 2N block problem
!
!     H2 = [ H_kin + V11   ,  V12      ]      S2 = [ S  0 ]
!          [ V12^T        ,  H_kin+V22]           [ 0  S ]
!
! where the three potential matrices come from PotElmnt with the
! v1/v2/v12 arrays of CalcPot_2x2 (the real KRb2 2x2 diabatic PES in the
! APH normalization V/(usys2*rho**2)); the kinetic energy matrices are
! state independent and are computed exactly as in SFunBas.
!
! The symmetric orthogonalization is done once on the (identical)
! diagonal overlap blocks and applied to both electronic states, so the
! 2N generalized eigenproblem reduces to a single DSYEV of a
! (2*neigen) x (2*neigen) matrix.
!
! The single-state (state 1) spectrum is computed in the same call from
! the same matrices, so the no-coupling limit can be verified within one
! run: with v12_scale = 0 and v22_same = .true. the 2N spectrum must be
! exactly doubly degenerate and equal to the single-state spectrum.
! A table (rho, mega, E1(k), E2(2k-1), E2(2k), splittings, state-2
! weights w2) is appended to 'sfun2x2_dump.txt'.
!
!   v12_scale : scaling of the coupling block (0 => decoupled)
!   v22_same  : if .true., V22 := V11 (two identical electronic states,
!               used for the degeneracy unit test)
!-----------------------------------------------------------------------
USE Numeric_Kinds_Module
USE Parms_Module
USE FileUnits_Module
USE Masses_Module
USE Converge_Module
USE Gaussb_Module
USE Convrsns_Module
USE AzBzCz_Module
USE Boundary_Module
USE Numbers_Module
IMPLICIT NONE
INTEGER jtot, parity, mega, nbasis, neigen, nang, Info, k
INTEGER marran, ibasis, iang, ithrho, mbad, jeigen, jbasis
INTEGER neigmin, megamin, iunit, nbasiss, megamax, n2, neigen2
INTEGER nangch(narran+1),  jrot(nbasis), chanl(nbasis)
LOGICAL loverlap, v22_same
REAL(KIND=WP_Kind) cpu1, cpu2, ctime, overhead
SAVE xptg, wptg, xpth, wpth, temp1, temp2, oscil
SAVE dpndx, doscdz, pn, nfac, fnorm, ScaleFac, iarran
REAL(Kind=WP_Kind) rho, ScaleFac
REAL(Kind=WP_Kind) weight(mxang)
REAL(Kind=WP_Kind) Ovrlp(nbasiss,nbasiss), hamil(nbasiss,nbasiss)
REAL(Kind=WP_Kind) philast(nang, nbasis), s(nbasiss,nbasiss), v(nang)
REAL(Kind=WP_Kind) vect(nbasiss, nbasiss), phipr(nang, nbasis)
REAL(Kind=WP_Kind) ThetAPH(nang), chiaph(nang)
REAL(Kind=WP_Kind) sthaph(nang), cthaph(nang), chitau(mxang)
REAL(Kind=WP_Kind) hnew(nbasiss*nbasiss), xvect(nbasiss*nbasiss)
REAL(Kind=WP_Kind) eigen(nbasiss), v12_scale
INTEGER, ALLOCATABLE:: iarran(:)
REAL(Kind=WP_Kind), ALLOCATABLE:: xptg(:,:), wptg(:,:)
REAL(Kind=WP_Kind), ALLOCATABLE:: xpth(:,:), wpth(:,:)
REAL(Kind=WP_Kind), ALLOCATABLE:: pn(:,:,:), dpndx(:,:,:), nfac(:)
REAL(Kind=WP_Kind), ALLOCATABLE:: oscil(:,:,:), doscdz(:,:,:), fnorm(:)
REAL(Kind=WP_Kind), ALLOCATABLE:: temp2(:), temp1(:)
REAL(Kind=WP_Kind), ALLOCATABLE:: h11(:,:), h22(:,:), h12(:,:), hkin(:,:)
REAL(Kind=WP_Kind), ALLOCATABLE:: h1m(:,:), h2(:,:), str2(:,:), w1(:,:), w2m(:,:)
REAL(Kind=WP_Kind), ALLOCATABLE:: mat1(:,:), mat2(:,:), xvec1(:,:)
REAL(Kind=WP_Kind), ALLOCATABLE:: eig1(:), eig2(:), w2(:), v1p(:), v2p(:), v12p(:)
REAL(Kind=WP_Kind), ALLOCATABLE:: vect1(:,:), workq(:), xvec2(:,:)
REAL(Kind=WP_Kind), ALLOCATABLE:: vsf1(:,:), vsf2(:,:), ovl2(:,:), ovlchk(:,:), tmp2(:,:)
REAL(Kind=WP_Kind) parityf2(mxang,narran), parityfn(mxang, narran)
REAL(Kind=WP_Kind) bfaci(mxang,narran), dzdthd(mxang,narran)
REAL(Kind=WP_Kind) dlnb(mxang,narran), dxdchi(mxang,narran)
REAL(Kind=WP_Kind) dxdtha(mxang,narran), dthddtha(mxang,narran)
REAL(Kind=WP_Kind) dthddchi(mxang,narran)
REAL(Kind=WP_Kind) :: pmin1, pmin2, maxdegen, maxdiff, maxdat
EXTERNAL popt, basis, primitiv, ovlaps, potelmnt, kintheta
EXTERNAL kinchi, CalcPot_2x2

!-----------------------------------------------------------------------
! Determine printing options.
!-----------------------------------------------------------------------
LOGICAL little, medium, full
INTEGER ithcall, ithsub
DATA ithcall /0/, ithsub /0/
DATA little /.false./, medium /.false./, full /.false./
CALL popt ('sfunbas ', little, medium, full, ithcall, ithsub)

IF(.NOT.ALLOCATED(xptg))THEN
   ALLOCATE(iarran(mxang))
   ALLOCATE(xptg(mxglegn, narran), wptg(mxglegn, narran))
   ALLOCATE(xpth(mxhermt, narran), wpth(mxhermt, narran))
   ALLOCATE(pn(0:mxl, mxang, narran), dpndx(0:mxl, mxang, narran), nfac(0:2*mxl))
   ALLOCATE(oscil(0:mxn, mxang, narran),doscdz(0:mxn, mxang, narran), fnorm(0:mxn))
   ALLOCATE(temp2(mxbasiss),temp1(mxbasiss))
ENDIF
ALLOCATE(h11(nbasiss,nbasiss), h22(nbasiss,nbasiss), h12(nbasiss,nbasiss))
ALLOCATE(hkin(nbasiss,nbasiss), h1m(nbasiss,nbasiss))
ALLOCATE(vect1(nbasiss,nbasiss), xvec1(nbasiss,nbasiss))
ALLOCATE(eig1(nbasiss), eig2(2*nbasiss), w2(2*nbasiss), v1p(nang), v2p(nang), v12p(nang))
ALLOCATE(h2(2*nbasiss,2*nbasiss), str2(2*nbasiss,2*nbasiss))
ALLOCATE(w1(nbasiss,nbasiss), w2m(2*nbasiss,2*nbasiss))
ALLOCATE(mat1(nbasiss,nbasiss), mat2(2*nbasiss,2*nbasiss))
ALLOCATE(workq(80*2*nbasiss), xvec2(2*nbasiss,2*nbasiss))
ALLOCATE(vsf1(nbasiss,nbasiss), vsf2(2*nbasiss,2*nbasiss))
ALLOCATE(ovl2(nbasiss,2*nbasiss), ovlchk(nbasiss,nbasiss), tmp2(nbasiss,2*nbasiss))

!-----------------------------------------------------------------------
! Legendre and oscillator basis (identical to SFunBas).
!-----------------------------------------------------------------------
CALL cpu_time(cpu1)
CALL basis (xptg, wptg, xpth, wpth, dpndx, doscdz, ThetAPH, sthaph, cthaph, chiaph, pn, oscil,   &
            nang, nangch, parity, parityfn, parityf2, bfaci, dzdthd, dlnb, dxdchi, dxdtha,       &
            dthddtha, dthddchi, rho, mega, megamin, 1, weight, nfac, fnorm, chitau, iarran)
CALL cpu_time(cpu2)
ctime = cpu2 - cpu1 - overhead
IF(Medium)WRITE(Out_Unit,*)'   Time in sfunbas_2x2(basis) = ', ctime

!-----------------------------------------------------------------------
! Primitive basis and overlap matrix (identical to SFunBas).
!-----------------------------------------------------------------------
CALL primitiv (nang, nbasis, mega, oscil, pn, parityfn, weight, phipr, wptg, wpth, bfaci)
CALL ovlaps (nang, nbasis, chanl, nangch, phipr, Ovrlp, rho, intwt, nbasiss, jrot)

!-----------------------------------------------------------------------
! 2x2 potential on the quadrature grid, then the three potential
! matrices (PotElmnt overwrites its output array).
!-----------------------------------------------------------------------
CALL CalcPot_2x2(nang, rho, ThetAPH, chiaph, v1p, v2p, v12p)
pmin1 = MINVAL(v1p)
pmin2 = MINVAL(v2p)
CALL potelmnt (nang, nbasis, chanl, nangch, phipr, v1p, h11, rho, intwt, nbasis, nbasiss, nbasiss, jrot)
IF(v22_same)THEN
   h22 = h11
ELSE
   CALL potelmnt (nang, nbasis, chanl, nangch, phipr, v2p, h22, rho, intwt, nbasis, nbasiss, nbasiss, jrot)
ENDIF
CALL potelmnt (nang, nbasis, chanl, nangch, phipr, v12p, h12, rho, intwt, nbasis, nbasiss, nbasiss, jrot)
h12 = h12*v12_scale

!-----------------------------------------------------------------------
! Kinetic energy matrices: state independent, computed exactly as in
! SFunBas (theta part, then chi part), accumulated into hkin.
!-----------------------------------------------------------------------
CALL kintheta (nang, nbasis, mega,  oscil, pn, parityfn,  weight, philast, wptg, wpth,     &
                       bfaci, dpndx, doscdz, dzdthd, dlnb, dxdtha, dthddtha)
ScaleFac=two/(usys*rho**2)
v=ScaleFac
CALL potelmnt (nang, nbasis, chanl, nangch, philast, v, hkin, rho, intwt, nbasis, nbasiss, nbasiss, jrot)
CALL kinchi   (nang, nbasis, mega, oscil, pn, parityfn, weight, philast, wptg, wpth,     &
               bfaci, dpndx, doscdz, dzdthd, dlnb, dxdchi, dthddchi, parityf2)
ScaleFac=ScaleFac/four
DO iang = 1, nang
   v(iang) = ScaleFac/(sthaph(iang))**2
ENDDO
CALL potelmnt  (nang, nbasis, chanl, nangch, philast, v, s, rho, intwt, nbasis, nbasiss, nbasiss, jrot)
hkin = hkin + s

!-----------------------------------------------------------------------
! Assemble the 2N x 2N Hamiltonian and the (block diagonal) 2N overlap.
! NOTE: the KRb2 PES documents the diabatic matrix as
!   | Vd11   Vd12 |
!   |-Vd12   Vd22 |
! (the paper quotes a complex Hermitian 2x2 PEM).  The real symmetric
! (Hermitian) convention H12 = H21^T is used here; for the eigenvalues of
! the potential matrix the two conventions are equivalent, only the
! propagation would feel the sign.
!-----------------------------------------------------------------------
n2 = 2*nbasiss
h2 = 0.d0
DO ibasis = 1, nbasiss
   DO jbasis = 1, nbasiss
      h2(ibasis, jbasis)                 = hkin(ibasis,jbasis) + h11(ibasis,jbasis)
      h2(nbasiss+ibasis, nbasiss+jbasis) = hkin(ibasis,jbasis) + h22(ibasis,jbasis)
      h2(ibasis, nbasiss+jbasis)         = h12(ibasis,jbasis)
      h2(nbasiss+ibasis, jbasis)         = h12(jbasis,ibasis)
   ENDDO
ENDDO

!-----------------------------------------------------------------------
! Symmetric orthogonalization of the overlap.  The 2N overlap is block
! diagonal with two identical blocks, so one decomposition serves both
! the single-state and the 2N problem.
!-----------------------------------------------------------------------
vect1 = Ovrlp
CALL DSYEV('V','U',nbasiss,vect1,nbasiss,eigen,workq,SIZE(workq),info)
IF(Info/=0)THEN
   WRITE(Out_Unit,*)'Diagonalizing the Overlaps (2x2): Info=',Info
   STOP 'SFunBas_2x2 (overlap)'
ENDIF
mbad=0
DO jbasis=1,nbasiss
   IF(eigen(jbasis)>eigmin)EXIT
   mbad=mbad+1
ENDDO
neigen=nbasiss-mbad
neigen2=2*neigen
DO jbasis=mbad+1,nbasiss
   eigen(jbasis)=one/sqrt(eigen(jbasis))
   jeigen=jbasis-mbad
   DO ibasis=1, nbasiss
      s(ibasis,jeigen)=vect1(ibasis,jbasis)*eigen(jbasis)
   ENDDO
ENDDO
WRITE(Out_Unit,'(3(A,I6))')'SFunBas_2x2: NBasiss=',nbasiss,' MBad=',mbad,' NEigen=',neigen

!-----------------------------------------------------------------------
! Single-state spectrum (state 1) - reference for the 2x2 check:
! mat1 = s^T (hkin + h11) s .
!-----------------------------------------------------------------------
h1m = hkin + h11
mat1 = 0.d0
w1(1:nbasiss,1:neigen) = MATMUL(h1m, s(1:nbasiss,1:neigen))
mat1(1:neigen,1:neigen) = MATMUL(TRANSPOSE(s(1:nbasiss,1:neigen)), w1(1:nbasiss,1:neigen))
xvec1 = mat1
CALL DSYEV('V','U',neigen,xvec1,nbasiss,eig1,workq,SIZE(workq),info)
IF(Info/=0) STOP 'SFunBas_2x2 (single-state DSYEV)'

!-----------------------------------------------------------------------
! 2N generalized eigenproblem: str2 = blkdiag(s,s),
! mat2 = str2^T H2 str2 (2*neigen x 2*neigen).
!-----------------------------------------------------------------------
str2 = 0.d0
DO jbasis = 1, neigen
   DO ibasis = 1, nbasiss
      str2(ibasis, jbasis)                = s(ibasis,jbasis)
      str2(nbasiss+ibasis, neigen+jbasis) = s(ibasis,jbasis)
   ENDDO
ENDDO
w2m = MATMUL(h2, str2)
mat2 = MATMUL(TRANSPOSE(str2), w2m)
xvec2 = mat2
CALL DSYEV('V','U',neigen2,xvec2,2*nbasiss,eig2,workq,SIZE(workq),info)
IF(Info/=0) STOP 'SFunBas_2x2 (2N DSYEV)'

!-----------------------------------------------------------------------
! Electronic composition of each 2N state.  The 2N overlap is block
! diagonal, so the two electronic components of an eigenvector are
! orthogonal in the S-metric and the weight of the second electronic
! state is simply the squared norm of its coefficients in the
! S-orthonormal basis (0..1; the pair sum is 1 for identical blocks).
!-----------------------------------------------------------------------
DO k = 1, neigen2
   w2(k) = SUM(xvec2(neigen+1:neigen2, k)**2)
ENDDO

!-----------------------------------------------------------------------
! DAT overlap projection: the N adiabatic (single-surface) functions
! expressed in the 2N diabatic surface-function basis.  Both sets are
! S-orthonormal, so the (neigen x 2*neigen) matrix
!     O(k,m) = < phi_k(adiabatic) | phi_m(2x2) >
! is an isometry onto the 2N space:  O O^T = I_neigen .  This is the
! single-step limit of the paper's nDAT=160 unitary sequence and is the
! matrix used to rotate the log-derivative matrix at r_DAT.
!-----------------------------------------------------------------------
vsf1 = MATMUL(s(1:nbasiss,1:neigen), xvec1(1:neigen,1:neigen))
vsf2 = MATMUL(str2, xvec2)
tmp2 = MATMUL(Ovrlp, vsf2(1:nbasiss,1:neigen2))
ovl2 = MATMUL(TRANSPOSE(vsf1), tmp2)
ovlchk = MATMUL(ovl2, TRANSPOSE(ovl2))
DO ibasis = 1, neigen
   ovlchk(ibasis,ibasis) = ovlchk(ibasis,ibasis) - one
ENDDO
maxdat = MAXVAL(ABS(ovlchk))

!-----------------------------------------------------------------------
! Check against the single-state spectrum.
!-----------------------------------------------------------------------
maxdegen = 0.d0; maxdiff = 0.d0
DO k = 1, neigen
   maxdegen = MAX(maxdegen, ABS(eig2(2*k-1)-eig2(2*k)))
   maxdiff  = MAX(maxdiff,  ABS(eig2(2*k-1)-eig1(k)))
ENDDO
OPEN(Unit=91, File='sfun2x2_dump.txt', Status='UNKNOWN', Position='APPEND')
WRITE(91,'(A)') '# SFunBas_2x2 (2N block surface functions)'
WRITE(91,'(A,ES20.12,A,I4,A,F7.3,A,I6,A,L2,A,ES13.5)') '# rho=', rho, '  mega=', mega, &
     '  v12_scale=', v12_scale, '  neigen=', neigen, '  v22_same=', v22_same, &
     '  potmin1(HA)=', pmin1
WRITE(91,'(A)') '#    k            E1(k) [Ha]            E2(2k-1)               E2(2k)         |dE(2k-1,2k)| |dE(2k-1)-E1(k)|  w2(2k-1) w2(2k)  w2pair'
DO k = 1, neigen
   WRITE(91,'(I6,2X,3ES22.14,2X,2ES14.4,2X,2F9.5,2X,F9.5)') k, eig1(k), eig2(2*k-1), eig2(2*k), &
        ABS(eig2(2*k-1)-eig2(2*k)), ABS(eig2(2*k-1)-eig1(k)), w2(2*k-1), w2(2*k), &
        w2(2*k-1)+w2(2*k)
ENDDO
WRITE(91,'(A,ES15.6,A,ES15.6)') '# MAX |E2(2k-1)-E2(2k)| = ', maxdegen, &
     '   MAX |E2(2k-1)-E1(k)| = ', maxdiff
WRITE(91,'(A,ES15.6)') '# DAT projection isometry error max |O O^T - I| = ', maxdat
FLUSH(91)
CLOSE(91)

!-----------------------------------------------------------------------
! Safety check analogous to SFunBas: the lowest surface energy must not
! lie below the potential minimum sampled on the grid.
!-----------------------------------------------------------------------
IF(MIN(pmin1,pmin2) > MINVAL(eig2))THEN
   WRITE(Out_Unit,*)'ERROR in SFunBas_2x2: lowest eigenvalue below potential minimum'
   WRITE(Out_Unit,*)'rho=', rho, ' potmin=', MIN(pmin1,pmin2), ' eigmin2=', MINVAL(eig2)
   STOP 'SFunBas_2x2: Problem Area'
ENDIF

DEALLOCATE(h11, h22, h12, hkin, h1m, vect1, xvec1, eig1, eig2, w2, v1p, v2p, v12p)
DEALLOCATE(h2, str2, w1, w2m, mat1, mat2, workq, xvec2)
DEALLOCATE(vsf1, vsf2, ovl2, ovlchk, tmp2)
RETURN
ENDSUBROUTINE SFunBas_2x2
