SUBROUTINE aphmatn(rho, n, w, w0, w1, w2, xksq, energy, rhoval, wtemp, ndim, diag, cor, asym)
!

!
!    5/ 4/1998  TRP  Modify a flawed IF statement at Stefano's urging.
!
USE FileUnits_Module
USE CSBasis_Module
USE Pbasis_MODULE
USE Numbers_Module
IMPLICIT NONE
LOGICAL little, medium, full, read_p, mid_zero
INTEGER n, nrho, ithcll, ithsub, i, nmodes, istart
INTEGER ndim, lam,megamax, jfact, lamp
INTEGER megacoup
REAL(Kind=WP_Kind) rholast, rhosave, rho, rho0, rho1, rho2
REAL(Kind=WP_Kind) energy(n), xksq(n), rhoval(*)
REAL(Kind=WP_Kind) w(n,n)
REAL(Kind=WP_Kind) w0(ndim,ndim,lammin:lammax)
REAL(Kind=WP_Kind) w1(ndim,ndim,lammin:lammax)
REAL(Kind=WP_Kind) w2(ndim,ndim,lammin:lammax)
REAL(Kind=WP_Kind) wtemp(ndim,ndim), diag(ndim,ndim,lammin:lammax)
REAL(Kind=WP_Kind) cor(ndim,ndim,lammin:lammax)
REAL(Kind=WP_Kind) asym(ndim,ndim,lammin:lammax)
SAVE
!-----------------------------------------------------------------------
! written by r. t pack. this routine transforms the wigner
! revised by X.L. check the pmat_unit writing, nrho=nrhosec
DATA nrho/3/
!-----------------------------------------------------------------------
! rho       hyperpsherical radius
! n         number of coupled equations
! w         potential matrix evaluated at rho. the potential matrix is
!           determined by a lagrange interpolation of w0, w1, w2.
! w0        stores potential matrix at rho-0
! w1        stores potential matrix at rho-1
! w2        stores potential matrix at rho-2
!-----------------------------------------------------------------------
!
! start of a new sector read in values of the potential.
!-----------------------------------------------------------------------
DATA rholast/-1.0d0/
DATA ithcll/0/,ithsub/0/
DATA little/.false./,medium/.false./,full/.false./

!      COMMON / basis / nlam(0:1000), lammin, lammax, megamin, megamax, megacoup
COMMON /totalj/ jfact
! revised by X.L.
megamax=maxmega
CALL popt('aphmatn ',little,medium,full,ithcll,ithsub)

IF(rholast==-1.d0)THEN
    rhosave = rho
ELSEIF(rho==rhosave)THEN
    rholast = -1.d0
    DO i = 1, 21
       rhoval(i) = -One
    ENDDO
ENDIF
read_p=.false.
CALL vsets(n*n,w,1,Zero)
!----------------------------------------------------------------------
! Read in potential Matrix elements.
!----------------------------------------------------------------------
istart = 0
DO lam=lammin,lammax
   CALL aphreadn(rho, nlam(lam), w0(1,1,lam), w1(1,1,lam),      &
                 w2(1,1,lam), nrho, nmodes, rhoval,             &
                 energy(istart+1),                              &
                 rho0, rho1, rho2,lam, read_p, lammin, megamin, &
                 megamax, mid_zero)
   IF(medium)THEN
      WRITE(Out_Unit,170) rho, lam
      WRITE(Out_Unit,'("w0, w1, w2=")')
      CALL MxOutD(w0, n, n, 1)
      CALL MxOutD(w1, n, n, 1)
      CALL MxOutD(w2, n, n, 1)
   ENDIF
!----------------------------------------------------------------------
! if needed, for lambda=1 add the diagonal asymetric_top term.
!----------------------------------------------------------------------
   IF(.NOT.mid_zero)THEN
      IF(lam==1.AND..NOT.cstest)THEN
          CALL diag_asym(nlam(lam), w0(1,1,lam), w1(1,1,lam),       &
                           w2(1,1,lam), rhoval, nrho, read_p)
      ENDIF
   ENDIF

!----------------------------------------------------------------------
! Interpret to find the potential at rho.
!----------------------------------------------------------------------
    CALL mat_interp(rho, nlam(lam), wtemp,                        &
                    w0(1,1,lam), w1(1,1,lam), w2(1,1,lam), xksq,  &
                    energy(istart+1), rho0, rho1, rho2)
!       WRITE(Msg_Unit,*)'after mat_interp',rho,wtemp(1,1),rho0,rho1,rho2
   IF(medium)THEN
       WRITE(Out_Unit,'("wtemp after the interpolation=")')
       CALL MxOutD(wtemp, ndim, ndim, 1)
   ENDIF

!----------------------------------------------------------------------
! read in the K=(A+B)/2 term and add K[J(J+1)-Jr(Jr+1)] to wtemp
! after scaling for the actual rho.
!----------------------------------------------------------------------
!         jfact=1
!         WRITE(Msg_Unit,*)' Warning: jfact is set up equal to 1'
   IF(jfact>0)THEN
      CALL readjtot(diag(1,1,lam), nlam(lam), lam, rhoval, read_p, lammin, megamin, megamax)
      CALL rhoscale(rho, diag(1,1,lam), wtemp, nlam(lam), nlam(lam), rhoval(2), .True.)
   ENDIF
   CALL mat_store(w, n, wtemp, nlam(lam), nlam(lam), lam, lam)
   IF(medium)THEN
      WRITE(Out_Unit,'("w=")')
      CALL MxOutD(w, n, n, 1)
   ENDIF

!----------------------------------------------------------------------
! read in the Coriolis term.
!----------------------------------------------------------------------
   lamp=lam+1
   IF(lamp<=lammax)THEN
      CALL readcor(cor(1,1,lam), nlam(lam), nlam(lamp), lam,        &
                   rhoval, read_p, lammin, megamin, megacoup)
      CALL rhoscale(rho, cor(1,1,lam), wtemp, nlam(lam), nlam(lamp), rhoval(2), .False.)
      CALL mat_store(w, n, wtemp, nlam(lamp), nlam(lam), lam, lamp)

!----------------------------------------------------------------------
! read in the Asymmetric Top term.
!----------------------------------------------------------------------
      lamp=lamp+1
      IF(lamp<=lammax)THEN
          CALL readasy(asym(1,1,lam), nlam(lam), nlam(lamp), lam,    &
                       rhoval, read_p, lammin, megamin, megacoup)
          CALL rhoscale(rho, asym(1,1,lam), wtemp, nlam(lam), nlam(lamp), rhoval(2), .False.)
          CALL mat_store(w, n, wtemp, nlam(lamp), nlam(lam), lam, lamp)
      ENDIF
   ENDIF
    istart = istart + nlam(lam)
ENDDO

IF(istart/=n)THEN
   WRITE(Msg_Unit,*)'Error istart/=n'
   WRITE(Msg_Unit,*)'istart,n=',istart,n
   WRITE(*,*)'DBG aphmatn: istart,n,lammin,lammax=',istart,n,lammin,lammax
   WRITE(*,*)'DBG aphmatn: nlam(lammin:lammax)=',(nlam(i),i=lammin,lammax)
   STOP 'aphmatn'
ENDIF

IF(read_p)THEN
   DO lam = lammax+1, megamax
      CALL aphskipn(nrho, rhoval, lam, pmat_unit, 0)
      IF(jfact>0)THEN
         CALL aphskipn(nrho, rhoval, lam, jtot_unit, 1)
      ENDIF
   ENDDO
   IF((lammax-lammin)>=1)THEN
      DO lam = lammax+1, megacoup
         CALL aphskipn(nrho, rhoval, lam, cor_unit, 1)
      ENDDO
      IF((lammax-lammin)>=2)THEN
         DO lam = lammax+1, megacoup
            CALL aphskipn(nrho, rhoval, lam, asym_unit, 1)
         ENDDO
      ENDIF
   ENDIF
ENDIF
rholast=rho

RETURN
170 FORMAT(1x,"the potential at r=",f10.5,' lamda=',i3)
ENDSUBROUTINE aphmatn
