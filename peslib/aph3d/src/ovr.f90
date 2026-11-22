SUBROUTINE ovr (id2, th2, ch2, f2, fbar, fbarp, id1, th1, ch1, f1, mega, energy, sovrlap, thetaval, chivals, naph, rho1, rho2,  &
                irho, ntheta, nchi, numnp1, numnp2, nel1, nel2, iunit2, iunit1, lammin, ovrlapmx, f3, IntChi,                     &
                IntTheta, igrid, xpt, wht, value, nquad,  h, p, hcoef)
USE Numbers_Module
USE FileUnits_Module
USE nstate_Module
USE qlam_Module
USE rhoovr_Module
!
!   $RCSfile: ovr.f,v $ $Revision: 1.15 $
!   $Date: 89/11/16 11:03:17 $
!   $State: Stable $
!

!                  P U R P O S E    O F    S U B R O U T I N E

!    This routine calculates overlap matrix elements between the
!    surface functions at the previous rho value (rho1) with the
!    surface functions at the current rho value (rho2).

!                  I N P U T    A R G U M E N T S

!    naph         Number of surface functions.
!    rho1         Previous rho at which surface functions were
!                 calculated.
!    rho2         Current rho at which surface functions were
!                 calculated.
!    irho         ith-rho.
!    f1           Surface functions calculated at rho1.
!    f2           Surface functions calculated at rho2.
!    id1          List of node numbers at rho1.
!    th1          List of theta indices at rho1.
!    ch1          List of chi indices at rho1
!    id2          List of node numbers at rho2.
!    th2          List of theta indices at rho2.
!    ch2          List for chi indices at rho2.
!    thetaval     List of theta values.
!    chivals      List of chi values.
!    ntheta       Number of theta values.
!    nchi         Number of chi values.
!    fbar         Temporary storage array.
!    fbarp        Temporary storage array.
!    mega         Current value of lambda.
!    energy       Temporary storage array.
!    numnp1
!    numnp2
!    nel1
!    nel2
!    iuinit2
!    iuinit1
!    lammin
!    ovrlampmx
!    f3
!    chic
!    IntTheta
!    igrid
!    xpt
!    wht
!    value
!    nquad
!    h
!    p
!    hcoef

!                 O U T P U T    A R G U M E N T S

!    sovrlap       Overlap matrix.

IMPLICIT NONE
LOGICAL little, medium, full, irdsurf
INTEGER   id1, th1, ch1, id2, th2, ch2, naph,i, irho, numnp1, numnp2, istart, mega,  ii, ntheta, nchi, nel1, nel2, naph2
INTEGER   mxelem, mxnod, ithcall, ithsub, node, k, ielem, igrid, iunit1, iunit2, lammin, j, nquad
INTEGER IntChi(numnp2),   IntTheta(numnp2)
REAL(Kind=WP_Kind) rho1, rho2, fbar, fbarp, thetaval, chivals, f1, f2,  energy, sovrlap, ovrlapmx, f3
REAL(Kind=WP_Kind) value, xpt, wht, h(9,nquad,nquad), p(2,9,nquad,nquad), hcoef(9, nquad, nquad)
DIMENSION fbar(naph, nquad, nquad), fbarp(naph, nquad, nquad),  id1(9, nel1), th1(numnp1), xpt(nquad), wht(nquad)
DIMENSION id2(9, nel2), th2(numnp2), ch2(numnp2), ch1(numnp1),  f2(numnp2,1), energy(naph), value(nquad, nquad)
DIMENSION sovrlap(naph,naph), igrid(*), thetaval(ntheta), chivals(nchi), f1(numnp1, naph)
DIMENSION ovrlapmx(naph), f3(numnp2)

!      E X T E R N A L S
EXTERNAL interpol, popt, rdwrsurf, matxelem, mxoutl, basisck, raph, scaleit
DATA irdsurf/.false./

!-----------------------------------------------------------------------
! Determine printing options.
!-----------------------------------------------------------------------
DATA ithcall /0/, ithsub /0/
DATA little /.false./, medium /.false./, full /.false./
CALL popt ('ovr     ', little, medium, full, ithcall, ithsub)

!-----------------------------------------------------------------------
! READ in surface functions for rho 1.
!-----------------------------------------------------------------------
      mxelem=1715        !Tempfix
          mxnod=1715         !Tempfix
      CALL rdwrsurf (id1, th1, thetaval, ch1, chivals, f1, energy,  nel1, &
      numnp1, naph, ntheta, nchi, iunit1, 1, mxnod, mxelem, mega, lammin)

!-----------------------------------------------------------------------
! READ in surface functions for rho 2.
!-----------------------------------------------------------------------
IF(irdsurf)THEN
      CALL rdwrsurf (id2, th2, thetaval, ch2, chivals, f2, energy, nel2, &
      numnp2, naph, ntheta, nchi, iunit2, 1, mxnod, mxelem, mega, lammin)
ENDIF

!-----------------------------------------------------------------------
! check to see IF the same grid has been used for both rho values.
!-----------------------------------------------------------------------
IF(numnp1/=numnp2.or.nel1/=nel2) GOTO 20
DO 15 ielem = 1, nel2
   DO 10 k = 1, 9
      IF(id1(k, ielem)/=id2(k, ielem)) GOTO 20
   10    CONTINUE
   15 CONTINUE
DO 16 node = 1,numnp2
   IF(th1(node)/=th2(node).or.ch1(node)/=ch2(node))GOTO 20
   16 CONTINUE

!-----------------------------------------------------------------------
! This is reached only IF the mesh at both rho values are identical.
! Hence we DO not need to interpolate.
!-----------------------------------------------------------------------
WRITE(Out_unit, 220)
GOTO 30

!-----------------------------------------------------------------------
! The meshes are not identical.
! Interpolate functions to put them on the same grid.
!-----------------------------------------------------------------------
   20 CALL interpol (f1, th1, ch1, id1, nel1, numnp1, th2, ch2, id2, nel2, &
                     numnp2, naph, ntheta, nchi, igrid, f3, f1, IntChi, IntTheta, energy)

!-----------------------------------------------------------------------
! Calculate Overlap matrix elements.
!-----------------------------------------------------------------------
 30   CONTINUE
IF(full)THEN
   WRITE(Out_unit,*)'Compare the first surface function at rho1 and ', &
   'rho2',rho1,rho2
   WRITE(Out_unit,*)'f1=',(f1(k,1),k=1,numnp1)
   WRITE(Out_unit,*)'f2=',(f2(k,1),k=1,numnp2)
ENDIF
rho1t = rho1
rho2t = rho2
DO 91 istart=1,naph,100
   CALL raph(f2,naph,numnp2,istart)
   naph2=min(naph+1-istart, 100)
CALL matxelem (nel2, naph, f1, th2, ch2, naph2, f2, numnp2, id2, fbar, fbarp, thetaval, chivals, 1, ntheta, nchi,     &
                     sovrlap(1,istart), rho1, rho2, xpt, wht, value,nquad, h, p, hcoef)

   91 CONTINUE

!-----------------------------------------------------------------------
! Scale the overlap integrals.
!-----------------------------------------------------------------------
CALL scaleit(naph*naph, sovrlap, 1, Two)
!-----------------------------------------------------------------------
! Perform a completeness check for the current basis.
!-----------------------------------------------------------------------
CALL basisck (sovrlap, ovrlapmx, naph)
!-----------------------------------------------------------------------
! Store the Overlap matrix on disk.
!-----------------------------------------------------------------------
WRITE(Out_unit,*) 'irho,rho1,rho2 = ',irho,rho1,rho2

WRITE(dvr19) irho-1,rho1,rho2
WRITE(dvr19) naph, rho1, rho2
DO j = 1, naph
   WRITE(dvr19) (sovrlap(i,j), i = 1, naph)
ENDDO

!-----------------------------------------------------------------------
! IF desired print the Overlap matrix.
!-----------------------------------------------------------------------
IF(medium)THEN
   DO ii=1,naph
      tee(ii)=ii
      lambda(ii)=mega
   ENDDO
   WRITE(Out_unit, 230) rho1, rho2
   CALL MxOutL(sovrlap, naph, naph, 0, 'aph  ', 'aph  ')
ENDIF

RETURN
!
  130 FORMAT(1x, 11f10.5)
  160 FORMAT(/, '***error***:jtot1/=jtot2.or.jtot1/=jtot:', 3i10)
  170 FORMAT(/, '***error***:mega1/=mega2.or.mega1/=mega:', 3i10)
  180 FORMAT(/, '***error***:nsym1/=nsym2.or.nsym1/=symmetry:',3i10)
  220 FORMAT(/, 'The grid for the two rho values are identical')
  230 FORMAT(/, 'overlap matrix elements: rho1, rho2, ',2es14.7)
END
