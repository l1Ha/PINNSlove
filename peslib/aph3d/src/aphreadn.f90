SUBROUTINE aphreadn(rho, n, w0, w1, w2, nrho, nmodes, rhoval, energy, rho0, &
                    rho1, rho2, lam, read_p, lammin, megamin, megamax,  mid_zero)
!
USE Numbers_Module
USE fuzzy_Module
USE FileUnits_Module
USE CSBasis_Module
USE centerho_Module
USE junk_MODULE
IMPLICIT NONE
LOGICAL little, medium, full, mid_zero, read_p, r_cstst, dbug
INTEGER n, nrho, nmodes, ithcll, ithsub, i, maxrho,lam
INTEGER megamax, lambda, lammin, megamin
REAL(Kind=WP_Kind) energy(n), rhoval(nrho), w0(n,n), w1(n,n), w2(n,n)
REAL(Kind=WP_Kind) rho, rho0, rho1, rho2
SAVE

!----------------------------------------------------------------------
! rho       hyperpsherical radius
! n         number of coupled equations
! w0        stores potential matrix at rho-0
! w1        stores potential matrix at rho-1
! w2        stores potential matrix at rho-2
! nrho      is the number of rhos per sector at which the potential
!           has been calculated.
! rhoval    an array of rho values of dimension nrho. These are the
!           rho values at which the potential was calculated.
! rho0      value of rho at which w0 was calculated.
! rho1      value of rho at which w1 was calculated.
! rho2      value of rho at which w2 was calculated.
!-----------------------------------------------------------------------
!
! start of a new sector read in values of the potential.
!-----------------------------------------------------------------------
DATA ithcll/0/, ithsub/0/
DATA little/.false./,medium/.false./,full/.false./
!      DATA fst_time /.true./

CALL popt('aphreadn', little, medium, full, ithcll, ithsub)
DBug=.false.
!-----------------------------------------------------------------------
! If rho>rhoval(nrho)THEN we are starting a new sector. Hence read
! in the rhovals for that sector and the matrices evaluated at those
! rhovalues.
!-----------------------------------------------------------------------
lambda=lam
IF(lambda==lammin.AND.lammin>megamin) lambda=megamin

IF(DBug)WRITE(Out_Unit,*)'Entering APHReadn',nmodes
IF(DBug)WRITE(Out_Unit,'(5es23.15,l2)')rho,(rhoval(i),i=1,nrho),fuzz,read_p
1 CONTINUE
IF((rho>rhoval(nrho).AND.ABS(rho-rhoval(nrho))>fuzz).or. read_p.or.newsect)THEN
   read_p=.true.
   newsect=.false.
!     2    READ(pmat_unit)rhocnt,nmodes,nrho,(rhoval(i),i=1,nrho),mid_zero
!          IF(fst_time)THEN
!           WRITE(Msg_Unit,*)' Warning: you are running without reading ',
!      >                     'the cstest variable from pmatrix file'
!           fst_time=.false.
!          ENDIF
   2 CONTINUE
   READ(pmat_unit,end=100,err=110)rhocnt,nmodes,nrho,(rhoval(i),i=1,nrho), mid_zero
   IF(DBug)Write(Out_Unit,*)rhocnt,nmodes,(rhoval(i),i=1,nrho), mid_zero
!     >                  mid_zero, r_cstst
! revised by X.L. no cstest variable in pmat_unit              
     r_cstst=cstest
   IF(DBug)WRITE(Out_Unit,*)'just after read:nmodes=',nmodes
   IF(r_cstst.neqv.cstest)THEN
      WRITE(Out_Unit,*)' cstest variable in .option is different',   &
                       'from that used to run abm program'
      STOP 'aphreadn'
   ENDIF
   IF(little)THEN
      WRITE(Out_Unit,*)' rhocnt,nmodes,nrho,mid_zero', '(rhoval(i),i=1,nrho) = ', &
                   rhocnt,nmodes,nrho,mid_zero,(rhoval(i),i=1,nrho)
      IF(DBug)WRITE(Out_Unit,*)'aphreadn'
   ENDIF
!------------------------------------------------------------------------
! Make sure that there are at least three points in the current sector
! so that we can interpolate.c
!------------------------------------------------------------------------
   IF(nrho<3)THEN
      WRITE(Out_Unit,*)'nrho is less than 3'
      STOP 'aphreadn'
   ENDIF

!------------------------------------------------------------------------
! Read in eigenenergies obtained from the surface function calculation.
!------------------------------------------------------------------------
   READ(pmat_unit,end=100,err=110)(energy(i),i=1,n)

!------------------------------------------------------------------------
! Read in the potential matrix at rho0.
!------------------------------------------------------------------------
   rho0=rhoval(1)
   IF(DBug)WRITE(Out_Unit,*)" rho0=",rho0
!         WRITE(Out_Unit,*)'aphreadn-1:nmodes=',nmodes
   CALL read_pmat(w0, n, nmodes, pmat_unit)

!------------------------------------------------------------------------
! Read in the potential matrix at rho1 unless those matrix elements are
! identically zero. If they are zero set them to zero instead of reading
! them in from the disk.
!------------------------------------------------------------------------
   rho1=rhoval(2)
   IF(DBug)WRITE(Out_Unit,*)" rho1=",rho0
   IF(.NOT.mid_zero)THEN
      IF(lambda==1.AND..NOT.cstest)THEN
         IF(DBug)WRITE(Out_Unit,*)'aphreadn-2:nmodes=',nmodes
         CALL read_pmat(w1, n, nmodes, pmat_unit)
      ELSE
         CALL vsets(n*n,w1,1,Zero)
      ENDIF
   ELSE
      CALL vsets(n*n,w1,1,Zero)
   ENDIF

!------------------------------------------------------------------------
! Read in the potential matrix at rho2.
!------------------------------------------------------------------------
   rho2=rhoval(3)
   IF(DBug)WRITE(Out_Unit,*)'aphreadn-3a:nmodes=',nmodes," rho2=",rho2
   CALL read_pmat(w2, n, nmodes, pmat_unit)
   IF(DBug)WRITE(Out_Unit,*)'aphreadn-3b:nmodes=',nmodes
   maxrho=3
   IF(DBug)THEN
      WRITE(Msg_Unit,*)'in read=',w2(1,1)
      IF(DBug)WRITE(Out_Unit,'("rho0=",ES15.7," rho1=",ES15.7," rho2=",ES15.7)')rho0,rho1,rho2
      IF(DBug)WRITE(Out_Unit,'("rhocent=",ES15.7," energy(1)=",ES15.7)')rhocnt,energy(1)
      !STOP
   ENDIF
!-----------------------------------------------------------------------
! if we are not starting in the first sector, position Pmatrx
! in the correct sector
!-----------------------------------------------------------------------
   !IF(rho<rho0)WRITE(Out_Unit,*)rho,rho0
   !IF(ABS(rho-rho0)>fuzz)WRITE(Out_Unit,*)ABS(rho-rho0),fuzz
   IF(rho<rho0.AND.ABS(rho-rho0)>fuzz)GOTO 2
!
   lambda=lambda+1
   IF(lambda-1<lam) GOTO 2
ENDIF
IF(little)THEN
   WRITE(Out_Unit,'(5es23.15,l2)')rho,(rhoval(i),i=1,nrho),fuzz,read_p
ENDIF
IF(DBug)WRITE(Out_Unit,*)'aphreadn-3c: nmodes=',nmodes
!-----------------------------------------------------------------------
! read in the potential matrix at rhoval(maxrho)
! This is used only if there are more that 3 rho values per sector.
!-----------------------------------------------------------------------
IF(rho>rhoval(maxrho).AND.ABS(rho-rhoval(nrho))>fuzz)THEN
   maxrho=maxrho+1
!-----------------------------------------------------------------------
! lambda is important to skip records if we are non starting in the
! first sector.
!-----------------------------------------------------------------------
    IF(maxrho>nrho)THEN
    DO i=lambda,megamax
       IF(DBug)WRITE(Out_Unit,*)'aphreadn-3c:nmodes=',nmodes
       CALL aphskipn(nrho, rhoval, i, pmat_unit, 0)
       IF(DBug)WRITE(Out_Unit,*)'aphreadn-3d:nmodes=',nmodes
    ENDDO
    lambda=lam
    GOTO 1
    ENDIF

!-----------------------------------------------------------------------
! replace w0 with w1
! replace w1 with w2
! read w2
!-----------------------------------------------------------------------
    rho0=rho1
    IF(DBug)WRITE(Out_Unit,*)'aphreadn-4a:nmodes=',nmodes
    CALL dcopy(n*n, w1, 1, w0, 1)
    IF(DBug)WRITE(Out_Unit,*)'aphreadn-4b:nmodes=',nmodes

    rho1=rho2

    IF(DBug)WRITE(Out_Unit,*)'aphreadn-4c:nmodes=',nmodes
    CALL dcopy(n*n, w2, 1, w1, 1)
    IF(DBug)WRITE(Out_Unit,*)'aphreadn-4d:nmodes=',nmodes

    rho2=rhoval(maxrho)
    IF(DBug)WRITE(Out_Unit,*)'aphreadn-4e:nmodes=',nmodes
    CALL read_pmat(w2, n, nmodes, pmat_unit)
ENDIF

IF(DBug)WRITE(Out_Unit,*)"Leaving APHReadn"
RETURN

100 CONTINUE
WRITE(Msg_Unit,*)' Error: EOF detected'
WRITE(Out_unit,*)' Error: EOF detected'
STOP "Stopped in APHReadn EOF detected"
110 CONTINUE
WRITE(Msg_Unit,*)' Error: ERR detected'
WRITE(Out_Unit,*)' Error: ERR detected'
STOP "Stopped in APHReadn ERR detected"
ENDSUBROUTINE APHReadn
