SUBROUTINE APHRead(rho, n, w, w0, w1, w2, jderiv, nrho, nmodes,   &
                    rhoval, energy, rho0, rho1, rho2, rholast)
!
USE Numbers_Module
USE FileUnits_MODULE
USE Narran_Module
USE nstate_MODULE
USE TotalEng_Module
USE Masses_Module
USE Qall_Module
USE fuzzy_MODULE
USE centerho_MODULE
IMPLICIT NONE
INTEGER n, ms, ithcll, ithsub, jderiv, nrho, i, j, maxrho, nmodes
LOGICAL little, medium, full, mid_zero
REAL(Kind=WP_Kind) rho, rholast, rho0, rho1, rho2
REAL(Kind=WP_Kind) energy(nstate), rhoval(21)
REAL(Kind=WP_Kind) w(n,n), w0(n,n), w1(n,n), w2(n,n)
!----------------------------------------------------------------------
! rho       hyperpsherical radius
! n         number of coupled equations
! w         potential matrix evaluated at rho. the potential matrix is
!           determined by a lagrange interpolation of w0, w1, w2.
! jderiv     calculates the j-derivative of w at rho.
! w0        stores potential matrix at rho-0
! w1        stores potential matrix at rho-1
! w2        stores potential matrix at rho-2
!-----------------------------------------------------------------------
!
! start of a new sector read in values of the potential.
!-----------------------------------------------------------------------
DATA ithcll/0/,ithsub/0/
DATA little/.false./,medium/.false./,full/.false./

CALL popt('aphread ',little,medium,full,ithcll,ithsub)
IF(jderiv<=0)THEN
1    CONTINUE
    IF((rho>rhoval(nrho).AND.ABS(rho-rhoval(nrho))>fuzz)     &
                                 .or.ABS(rho-rholast)<fuzz)THEN
2          READ(pmat_unit)rhocnt,nmodes,nrho,(rhoval(i),i=1,nrho), mid_zero
    IF(little)THEN
        WRITE(Out_Unit,*)' rhocnt,nmodes,nrho,', '(rhoval(i),i=1,nrho) = ',   &
                      rhocnt,nmodes,nrho,(rhoval(i),i=1,nrho)
    ENDIF
    IF(nrho<3)GOTO 88
    READ(pmat_unit)(energy(i),i=1,n)
    rho0=rhoval(1)
    DO 10 j=1,nmodes
        IF(j<=n)THEN
            READ(pmat_unit)(w0(i,j),i=j,n)
            DO 9 i = j,n
                w0(j,i)=w0(i,j)
9             CONTINUE
        ENDIF
        IF(j>n)READ(pmat_unit)
10       CONTINUE

    rho1=rhoval(2)
    IF(.NOT.mid_zero)THEN
        DO j=1,nmodes
            IF(j<=n)THEN
                READ(pmat_unit)(w1(i,j),i=j,n)
                DO i =j,n
                w1(j,i)=w1(i,j)
                ENDDO
            ENDIF
            IF(j>n)READ(pmat_unit)
        ENDDO
    ELSE
        CALL vsets(n*n,w1,1,Zero)
    ENDIF

    rho2=rhoval(nrho)
    DO 30 j=1,nmodes
        IF(j<=n)THEN
            READ(pmat_unit)(w2(i,j),i=j,n)
            DO 29 i =j,n
                w2(j,i)=w2(i,j)
29             CONTINUE
        ENDIF
        IF(j>n)READ(pmat_unit)
30       CONTINUE
    maxrho=3
    IF(medium)WRITE(Out_Unit,*)rho0,rho1,rho2,rhocnt,energy(1)
! if we are not starting in the first sector, position Pmatrx
! in the correct sector
    IF(rholast<0.0d0.AND. (rho<rho0.AND.ABS(rho-rho0)>fuzz))GOTO 2
    ENDIF
!-----------------------------------------------------------------------
! read in the potential matrix at rhoval(maxrho)
!-----------------------------------------------------------------------
40    IF(rho>rhoval(maxrho).AND.ABS(rho-rhoval(nrho))>fuzz)THEN
    maxrho=maxrho+1
    WRITE(Out_Unit,*)'maxrho=',maxrho
    IF(maxrho>nrho)GOTO 1
    rho0=rho1
    rho1=rho2
    rho2=rhoval(maxrho)
    DO 50 j=1,nmodes
        IF(j<=n)THEN
            DO 60 i=1,n
                w0(i,j)=w1(i,j)
                w1(i,j)=w2(i,j)
60              CONTINUE
            READ(pmat_unit)(w2(i,j),i=1,n)
        ELSE
            READ(pmat_unit)
        ENDIF
50       CONTINUE
    ENDIF
ENDIF

RETURN
88 WRITE(Out_Unit,85)
STOP 'aphread'
85 FORMAT(//,'nrho is less than 3')

ENDSUBROUTINE APHRead
