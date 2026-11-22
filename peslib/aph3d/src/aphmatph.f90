SUBROUTINE aphmatph(rho, n, g0, g1, g2)
!
! interpolates dipole moment function,
! reads lower state wave function,
! determines dipole operator clebsch gordon algebra,
! combines these three elements to create inhomogeneity gamma
!
!              -m. braunstein, Jan. 1991
!
USE Numeric_Kinds_Module
USE FileUnits_MODULE
USE Narran_Module
USE TotalEng_Module
USE Masses_Module
USE Qall_Module
USE fuzzy_MODULE
USE centerho_MODULE
USE Numbers_Module
USE constnts_MODULE
USE dip_MODULE

IMPLICIT NONE
INTEGER n, nrho, nmodes, i, j, maxrho, npointemp, ia, ib, ic
LOGICAL little, medium
REAL(Kind=WP_Kind) rhoval(21), rholast, wig, rho, rho0, rho1, rho2, sum
REAL(Kind=WP_Kind) a0, a1, a2, rhop, cleb, ra, rb, rc, fact, efield
REAL(Kind=WP_Kind) g0(n,n),g1(n,n),g2(n,n),dmu(nstate,nstate),psi(nstate), gamma(nstate)

!--------------------------------------------------------------------
! rho       hyperpsherical radius
! n         number of coupled equations
! dmu       dipole moment matrix evaluated at rho.
!           the dipole moment matrix is
!           determined by a lagrange interpolation of w0, w1, w2.
! g0        stores dipole moment matrix at rho-0
! g1        stores dipole moment  matrix at rho-1
! g2        stores dipole moment  matrix at rho-2
! psi       wave function for lower state
! gamma     inhomogeneous term
!---------------------------------------------------------------------
! start of a new sector READ in values of the dipole matrix element
!---------------------------------------------------------------------
DATA little/.false./,medium/.false./
DATA rhoval/21*-1.d0/,nrho/1/
DATA rholast/-1.d0/
!
DATA wig/0.d0/
!
! READ control
!
WRITE(Out_Unit,*)'rho,rholast',rho,rholast
1 IF((rho>rhoval(nrho).AND.ABS(rho-rhoval(nrho))>fuzz) .or.ABS(rho-rholast)<fuzz)THEN
!
! READ from dipole moment file, tmutp
!
2  READ(tmutp_unit)rhocnt,nmodes,nrho,(rhoval(i),i=1,nrho)
!
IF(little)THEN
    WRITE(Out_Unit,*)' rhocnt,nmodes,nrho,(rhoval(i),i=1,nrho) = ', rhocnt,nmodes,nrho,(rhoval(i),i=1,nrho)
ENDIF
WRITE(Out_Unit,*)'rhoval = ',(rhoval(i),i=1,nrho)
!
rho0=rhoval(1)
DO 10 i=1,nsfuncu
READ(tmutp_unit)(g0(i,j),j=1,nsfuncl)
10 CONTINUE
!
rho1=rhoval(2)
DO 20 i=1,nsfuncu
READ(tmutp_unit)(g1(i,j),j=1,nsfuncl)
20 CONTINUE
!
rho2=rhoval(3)
DO 30 i=1,nsfuncu
READ(tmutp_unit)(g2(i,j),j=1,nsfuncl)
30 CONTINUE
maxrho=3
IF(medium)WRITE(Out_Unit,*)rho0,rho1,rho2,rhocnt
! if we are not starting in the first sector, position Pmatrx
! in the correct sector
IF(rholast<0.0d0.AND.(rho<rho0.AND.ABS(rho-rho0)>fuzz))GOTO 2
ENDIF
WRITE(Out_Unit,*)g0(1,1),g1(1,1),g2(1,1),rho0,rho1,rho2
!---------------------------------------------------------------------
! READ in the potential matrix at rhoval(maxrho)
!---------------------------------------------------------------------
!   40 IF(rho>rhoval(maxrho).AND.ABS(rho-rhoval(nrho))>fuzz)THEN
!
!      WRITE(Out_Unit,*) 'aphmatph - error'
!
!      maxrho=maxrho+1
!      IF(maxrho>nrho)GOTO 1
!      rho0=rho1
!      rho1=rho2
!      rho2=rhoval(maxrho)
!      DO 50 j=1,nmodes
!      IF(j<=n)THEN
!      DO 60 i=1,n
!      g0(i,j)=g1(i,j)
!  60  g1(i,j)=g2(i,j)
!      READ(tmutp_unit)(g2(i,j),i=1,n)
!      ELSE
!      READ(tmutp_unit)
!      ENDIF
!   50 CONTINUE
!      ENDIF
!--------------------------------------------------------------------
! set lagrange constants to determine electronic dipole moment at rho.
!--------------------------------------------------------------------
IF(ABS(rho1-rho2)<fuzz)THEN
    a0=(rho-rho1)/(rho0-rho1)
    a1=(rho-rho0)/(rho1-rho0)
    a2=0.0d0
ELSE
IF(ABS(rho0-rho1)<fuzz)THEN
    a0=0.0d0
    a1=(rho-rho2)/(rho1-rho2)
    a2=(rho-rho1)/(rho2-rho1)
ELSE
    a0=(rho-rho1)*(rho-rho2)/(rho0-rho1)/(rho0-rho2)
    a1=(rho-rho0)*(rho-rho2)/(rho1-rho0)/(rho1-rho2)
    a2=(rho-rho0)*(rho-rho1)/(rho2-rho0)/(rho2-rho1)
ENDIF
ENDIF
!---------------------------------------------------------------------
! interpolate electronic dipole moment at rho
!---------------------------------------------------------------------
DO i=1,nsfuncu
   DO j=1,nsfuncl
      dmu(i,j)=a0*g0(i,j)+a1*g1(i,j)+a2*g2(i,j)
   ENDDO
ENDDO
!----------------------------------------------------------------------
! READ in wave function of lower state
!----------------------------------------------------------------------
READ(inho_unit) npointemp,rhop,(psi(i),i=1,nsfuncl)
!---------------------------------------------------------------------
! find dipole moment clebsch gordon constant
!---------------------------------------------------------------------
IF(wig==0.d0)THEN
!
WRITE(Out_Unit,*)'start interpolation of dipole matrix element'
!
! DO this once
!
IF(prf/=pri+prd)THEN
    STOP 'check parities - aphmatph'
ENDIF
!
ia=0
ib=0
ic=0
IF(lamd==0)ia=1
IF(lamf==0)ib=1
IF(lami==0)ic=1
!
wig=cleb(ji,jd,jf,mi,md,mf)*sqrt(REAL((2*jd+1)*(2*ji+1),WP_Kind)&
   /(16.d0*pi*pi*REAL((2*jf+1)*(1+ia)*(1+ib)*(1+ic),WP_Kind)))
!
ra=1.d0
rb=1.d0
rc=1.d0
IF(mod(jf+prf,2)/=0)ra=-ra
IF(mod(ji+lami+pri,2)/=0)rb=-rb
IF(mod(jd+lamd+prd,2)/=0)rc=-rc
fact=0.d0
IF(lamf==lami+lamd)THEN
fact=cleb(ji,jd,jf,lami,lamd,lamf)
ENDIF
IF(lami==0.AND.lamf==0.AND.lamd==0)THEN
fact=fact+ra*cleb(ji,jd,jf,0,0,0)
ENDIF
IF(lamd==lami+lamf)THEN
fact=fact+rb*cleb(ji,jd,jf,-lami,lamd,lamf)
ENDIF
IF(lamf==lami-lamd)THEN
fact=fact+rc*cleb(ji,jd,jf,lami,-lamd,lamf)
ENDIF
wig=wig*fact
!
! evaluate value of the field, multiply by a constant
!
efield=sqrt((8.d0*pi*Etot)/clight)
efield=-efield*sqrt((8.d0*pi*pi)/3.d0)
!
! account for body frame dipole axis
!
IF(iaxis==1)efield=-1.d0*efield
IF(iaxis==2)THEN
efield=-1.d0*efield
WRITE(Out_Unit,*)'y-axis, remember factor of i in S matrix'
ENDIF
!
! multiply constants together
!
wig=wig*efield*usys2
!
! end of wig evaluation
!
ENDIF
!
! multiply lower state wave function by these constants
!
DO 90 i=1,nsfuncl
psi(i)=psi(i)*wig
90    CONTINUE
!
! multiply interpolated dipole moment by lower state wave function
! to obtain inhomogeneity at rho
!
DO 100 i=1,nsfuncu
sum=0.d0
DO 110 j=1,nsfuncl
    sum=sum+dmu(i,j)*psi(j)
110     CONTINUE
gamma(i)=sum
100   CONTINUE
!
! WRITE matrix to file
!
WRITE(Gamma_Unit)(gamma(i),i=1,nsfuncu)
!
! WRITE output to unit 6 for debug
!
WRITE(Out_Unit,*)'rho = ',rho
WRITE(Out_Unit,*)'psi(1) = ',psi(1)
WRITE(Out_Unit,*)'wig = ', wig
WRITE(Out_Unit,*)'dmu(1,1) =',dmu(1,1)
WRITE(Out_Unit,*)'gamma(1) = ',gamma(1)
!
rholast=rho
!
RETURN
80  FORMAT(//,'maxrho,nrho,rho,rho0,rho1,rho2: ',2i5,4f10.5)
81  FORMAT(//,'rho0=',f10.5)
82  FORMAT(//,'rho1=',f10.5)
83  FORMAT(//,'rho2=',f10.5)
84  FORMAT(//,'rho=',f10.5)
85  FORMAT(//,'nrho is less than 3')
ENDSUBROUTINE aphmatph
