SUBROUTINE readasy(asym, n1, n2, lam, rhoval, read_p, lammin,  megamin, megamax)

!
USE FileUnits_MODULE
USE fuzzy_MODULE
USE Masses_Module
IMPLICIT NONE
!
CHARACTER(LEN=11) ident
INTEGER mega, megap, nmodes, nmodesp, lam, j, i, n1, n2, megamax
INTEGER lambda, lammin, megamin
REAL(Kind=WP_Kind) asym(n1,n2), rhoval(*), rhocent, rhofst, prefac
REAL(Kind=WP_Kind) xtrafac, f_cor, f_asy
LOGICAL little, read_p
COMMON /xtra/ xtrafac, f_cor(0:1000), f_asy(0:1000)

!-----------------------------------------------------------------------
! f_asy:   asymetric_top factor
!          = 1/(2mu)*[jtot(jtot+1) - mega(mega+1)]**(1/2)* mega =  0
!            *[jtot(jtot+1) - (mega+1)(mega+2)]**(1/2)/4*
!            *[1 + (-1)**(jtot+parity)]/2**(1/2)
!          = 1/(2mu)*[jtot(jtot+1) - mega(mega+1)]**(1/2)*
!            *[jtot(jtot+1) - (mega+1)(mega+2)]**(1/2)/4   mega/=0
!-----------------------------------------------------------------------
lambda=lam
prefac=usys2*f_asy(lam)
1    IF(lambda==lammin.AND.lammin>megamin) lambda=megamin
IF(read_p)THEN
3       READ(asym_unit)rhocent,rhofst,mega,megap,nmodes,nmodesp,ident

    IF(megap/=lambda+2) goto 120
    IF(n1>nmodes) goto 130
    IF(n2>nmodesp) goto 140
    IF(little)THEN
    WRITE(Out_Unit,*)' rhocent,rhoval(1),mega,megap', ', nmodes,nmodesp,ident=',  &
                         rhocent,rhofst,mega,megap,nmodes, nmodesp,ident
    ENDIF

    DO 2 j=1,nmodesp
    IF(j<=n2)THEN
        READ(asym_unit,err=100)(asym(i,j),i=1,n1)
        DO 4 i=1,n1
            asym(i,j) = asym(i,j)*prefac
4             CONTINUE
    ELSE
        READ(asym_unit,err=100)
    ENDIF
2     CONTINUE

    lambda=lambda+1
    IF(lambda-1<lam) GOTO 3

ENDIF
!-----------------------------------------------------------------------

IF(rhofst<rhoval(1).AND.ABS(rhofst-rhoval(1))>fuzz)THEN
    DO i=lambda,megamax-2
       CALL aphskipn(1, rhoval, i, asym_unit, 1)
    ENDDO
    lambda=lam
    GOTO 1
ENDIF

RETURN
!-----------------------------------------------------------------------
120   WRITE(Msg_Unit,*)' Warning: megap/=(lambda+2). megap,',   &
                  ' (lambda+2)=',megap,lambda+2
STOP 'readasy'
130   WRITE(Msg_Unit,*)' Warning: n1>nmodes. n1, nmodes=', n1,nmodes
STOP 'readasy'
140   WRITE(Msg_Unit,*)' Warning: n2>nmodesp. n2, nmodesp=', n2,nmodesp
STOP 'readasy'
100   WRITE(Msg_Unit,*)' Error reading asym_top_ABM file'
WRITE(Msg_Unit,*)' n1, nmodes=',n1,nmodes
STOP 'readasy'

END
