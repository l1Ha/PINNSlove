SUBROUTINE readjtot(diag, n, lam, rhoval, read_p, lammin, megamin, megamax)
  USE FileUnits_MODULE
  USE Fuzzy_MODULE
  USE Masses_Module
  IMPLICIT NONE
!
  CHARACTER(LEN=11) ident
  INTEGER mega, megap, nmodes, nmodesp, lam, j, i, n, megamax, lambda, jfact, lammin, megamin
  REAL(Kind=WP_Kind) diag(n,n), rhoval(*), rhocent, rhofst, prefac, a, b
  LOGICAL little, read_p
  COMMON /totalj/ jfact

!-----------------------------------------------------------------------
! prefac = 2*mu* [J*(J+1)-Jref*(Jref+1)]
diag = (a+b)/2
!-----------------------------------------------------------------------
  prefac=usys2*jfact
  lambda=lam
1    IF(lambda==lammin.AND.lammin>megamin) lambda=megamin
  IF(read_p)THEN
3       READ(jtot_unit)rhocent,rhofst,mega,megap,nmodes, nmodesp,ident

     IF(mega/=lambda) goto 10
     IF(little)THEN
        WRITE(Out_Unit,*)' rhocent,rhoval,mega,megap',  &
                         ', nmodes,nmodesp,ident=',        &
                    rhocent,rhofst,mega,megap,nmodes,      &
                    nmodesp,ident
     ENDIF

     DO 2 j=1,nmodes
        IF(j<=n)THEN
           READ(jtot_unit,err=100)(diag(i,j),i=j,n)
           DO 5 i=j,n
              diag(i,j) = diag(i,j)*prefac
              diag(j,i) = diag(i,j)
5           CONTINUE
        ELSE
           READ(jtot_unit,err=100)
        ENDIF
2     CONTINUE
     lambda=lambda+1
     IF(lambda-1<lam) GOTO 3

  ENDIF

!-----------------------------------------------------------------------
  IF(rhofst<rhoval(1).AND.ABS(rhofst-rhoval(1))>fuzz)THEN
     DO i=lambda,megamax
        CALL aphskipn(1, rhoval, i, jtot_unit, 1)
     ENDDO
     lambda=lam
     GOTO 1
  ENDIF

  RETURN
10    WRITE(Msg_Unit,*)' Warning: mega/=lambda. mega, lambda=',  mega,lambda
  STOP 'readiagj'
100   WRITE(Msg_Unit,*)' Error reading diag_jtot_ABM file'
  STOP 'readiagj'

  END
