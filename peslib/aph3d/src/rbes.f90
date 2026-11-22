SUBROUTINE rbes (n,z,zjn,zjnp,zyn,zynp)
  USE Numeric_Kinds_Module
  USE Numbers_Module
!
!-------------------------------------------------------------------
! this routine is from roy gordon's qcpe program.
! <rbes> -- computes ricatti bessel functions, reg + irreg.  + deriv
! l -- order of bessel function   (INTEGER(KIND=IW_Kind))
! x -- argument of bessel function    (REAL(Kind=WP_Kind))
! zjn -- returned value of regular bessel function
! zjnp -- returned value of deriv. of reg. bessel fn.
! zyn -- returned value of irregular bessel function
! zynp -- returned value of deriv. of irreg. bessel fn.
! x, zjn, zjnp, zyn, zynp -- all REAL(Kind=WP_Kind)
!-------------------------------------------------------------------
  USE FileUnits_Common_Module
  IMPLICIT NONE
  INTEGER n
  REAL(Kind=WP_Kind) sz, z, cz
  REAL(Kind=WP_Kind) zinv, aj, bj, ay, by, zjn, zjnp, zyn, zynp, xn
  REAL(Kind=WP_Kind) delf, factor, zsq, top, c, u, term, du, dterm
  REAL(Kind=WP_Kind) xj, d2np3, den, dfact, const
  sz=sin(z)
  cz=cos(z)
  zinv=one/z
  aj=sz
  bj=sz*zinv-cz
  ay=-cz
  by=-sz-cz*zinv
  !IF(n-1) 10,20,30
  IF(n<1)THEN
!-------------------------------------------------------------------
! functions of order zero
!-------------------------------------------------------------------
10 zjn=aj
  zjnp=cz
  zyn=ay
  zynp=sz
  RETURN
  ELSEIF(n==1)THEN
!-------------------------------------------------------------------
! functions of order one
!-------------------------------------------------------------------
20 zjn=bj
  zjnp=aj-zinv*bj
  zyn=by
  zynp=ay-zinv*by
  RETURN
  ELSE
!-------------------------------------------------------------------
! functions for orders greater than one
!-------------------------------------------------------------------
30 xn=n
  delf=zinv+zinv
  factor=delf+zinv
  IF(z>xn) GOTO 60
  zsq=z*z
  top=xn+xn
  xj=three
  c=zsq/three
40 zyn=by*factor-ay
  ay=by
  by=zyn
  factor=factor+delf
  xj=xj+two
  c=c*(z/xj)
  IF(xj<top) GOTO 40
  zynp=ay-zinv*by*xn
  factor=-half*zsq
  u=one
  term=one
  du=one
  dterm=one
  d2np3=top+three
  den=d2np3
  dfact=zero
50 dfact=dfact+one
  term=term*(factor/(dfact*den))
  u=u+term
  den=den+two
  dterm=dterm*(factor/(dfact*den))
  du=du+dterm
  IF(ABS(term)>1.d-9) GOTO 50
  zjn=u*c
  zjnp=(xn+one)*zjn*zinv-(z*du/d2np3)*c
  RETURN
60 const=zinv*xn
  top=const+const
70 aj=factor*bj-aj
  ay=factor*by-ay
  factor=factor+delf
  IF(factor>top) GOTO 80
  bj=factor*aj-bj
  by=factor*ay-by
  factor=factor+delf
  IF(factor<top) GOTO 70
  zjn=bj
  zyn=by
  zjnp=aj-const*bj
  zynp=ay-const*by
  RETURN
80 zjn=aj
  zyn=ay
  zjnp=bj-const*aj
  zynp=by-const*ay
  ENDIF
  RETURN
!----------------***end-rbes***---------------------------------------
  END
