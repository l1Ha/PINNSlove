SUBROUTINE printm (printit, itype, name, array, n, m, nmax )
USE Numeric_Kinds_Module
!
!     <printm> -- driver for mprint, cprint, iprint
!
IMPLICIT NONE
INTEGER itype, i, j, nl, nc, it, nmax, mp, nh, m, n, l
CHARACTER(LEN=4) lx(10), lv(4), name
REAL(Kind=WP_Kind) array(*)
EXTERNAL printit
DATA lx / "1", "2", "3", "4", "5", "6", "7", "8", "9", "0" /
it=itype*nmax/8
lv(1)=name
mp=10
IF(itype==4)mp=22
nh=1
nl=m
i=0
j=0
   1  nc=min0(nl,mp)
CALL printit(lv(1), array(nh), n, nc, nmax)
nl=nl-nc
IF(nl<=0)RETURN
i=i+1
lv(4)=lx(mod(i,10)+1)
nh=nh+it*nc
l=i/10
IF(l==j)GOTO 1
lv(3)=lx(l)
j=j+1
GOTO 1
END
