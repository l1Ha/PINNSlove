SUBROUTINE factor
!
USE FileUnits_Module
USE factrl_Module
USE vector_Module
USE vwaste_Module
INTEGER I
REAL(Kind=WP_Kind) a
!
DO 10 i=2, 64
   vscr(i)=log(fpn(i))
 10   CONTINUE

DO 20 i=65, 101
   vscr(i)=log(fpn(i-64)+64d0)
 20   CONTINUE

ftl(1)=0d0
DO 30 i=2, 101
   a=i-1
   ftl(i)=ftl(i-1)+vscr(i)
 30   CONTINUE

RETURN
END
