SUBROUTINE lgrst(neig,ovlp,ilarg,jlarg,iskip,jskip)
USE Numeric_Kinds_Module
IMPLICIT NONE
INTEGER i,j, neig, ilarg, jlarg, iskip(neig), jskip(neig)
REAL(Kind=WP_Kind) ovlp(neig,neig), larg

larg=-10.d0
DO 2 i=1,neig
   IF(iskip(i)/=0)GOTO 2
   DO 1 j=1,neig
   IF(jskip(j)/=0)GOTO 1
   IF(ABS(ovlp(i,j))>larg)THEN
      ilarg=i
      jlarg=j
      larg=ABS(ovlp(i,j))
   ENDIF
   1  CONTINUE
2  CONTINUE

iskip(ilarg)=1
jskip(jlarg)=1

RETURN
ENDSUBROUTINE lgrst
