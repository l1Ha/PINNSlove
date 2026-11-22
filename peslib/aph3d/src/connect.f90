SUBROUTINE connect(neig,nrho,jrho,ovlp,track)
USE Numeric_Kinds_Module
IMPLICIT NONE
INTEGER jrho,i,neig,nrho,track(neig,NRho)
INTEGER ilarg, jlarg, iskip(10000), jskip(10000)
REAL(Kind=WP_Kind) ovlp(neig,neig)

iskip=0
jskip=0

DO i=1,neig
   CALL lgrst(neig,ovlp,ilarg,jlarg,iskip,jskip)
   track(jlarg,jrho)=track(ilarg,jrho-1)
ENDDO

RETURN
ENDSUBROUTINE connect
