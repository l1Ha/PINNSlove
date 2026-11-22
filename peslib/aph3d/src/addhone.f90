SUBROUTINE addhone ( n, a, lda, h)
USE FileUnits_Module
!-----------------------------------------------------------------------
!   This routine adds the identity matrix to matrix a.
!-----------------------------------------------------------------------
!trp>>
!    9/15/1997 TRP added lda argument so that matrix a may
!                  be any size.
!<<trp
IMPLICIT NONE
INTEGER n, i, lda
!trp>>
!   the following line modified as noted above 9/15/97.
!     DIMENSION a(n,n)
REAL(Kind=WP_Kind)  h, a(lda,*)
!<<trp


!-----------------------------------------------------------------------
!  Add h*Identity matrix.
!-----------------------------------------------------------------------
DO i = 1, n
    a(i,i) = a(i,i)+h
ENDDO

RETURN
ENDSUBROUTINE addhone
