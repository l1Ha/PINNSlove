SUBROUTINE tdiag (v, e, t, subd, n)
USE Numeric_Kinds_Module
USE FileUnits_Module
USE Numbers_Module
!     <tdiag>--calls tred3, tql2, and trbak3 to diagonalize a matrix
!     v(n*(n+1)/2) -- matrix to be diagonalized, must be symmetric
!     e(n) -- output, the eigenvalues of <v>
!     t(n,n) -- output, the eigenvectors of <v>
!     n  -- appropriate dimension variables   (INTEGER)
!     subd(n) -- scratch array used by tred2 and tql2
!
IMPLICIT NONE
INTEGER n, nv, i, j, ierr
EXTERNAL tred3, tql2, trbak3
REAL(Kind=WP_Kind) e(n), t(n,n), subd(n), v(1)

WRITE(Out_Unit,'("TDIAG and the routine that it calls should be eliminated by calling a LAPACK routine")')
nv=(n*(n+1))/2
CALL tred3(n,nv,v,e,subd,subd)

DO i=1,n
   DO j=1,n
      t(j,i)=zero
   ENDDO
   t(i,i)=one
ENDDO

CALL tql2 (n, n, e, subd, t, ierr)

IF(ierr==0)THEN
    CALL trbak3 (n, n, nv, v, n, t)
    RETURN
ELSE
    WRITE(Out_Unit,40)ierr
    STOP 'tdiag'
ENDIF

!----------------***end-tdiag***---------------------------------------
RETURN
   40 FORMAT( ' diagonalization error #', i3, ' returned by tql2' )
END
