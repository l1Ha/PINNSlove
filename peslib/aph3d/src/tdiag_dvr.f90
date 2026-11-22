SUBROUTINE Tdiag_DVR ( v, e, t, subd, n, nm, nv )
USE Numeric_Kinds_Module
USE FileUnits_Module
!     <tdiag>--calls tred3, tql2, and trbak3 to diagonalize a matrix
!     v(n*(n+1)/2) -- matrix to be diagonalized, must be symmetric
!     e(n) -- output, the eigenvalues of <v>
!     t(n,n) -- output, the eigenvectors of <v>
!     n  -- appropriate dimension variables   (INTEGER*4)
!     subd(n) -- scratch array used by tred2 and tql2
!
   10 FORMAT ( ' diagonalization error #', i3, ' returned by tql2' )
INTEGER n, nm, nv, i, j, ierr
DIMENSION e(n), t(nm,n), subd(n), v(nv)
REAL(Kind=WP_Kind) v, e, t, subd
WRITE(Out_Unit,*)'TDIAG should be eliminated:, n,nv=',n,nv
CALL tred3(n,nv,v,e,subd,subd)
DO i=1,n
   DO j=1,n
      t(j,i)=0.d0
   ENDDO
   t(i,i)=1.d0
ENDDO
CALL tql2(nm,n,e,subd,t,ierr)

IF(ierr/=0)THEN
   WRITE(Out_Unit,*)'Error returned from tql2', ierr
   WRITE(Msg_Unit,*)'Error returned from tql2', ierr
   RETURN
ENDIF

CALL trbak3(nm,n,nv,v,n,t)

RETURN
ENDSUBROUTINE Tdiag_DVR

