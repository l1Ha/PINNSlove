SUBROUTINE tdiagrw ( v, e, t, subd, n, nmax, nvmax )
USE Numeric_Kinds_Module
USE FileUnits_Module
!     <tdiag>--calls tred3, tql2, and trbak3 to diagonalize a matrix
!     v(n*(n+1)/2) -- matrix to be diagonalized, must be symmetric
!     e(n) -- output, the eigenvalues of <v>
!     t(n,n) -- output, the eigenvectors of <v>
!     n  -- appropriate dimension variables   (INTEGER)
!     subd(n) -- scratch array used by tred2 and tql2
!
   10 FORMAT( ' diagonalization error #', i3, ' returned by tql2' )
INTEGER I, IErr, J, N, Nmax, Nv, NvMax
REAL(Kind=WP_Kind) e(nmax), t(nmax,nmax), subd(nmax), v(nvmax)
nv=(n*(n+1))/2
CALL tred3rw(n,nv,v,e,subd,subd,nmax,nvmax)
DO 30 i=1,n
DO 20 j=1,n
t(j,i)=0.d0
   20 CONTINUE
t(i,i)=1.d0
   30 CONTINUE
CALL tql2rw(n,n,e,subd,t,ierr,nmax,nvmax)
IF(ierr/=0) GOTO 40
CALL trbak3rw(n,n,nv,v,n,t,nmax,nvmax)
RETURN
   40 CONTINUE
IF( ierr /= 0 ) WRITE(Out_Unit,10) ierr
RETURN
!  --------------***END-tdiag***---------------------------------------
END
