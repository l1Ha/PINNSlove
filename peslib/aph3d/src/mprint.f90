SUBROUTINE mprint ( name, a, n, m, nmax )
USE Numeric_Kinds_Module
USE FileUnits_Module
!     <mprint> -- prints out a REAL(Kind=WP_Kind) matrix a(n,m), m<13
!     name -- 4 CHARACTER identification  (INTEGER), eg, 4hmtrx
!     a(nmax,m) -- array to be printed   (REAL(Kind=WP_Kind))
!     n -- number of rows to be printed
!     m -- number of columns to be printed, <=10 for effective results
!     nmax -- row dimension of a       (INTEGER)
CHARACTER(LEN=4)Name(1)
INTEGER n, m, nmax, i, j
REAL(Kind=WP_Kind) a(nmax,m)
 1000 FORMAT( /, a4, 5x, 10g12.4 )
 1001 FORMAT( 10x, 10g12.4 )
print 1000, name(1), (a(1,i), i = 1, m )
IF( n== 1 ) RETURN
DO 1 i = 2, n
print 1001, (a(i,j), j = 1, m )
    1 CONTINUE
RETURN
END
