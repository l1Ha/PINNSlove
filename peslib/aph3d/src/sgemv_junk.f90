SUBROUTINE sgemv_junk(m, n, a, ia, x, ix, y, iy, job)
USE Numeric_Kinds_Module
IMPLICIT NONE
!***beginprologue  sgemv_junk
!***datewritten            (yymmdd)
!***revisiondate  831207   (yymmdd)
!***categoryno.  d1b6
!***keywordsmatrix multiplication,vector
!***authorjordan, tom, los alamos national laboratory
!***purposeperforms multiplication of a matrix a times a vector x or a
!  transposed times x.  four forms of accumulation are
!  provided.
!***description
!
!  given storage array a, this routine can be
!  t
!  used to compute a x=y or ax=y.  in addition, four forms
!  of accumulation are provided for each form of transposition.
!  these are ax=y, y+ax=y, -ax=y and y-ax=y.  this routine is
!  designed to provide optimal efficiency on the cray-1 and
!  near optimal efficiency on the cdc 7600.
!
!  -arguments-
!
!  on input
!  m = the number of rows of array a, IF matrix a is to be used.
!  = the number of columns of array a, IF a(transpose)is used.
!
!  n = the number of columns of array a, IF matrix a is used.
!  = the number of rows of array a, IF a(transpose) is used.
!  = the number of elements of the vector x.
!
!  a = the input matrix a, dimensioned a(ia,n).
!  ia = the leading dimension of a.
!  x = the input vector.
!  ix = the increment between elements of x.
!  y = the product vector.  (output)
!  iy = the spacing between elements of the result vector.
!  job = +1  y =   + a*x
!  +2  y = y + a*x
!  +3  y =   + a(transpose)*x
!  +4  y = y + a(transpose)*x
!
!  IF the job is negative, THEN subtraction is performed instead
!  of addition.
!
!  on output
!  y = the result vector as defined by job.  y may be a vector
!  with arbitrary but uniform spacing, iy, of elements.
!
!***references(NONE)
!***routinescalled  saxpy_junk
!***endprologue  sgemv_junk
INTEGER m, n, ia, ix, iy, job, ii, j, ij
REAL(Kind=WP_Kind) a(ia*n*n,1), x(ix,n*n), y(iy,n*n)
INTRINSIC iabs, mod
!***firstexecutable statement  sgemv_junk
IF(n <= 0) RETURN
ii = 1
ij = ia
IF(((iabs(job)-1)/2) == 0) GOTO 1
ij = 1
ii = ia
    1 CONTINUE
IF(mod(iabs(job)-1,2) /= 0) GOTO 4
!
!   set y = 0.
!
DO 2 j=1, m
   y(1,j) = 0.d0
    2 CONTINUE
    4 CONTINUE
!
!  accumulate columns of a to y
!
IF(job < 0) GOTO 5
DO 3 j=1, n
   CALL saxpy_junk (m, x(1,j), a(1+(j-1)*ij,1), ii, y(1,1), iy)
    3 CONTINUE
RETURN
    5 CONTINUE
DO  6  j=1, n
   CALL saxpy_junk (m, -x(1,j), a(1+(j-1)*ij,1), ii, y(1,1), iy)
    6 CONTINUE
RETURN
END
