SUBROUTINE sgemm_junk(m, n, l, a, ia, b, ib, c, ic, jtrpos, job)
!
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
!***beginprologue  sgemm_junk
!***datewritten            (yymmdd)
!***revisiondate  831207   (yymmdd)
!***categoryno.  d1b6
!***keywordsmatrix multiplication
!***authorjordan, tom, los alamos national laboratory
!***purposeperforms matrix multiplication with all forms
!  of transposition of the input and output matrices
!  and with four forms of accumulation.
!***description
!
!  given storage arrays a,b,c, this routine may
!  t
!  be used to compute ab=c, ab=c ,----,a b =c
!  (see the jtrpos parameter).  in addition, for each form
!  of transposition, the following four forms of accumulation
!  are possible:  ab=c, c+ab=c, -ab=c and c-ab=c (see
!  the job parameter).  much explicit transposition can be
!  avoided by using this routine and the transposition is
!  achieved at no cost.
!  this is a general routine for multiplying two arbitrary matrices
!  a and b with four forms of accumulation (see the job parameter)
!  and with any of the combination of transposes (see the jtrpos
!  parameter) including the result matrix c.
!
!  -arguments-
!
!  on input
!  m = the number of rows (columns) of a (a(transpose)) and
!  the number of rows (columns) of c (c(transpose)).
!  n = the number of colunms (rows) of a (a(transpose)) and the
!  the number of rows (columns) of b (b(transpose)).
!  l = the number of columns (rows) of b (b(transpose)) and the
!  number of columns (rows) of c (c(transpose)).
!  a, b, c = a 2-d array for a (a(transpose), etc.
!  ia, ib, ic = the first dimension of a, b, c in the calling
!  program.
!  jtrpos = a transposition indicator defined as follows:
!
!  value                transposition
!  0                  a   *b   =c
!  1                  a   *b   =c(t)
!  2                  a   *b(t)=c
!  3                  a   *b(t)=c(t)
!  4                  a(t)*b   =c
!  5                  a(t)*b   =c(t)
!  6                  a(t)*b(t)=c
!  7                  a(t)*b(t)=c(t)
!
!  job = a task designator
!  -1   c =  -ab
!  -2   c = c-ab
!  +1   c =  +ab
!  +2   c = c+ab
!
!  on output
!  c = a 2-d array containing the specified product
!  c (c(transpose)).
!
!***references(NONE)
!***routinescalled  sgemv_junk
!***specialimplementation  cftmath
!***endprologue  sgemm_junk
INTEGER m, n, l, ia, ib, ic, jtrpos, job, jb, isign, jump, j
REAL(Kind=WP_Kind) a(ia,n), b(ib,l), c(ic,l)
!***firstexecutable statement  sgemm_junk
IF (l <= 0) RETURN
jb = isign(iabs(job)+2*(jtrpos/4),job)
jump = jtrpos+1
!GOTO (10, 20, 30, 40, 50, 60, 70, 80),jump
IF(Jump==1)goto 10
IF(Jump==2)goto 20
IF(Jump==3)goto 30
IF(Jump==4)goto 40
IF(Jump==5)goto 50
IF(Jump==6)goto 60
IF(Jump==7)goto 70
IF(Jump==8)goto 80
 10   CONTINUE
!
!**   routine forms the matrix product c=ab, the routine calls a cal
!  kernel to multiply a matrix times a vector -
!input
!  m = the number of rows of a and c.
!  n = the number of columns of a and rows of b
!  l = the number of columns of b and c.
!  a,b,c = 2-d arrays containing the matrices a,b and c
!  ia,ib,ic= the first dimensions of arrays a,b and c
!output
!  c = the product matrix c=ab
!
DO  15    j=1, l
   CALL sgemv_junk (m, n, a, ia, b(1,j), 1, c(1,j), 1, jb)
 15   CONTINUE
RETURN
 50   CONTINUE
!
!**   this routine forms the matrix product a(transpose)b=c
!input
!  m = the number of columns of a and rows of c.
!  n = the number of rows of a and b.
!  l = the number of columns of b and c
!  a,b,c = 2-d arrays containing matrices a,b and c.
!  ia,ib,ic = the first dimensions of arrays a,b,c.
!output
!  c = the product matrix c=a(transpose)b
!
DO j=1, l
   CALL sgemv_junk (m, n, a, ia, b(1,j), 1, c(1,j), 1, jb)
ENDDO
RETURN
 30   CONTINUE
!
!cc   this routine forms the matrix product c=ab(transpose)
!
!input
!  m = the number of rows of a and c.
!  n = the number of columns of a and c.
!  l = the number of rows of b and columns of c.
!  a,b,c= 2-d arrays containing the matrices a,b and c.
!  ia,ib,ic= the first dimension of arrays a,b,c.
!output
!  c = the product matrix c = ab(transpose)
!
DO 35    j=1, l
   CALL sgemv_junk (m, n, a, ia, b(j,1), ib, c(1,j), 1, jb)
 35   CONTINUE
RETURN
 20   CONTINUE
!
!cc   this routine forms the matrix product ab and stores its transpose
!  in c.
!
!input
!  m = the number of rows of a and columns of c.
!  n = the number of columns of a and rows of b.
!  l = the number of columns of b and rows of c.
!  a,b,c = 2-d arrays containing matrices a,b and c.
!  ia,ib,ic = the first dimension of arrays a,b,c.
!output
!  c  contains the transpose of ab.
!
DO 25    j=1, l
   CALL sgemv_junk (m, n, a, ia, b(1,j), 1, c(j,1), ic, jb)
 25   CONTINUE
RETURN
 70   CONTINUE
!
!cc   this routine forms the matrix product c=a(transpose)b(transpose)
!
!input
!  m = the number of columns of a and rows of c.
!  n = the number of rows of a and columns of b.
!  l = the number of rows of b and columns of c
!  a,b,c = 2-d arrays containing matrices a,b and c.
!  ia,ib,ic= the first dimensions of array a,b,c.
!output
!  c = the product matrix c=a(transpose)b(transpose)
!
DO 75    j=1, l
   CALL sgemv_junk (m, n, a, ia, b(j,1), ib, c(1,j), 1, jb)
 75   CONTINUE
RETURN
 60   CONTINUE
!
!**   this routine forms the matrix product a(transpose)b and stores
!  its transpose in c.
!input
!  m = the number of columns of a and c.
!  n = the number of rows of a and b.
!  l = the number of columns of b and rows of c.
!  a,b,c = 2-d arrays containing matrices a,b and c.
!  ia,ib,ic = the first dimensions of arrays a,b,c.
!output
!  c contains the transpose of the product matrix a(transpose)b
!
DO j=1, l
   CALL sgemv_junk (m, n, a, ia, b(1,j), 1, c(j,1), ic, jb)
ENDDO
RETURN
 40   CONTINUE
!
!**   this routine forms the matrix product ab(transpose) and stores
!  the its transpose in c
!input
!  m = the number of rows of a and columns of c.
!  n = the number of columns of a and b.
!  l = the number of rows of b and c
!  a,b,c = 2-d arrays containing matrices a,b and c.
!  ia,ib,ic = the first dimensions of arrays a,b,c.
!output
!  c contains the transpose of the product ab(transpose)
!
DO j=1, l
   CALL sgemv_junk (m, n, a, ia, b(j,1), ib, c(j,1), ic, jb)
ENDDO
RETURN
 80   CONTINUE
!
!**   this routine forms the matrix product a(transpose)b(transpose) and
!  stores its transpose in c.
!input
!  m = the number of columns of a and c.
!  n = the number of rows of a and columns of b.
!  l = the number of rows of b and c.
!  a,b,c = 2-d arrays containing matrices a,b and c.
!  ia,ib,ic = the first dimensions of arrays a,b,c.
!output
!  c contains the transpose of the product a(transpose)b(transpose)
DO j=1, l
   CALL sgemv_junk (m, n, a, ia, b(j,1), ib, c(j,1), ic, jb)
ENDDO
RETURN
END
