SUBROUTINE MxOutL(a, nrows, ncols, ms, mrows, mcolmn)
USE FileUnits_Module
USE Narran_Module
USE nstate_Module
USE qlam_Module
USE Qall_Module
!
!-----------------------------------------------------------------------
!  matrix printing routine which labels the columns and rows.
!  The matrix can have:
!  1. Space-fixed rows or columns.
!  2. Projectile body-fixed rows or columns.
!  3. Molecular body-fixed rows or columns.
!  4. APH body-fixed rows or columns.
!
!  INPUT PARAMETERS.
!  a       matrix to be printed.
!  nrows   number of rows in matrix a.
!  ncols   number of columns in matrix a.
!  ms      storage option.
!  ms=0  the matrix is stored as a general matrix.
!  ms=1  the matrix is stored as a symmetric matrix.
!  mrows   label for the rows.
!  mrows='space' space-fixed basis for rows.
!  mrows='probf' projectile body-fixed basis for rows.
!  mrows='molbf' molecular body-fixed basis for rows.
!  mrows='aph  ' APH body-fixed basis for rows.
!  mcolmn  label for the columns.
!  mcols='space' space-fixed basis for columns.
!  mcols='probf' projectile body-fixed basis for columns.
!  mcols='molbf' molecular body-fixed basis for columns.
!  mcols='aph  ' APH body-fixed basis for columns.
!
!  OUTPUT PARAMETERS.
!  NONE: This routine does not modify any input parameters.
!
!-----------------------------------------------------------------------
IMPLICIT NONE
INTEGER j1, j2, i, j , nrows, ncols, ms, index, ind
CHARACTER(LEN=1) alabel(narran)
CHARACTER(LEN=5) mrows, mcolmn
CHARACTER(LEN=12) lspace, lprobf, lmolbf, laph
REAL(Kind=WP_Kind) a

!  D I M E N S I O N S

DIMENSION a(1), index(1:10)
DATA lspace /'ch  v  j  l '/
DATA lprobf /'ch  v  j meg'/
DATA lmolbf /'ch  v  l lam'/
DATA   laph /'  t  lam    '/
DATA alabel/'a', 'b', 'c'/
!-----------------------------------------------------------------------
!  Print the matrix with space-fixed labels for the columns.
!  The matrix can have:
!  1. Space-fixed rows.
!  2. Projectile body-fixed rows.
!  3. Molecular body-fixed rows.
!  4. APH body-fixed rows.
!-----------------------------------------------------------------------
IF(mcolmn == 'space')THEN
   i=1
   j1=1
   j2=10
   IF(j2 > ncols) j2=ncols
   DO 11 j1=1, ncols, 10
      j2=j1+9
      IF(j2 > ncols) j2=ncols
      WRITE(Out_Unit,90) j1, j2
      WRITE(Out_Unit,60) (lspace, j=j1, j2)
!-----------------------------------------------------------------------
!  The rows of this matrix has space-fixed labels.
!-----------------------------------------------------------------------
      IF(mrows == 'space')THEN
         WRITE(Out_Unit,20)lspace, (alabel(kchan(j)), mvib3(j), jrot3(j), lorb3(j), j=j1, j2)
         DO 21 i=1, nrows
            CALL indxer(index, ms, i, j1, nrows, ncols)
            WRITE(Out_Unit,70) alabel(kchan(i)), mvib3(i), jrot3(i), lorb3(i), (a(index(ind)), ind=1, j2-j1+1)
 21            CONTINUE
!-----------------------------------------------------------------------
!  The rows of this matrix has projectile body-fixed labels.
!-----------------------------------------------------------------------
      ELSEIF(mrows == 'probf')THEN
         WRITE(Out_Unit,20)lprobf, (alabel(kchan(j)), mvib3(j),jrot3(j), lorb3(j), j=j1, j2)
         DO 31 i=1, nrows
            CALL indxer(index, ms, i, j1, nrows, ncols)
            WRITE(Out_Unit,70) alabel(kchan(i)), mvib3(i), jrot3(i), mega3(i), (a(index(ind)), ind=1, j2-j1+1)
 31            CONTINUE
!-----------------------------------------------------------------------
!  The rows of this matrix has molecular body-fixed labels.
!-----------------------------------------------------------------------
      ELSEIF(mrows == 'molbf')THEN
         WRITE(Out_Unit,20)lmolbf, (alabel(kchan(j)), mvib3(j),jrot3(j), lorb3(j), j=j1, j2)
         DO i=1, nrows
            CALL indxer(index, ms, i, j1, nrows, ncols)
            WRITE(Out_Unit,70) alabel(kchan(i)),mvib3(i), lorb3(i), mega3(i), (a(index(ind)), ind=1, j2-j1+1)
         ENDDO
!-----------------------------------------------------------------------
!  The rows of this matrix has APH body-fixed labels.
!-----------------------------------------------------------------------
      ELSEIF(mrows == 'aph  ')THEN
         WRITE(Out_Unit,30)laph, (alabel(kchan(j)), mvib3(j),jrot3(j), mega3(j), j=j1, j2)
         DO 51 i=1, nrows
            WRITE(Out_Unit,80) tee(i), lambda(i), (a(index(ind)), ind=1, j2-j1+1)
 51            CONTINUE
      ENDIF
 11      CONTINUE
   RETURN
ENDIF
!-----------------------------------------------------------------------
!  Print the matrix with projectile body-fixed labels for the columns.
!  The matrix can have:
!  1. Space-fixed rows.
!  2. Projectile body-fixed rows.
!  3. Molecular body-fixed rows.
!  4. APH body-fixed rows.
!-----------------------------------------------------------------------
IF(mcolmn == 'probf')THEN
   i=1
   j1=1
   j2=10
   IF(j2 > ncols) j2=ncols
   DO 12 j1=1, ncols, 10
      j2=j1+9
      IF(j2 > ncols) j2=ncols
      WRITE(Out_Unit,90) j1, j2
      WRITE(Out_Unit,60) (lprobf, j=j1, j2)
!-----------------------------------------------------------------------
!  The rows of this matrix has space-fixed labels.
!-----------------------------------------------------------------------
      IF(mrows == 'space')THEN
         WRITE(Out_Unit,20)lspace, (alabel(kchan(j)), mvib3(j), jrot3(j), lorb3(j), j=j1, j2)
         DO 22 i=1, nrows
            CALL indxer(index, ms, i, j1, nrows, ncols)
            WRITE(Out_Unit,70) alabel(kchan(i)), mvib3(i), jrot3(i), lorb3(i),(a(index(ind)), ind=1, j2-j1+1)
 22            CONTINUE
!-----------------------------------------------------------------------
!  The rows of this matrix has projectile body-fixed labels.
!-----------------------------------------------------------------------
      ELSEIF(mrows == 'probf')THEN
         WRITE(Out_Unit,20)lprobf, (alabel(kchan(j)), mvib3(j),jrot3(j), lorb3(j), j=j1, j2)
         DO 32 i=1, nrows
            CALL indxer(index, ms, i, j1, nrows, ncols)
            WRITE(Out_Unit,70) alabel(kchan(i)),mvib3(i), jrot3(i), mega3(i),(a(index(ind)), ind=1, j2-j1+1)
 32            CONTINUE
!-----------------------------------------------------------------------
!  The rows of this matrix has molecular body-fixed labels.
!-----------------------------------------------------------------------
      ELSEIF(mrows == 'molbf')THEN
         WRITE(Out_Unit,20)lmolbf, (alabel(kchan(j)), mvib3(j),jrot3(j), lorb3(j), j=j1, j2)
         DO i=1, nrows
            CALL indxer(index, ms, i, j1, nrows, ncols)
            WRITE(Out_Unit,70) alabel(kchan(i)),mvib3(i), lorb3(i), mega3(i), (a(index(ind)), ind=1, j2-j1+1)
         ENDDO
!-----------------------------------------------------------------------
!  The rows of this matrix has APH body-fixed labels.
!-----------------------------------------------------------------------
      ELSEIF(mrows == 'aph  ')THEN
         WRITE(Out_Unit,30)laph, (alabel(kchan(j)), mvib3(j), jrot3(j), mega3(j), j=j1, j2)
         DO 52 i=1, nrows
            WRITE(Out_Unit,80) tee(i), lambda(i), (a(index(ind)), ind=1, j2-j1+1)
 52            CONTINUE
      ENDIF
 12      CONTINUE
   RETURN
ENDIF
!-----------------------------------------------------------------------
!  Print the matrix with molecule body-fixed labels for the columns.
!  The matrix can have:
!  1. Space-fixed rows.
!  2. Projectile body-fixed rows.
!  3. Molecular body-fixed rows.
!  4. APH body-fixed rows.
!-----------------------------------------------------------------------
IF(mcolmn == 'molbf')THEN
   i=1
   j1=1
   j2=10
   IF(j2 > ncols) j2=ncols
   DO 13 j1=1, ncols, 10
      j2=j1+9
      IF(j2 > ncols) j2=ncols
      WRITE(Out_Unit,90) j1, j2
      WRITE(Out_Unit,60) (lmolbf, j=j1, j2)
!-----------------------------------------------------------------------
!  The rows of this matrix has space-fixed labels.
!-----------------------------------------------------------------------
      IF(mrows == 'space')THEN
         WRITE(Out_Unit,20)lspace, (alabel(kchan(j)), mvib3(j),jrot3(j), lorb3(j), j=j1, j2)
         DO 23 i=1, nrows
            CALL indxer(index, ms, i, j1, nrows, ncols)
            WRITE(Out_Unit,70) alabel(kchan(i)), mvib3(i), jrot3(i), lorb3(i), (a(index(ind)), ind=1, j2-j1+1)
 23            CONTINUE
!-----------------------------------------------------------------------
!  The rows of this matrix has projectile body-fixed labels.
!-----------------------------------------------------------------------
      ELSEIF(mrows == 'probf')THEN
         WRITE(Out_Unit,20)lprobf, (alabel(kchan(j)), mvib3(j),jrot3(j), lorb3(j), j=j1, j2)
         DO 33 i=1, nrows
            CALL indxer(index, ms, i, j1, nrows, ncols)
            WRITE(Out_Unit,70) alabel(kchan(i)), mvib3(i), jrot3(i), mega3(i), (a(index(ind)), ind=1, j2-j1+1)
 33            CONTINUE
!-----------------------------------------------------------------------
!  The rows of this matrix has molecular body-fixed labels.
!-----------------------------------------------------------------------
      ELSEIF(mrows == 'molbf')THEN
         WRITE(Out_Unit,20)lmolbf, (alabel(kchan(j)), mvib3(j),jrot3(j), lorb3(j), j=j1, j2)
         DO i=1, nrows
            CALL indxer(index, ms, i, j1, nrows, ncols)
            WRITE(Out_Unit,70) alabel(kchan(i)),mvib3(i), lorb3(i), mega3(i), (a(index(ind)), ind=1, j2-j1+1)
         ENDDO
!-----------------------------------------------------------------------
!  The rows of this matrix has APH body-fixed labels.
!-----------------------------------------------------------------------
      ELSEIF(mrows == 'aph  ')THEN
         WRITE(Out_Unit,30)laph, (alabel(kchan(j)), mvib3(j),jrot3(j), mega3(j), j=j1, j2)
         DO i=1, nrows
            WRITE(Out_Unit,80) tee(i), lambda(i), (a(index(ind)), ind=1, j2-j1+1)
         ENDDO
      ENDIF
 13      CONTINUE
   RETURN
ENDIF
!-----------------------------------------------------------------------
!  Print the matrix with APH labels for the columns.
!  The matrix can have:
!  1. Space-fixed rows.
!  2. Projectile body-fixed rows.
!  3. Molecular body-fixed rows.
!  4. APH body-fixed rows.
!-----------------------------------------------------------------------
IF(mcolmn == 'aph  ')THEN
   i=1
   j1=1
   j2=10
   IF(j2 > ncols) j2=ncols
   DO 14 j1=1, ncols, 10
      j2=j1+9
      IF(j2 > ncols) j2=ncols
      WRITE(Out_Unit,90) j1, j2
      WRITE(Out_Unit,60) (laph, j=j1, j2)
!-----------------------------------------------------------------------
!  The rows of this matrix has space-fixed labels.
!-----------------------------------------------------------------------
      IF(mrows == 'space')THEN
         WRITE(Out_Unit,50)lspace, (tee(j), lambda(j), j=j1, j2)
         DO 24 i=1, nrows
            CALL indxer(index, ms, i, j1, nrows, ncols)
            WRITE(Out_Unit,70) alabel(kchan(i)), mvib3(i), jrot3(i), lorb3(i), (a(index(ind)), ind=1, j2-j1+1)
 24            CONTINUE
!-----------------------------------------------------------------------
!  The rows of this matrix has projectile body-fixed labels.
!-----------------------------------------------------------------------
      ELSEIF(mrows == 'probf')THEN
         WRITE(Out_Unit,50)lprobf,(tee(j), lambda(j), j=j1, j2)
         DO 34 i=1, nrows
            CALL indxer(index, ms, i, j1, nrows, ncols)
            WRITE(Out_Unit,70) alabel(kchan(i)), mvib3(i), jrot3(i), mega3(i), (a(index(ind)), ind=1, j2-j1+1)
 34            CONTINUE
!-----------------------------------------------------------------------
!  The rows of this matrix has molecular body-fixed labels.
!-----------------------------------------------------------------------
      ELSEIF(mrows == 'molbf')THEN
         WRITE(Out_Unit,50)lmolbf, (tee(j), lambda(j), j=j1, j2)
         DO i=1, nrows
            CALL indxer(index, ms, i, j1, nrows, ncols)
            WRITE(Out_Unit,70) alabel(kchan(i)), mvib3(i), lorb3(i), mega3(i), (a(index(ind)), ind=1, j2-j1+1)
         ENDDO
!-----------------------------------------------------------------------
!  The rows of this matrix has APH body-fixed labels.
!-----------------------------------------------------------------------
      ELSEIF(mrows == 'aph  ')THEN
         WRITE(Out_Unit,50)laph,(tee(j), lambda(j), j=j1, j2)
         DO 54 i=1, nrows
            CALL indxer(index, ms, i, j1, nrows, ncols)
            WRITE(Out_Unit,80) tee(i), lambda(i), (a(index(ind)), ind=1, j2-j1+1)
 54            CONTINUE
      ENDIF
 14      CONTINUE
   RETURN
ENDIF
!-----------------------------------------------------------------------
!  STOP execution IF column label is not correct.
!-----------------------------------------------------------------------
WRITE(Out_Unit,*)'Execution was stopped in routine MXOUTL'
WRITE(Out_Unit,*)'Matrix column label not valid'
STOP 'mxoutl'
!----------------***END-mxoutl***---------------------------------------
!
 20   FORMAT(1x, a12, 10(1x, a1, 3i3, 1x))
 30   FORMAT(1x, a12, 10(1x, a1, 3i3, 1x))
 40   FORMAT(1x, a12, 10(1x, a1, 3i3, 1x))
 50   FORMAT(1x, a12, 10(1x, i3, 1x, i3, 4x))
 60   FORMAT(/, 12x, 10a12)
 70   FORMAT(1x, 1x, a1, 3i3, 10(1pe12.4))
 80   FORMAT(1x, 1x, 3x, 2i3, 10(1pe12.4))
 90   FORMAT(/, '  columns ', i4, ' through ', i4)
END
