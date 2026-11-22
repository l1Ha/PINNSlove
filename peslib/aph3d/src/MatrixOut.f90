SUBROUTINE MatrixOut(Matrix, NRows, NCols, MatrixName, MatrixHeading,  &
                       RowLbl, ColLbl, LabelFlag, LabelNames, FullMat)
USE Numeric_Kinds_Module
USE Fileunits_Module
IMPLICIT NONE
LOGICAL LabelFlag, FullMat
INTEGER NRows, NCols, Row, Col, NRowP, NColP, ColB, ColE
CHARACTER(LEN=12)  Name
CHARACTER(LEN=12), PARAMETER:: Blnk12='            '
CHARACTER(LEN=*) MatrixHeading
CHARACTER(LEN=*) MatrixName
CHARACTER(LEN=*) LabelNames
CHARACTER(LEN=11)  RowLbl(NRows), ColLbl(NCols), Lbl
CHARACTER(LEN=12), PARAMETER:: Blnk11='          '
REAL(Kind=WP_Kind) Matrix(NRows,NCols)

Name=Blnk12
Name=MatrixName(1:Min(Len_Trim(MatrixName),12))
Lbl=Blnk11
Lbl=Labelnames(1:Min(Len_Trim(LabelNames),11))

WRITE(Out_Unit,*)MatrixHeading
WRITE(Out_Unit,*)Name,' is a',NRows,' by',NCols,' matrix.'

IF(.NOT.FullMat)THEN
   IF(LabelFlag)THEN
      NRowP=Min(8,NRows)
      NColP=Min(8,NCols)
   ELSE
      NRowP=Min(9,NRows)
      NColP=Min(9,NCols)
   ENDIF
   IF(NCols>9.or.NRows>9)THEN
      WRITE(Out_Unit,*)'Printing ',NRowP,' Rows and ',NColP, ' Columns of matrix ', Name
   ENDIF
ELSE
    NRowP=NRows
    NColP=NCols
ENDIF

Row=1
IF(LabelFlag)THEN
   DO ColB=1,NColP,8
      ColE=Min(ColB+7,NColP)
      Row=1
      WRITE(Out_Unit,250)Name, (Col,Col=ColB,ColE)
      WRITE(Out_Unit,260)Lbl, (ColLbl(Col),Col=ColB,ColE)
      WRITE(Out_Unit,270)Row, RowLbl(Row),(Matrix(1,Col),Col=ColB,ColE)
      DO Row = 2, NRowP
          WRITE(Out_Unit,270)Row, RowLbl(Row),(Matrix(Row,Col),Col=ColB,ColE)
      ENDDO
   ENDDO
ELSE
   DO ColB=1,NColP,9
      ColE=Min(ColB+8,NColP)
      Row=1
      WRITE(Out_Unit,230)Name, (Col,Col=ColB,ColE)
      WRITE(Out_Unit,240)Row, (Matrix(1,Col),Col=ColB,ColE)
      DO Row = 2, NRowP
          WRITE(Out_Unit,240)Row, (Matrix(Row,Col),Col=ColB,ColE)
      ENDDO
   ENDDO
ENDIF
!
!
230  FORMAT(1x, A12,2x,9(3x, I3, 3x, 'Col'))
240  FORMAT(1x, 5x,'Row ', I3, 2x, 9(ES12.4))
250  FORMAT(1x, A12,9x,8(3x, I3, 3x, 'Col'))
260  FORMAT(1x, 8x, A11,3x,8(A11,1x))
270  FORMAT(1x, 'Row ', I3, 1x, A11,2x, 8(ES12.4))
END
