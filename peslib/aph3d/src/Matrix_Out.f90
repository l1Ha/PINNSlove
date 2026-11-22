SUBROUTINE Matrix_Out (Matrix, NRows, NCols, MatrixName, MatrixHeading, Print_Flag)
USE Numeric_Kinds_Module, ONLY: IW_Kind, WP_Kind     ! Numeric Kinds
USE FileUnits_Asymptotic_Module, ONLY: Out_Unit, Msg_Unit         ! Contains fileunit numbers
USE MatrixLabels_Module, ONLY: LabelRows, LabelColumns, NQOpen, NQLabels, QStates, RowNames, ColumnNames      ! Use matrixlabels module for printing quantum labels
IMPLICIT NONE
! General Matrix printing routine.
!
!     Required Input <=====
!        Matrix     Matrix of dimension (NRows,NCols) to be printed
!        NRows      Number of rows in Matrix
!        NCols      Number of columns in Matrix
!        MatrixName Actual name of the matrix to be printed
!        MatrixHeading  Short Description of the matrix
!        Print_Flag Controls the amount of printed output
!           Print_Flag='None'    Nothing is printed
!           Print_Flag='Little'  The first row will be printed
!           Print_Flag='Medium'  A 4-by4 matrix submatrix will be printed.
!           Print_Flag='Full'    The Full Matrix will be printed
!
! This routine is called by:
!      All routines needing Matrix output.
INTEGER(KIND=IW_Kind)               :: Row, Col, NRowP, NColP, ColB, ColE, IthQuant
INTEGER(KIND=IW_Kind), INTENT(IN)   :: NRows, NCols
INTEGER(KIND=IW_Kind), PARAMETER    :: PrintCol=9
CHARACTER(LEN=3),      PARAMETER    :: ColN='Col'
CHARACTER(LEN=6),      INTENT(IN)   :: Print_Flag
CHARACTER(LEN=*),      INTENT(IN)   :: MatrixHeading, MatrixName
CHARACTER(LEN=12)                   :: Name
CHARACTER(LEN=12),     PARAMETER    :: Blnk6 ="      "
CHARACTER(LEN=12),     PARAMETER    :: Blnk12="            "
CHARACTER(LEN=15)                   :: LabelValues(PrintCol)
CHARACTER(LEN=21),     PARAMETER    :: ProcName='Matrix_Out'
REAL(KIND=WP_Kind),    INTENT(IN)   :: Matrix(NRows,NCols)

CALL PoptAsy(ProcName, Print_Flag)
IF(Print_Flag=='None')RETURN                                 ! Return if matrix printing is not desired

IF(TRIM(MatrixName)=="")THEN
   WRITE(Out_Unit,*)"Warning MatrixName not defined"
   WRITE(Msg_Unit,*)"Warning MatrixName not defined"
ELSEIF(TRIM(MatrixHeading)=="")THEN
   WRITE(Out_Unit,*)"Warning MatrixHeading not defined"
   WRITE(Msg_Unit,*)"Warning MatrixHeading not defined"
ELSEIF(LEN_TRIM(Print_Flag)==0)THEN
   WRITE(Out_Unit,*)"Warning PrintFlag not defined"
   WRITE(Msg_Unit,*)"Warning PrintFlag not defined"
ELSEIF(NRows<=0.or.NCols<0)THEN
   WRITE(Out_Unit,*)"ERROR: NRows<=0.or.NCols<0", NRows, NCols
   WRITE(Msg_Unit,*)"ERROR: NRows<=0.or.NCols<0", NRows, NCols
   STOP "Error in Matrix_Out"
ENDIF

Name=Blnk12
Name=MatrixName(1:Min(Len_Trim(MatrixName),12))              ! Trim Matrix name if it is too long

WRITE(Out_Unit,*)
WRITE(Out_Unit,*)MatrixHeading                               ! Print Heading
WRITE(Out_Unit,*)Name,' is a',NRows,' by',NCols,' matrix.'   ! Print Number of rows and columns

!
! Determine number of rows and columns of matrix to print
!
IF(Print_Flag=='Medium')THEN
   NRowP=Min(PrintCol,NRows)
   NColP=Min(PrintCol,NCols)
   IF(NCols>PrintCol.or.NRows>PrintCol)THEN
      WRITE(Out_Unit,*)'Printing ',NRowP, ' Rows and ',NColP, ' Columns of ', Name
   ENDIF
ELSEIF(Print_Flag=='Little')THEN
   NRowP=1
   NColP=Min(PrintCol,NCols)
ELSE
   NRowP=NRows
   NColP=NCols
ENDIF

!
!  Print matrix
!
Row=1
DO ColB=1,NColP,PrintCol
   ColE=Min(ColB+PrintCol-1,NColP)
   DO Col=ColB,ColE
      IF(LabelColumns)THEN
          WRITE(LabelValues(Col+1-ColB),'(5I3)')(QStates(IthQuant,Col),IthQuant=1,NQlabels)            ! Get column values
      ENDIF
   ENDDO
   IF(LabelRows.AND.LabelColumns)THEN
      WRITE(Out_Unit,'(15x, 9(5x, A15, 3x))')(ColumnNames,Col=ColB,ColE)                              ! Print column names
      WRITE(Out_Unit,'(A15,9(5x, A15, 3x))')RowNames, (LabelValues(Col+1-ColB),Col=ColB,ColE)
   ELSEIF(LabelRows.AND.(.NOT.LabelColumns))THEN
      WRITE(Out_Unit,'(A15, 9(9x, A3,I3, 8x))')RowNames, (Coln,Col,Col=ColB,ColE)
   ELSEIF((.NOT.LabelRows).AND.LabelColumns)THEN
      WRITE(Out_Unit,'(15x, 9(5x, A15, 3x))')(ColumnNames,Col=ColB,ColE)                              ! Print column names
      WRITE(Out_Unit,'(15x, 9(5x, A15, 3x))')(LabelValues(Col+1-ColB),Col=ColB,ColE)
   ELSEIF((.NOT.LabelRows).AND.(.NOT.LabelColumns))THEN
      WRITE(Out_Unit,'(15x, 9(9x, A3,I3, 8x))')(Coln,Col,Col=ColB,ColE)
   ENDIF
   DO Row = 1, NRowP
      IF(LabelRows)THEN
         WRITE(LabelValues(1),'(5I3)')(QStates(IthQuant,Row),IthQuant=1,NQlabels)                     ! Get row rabel
         WRITE(Out_Unit,'(A15, 9(ES23.15))')LabelValues(1), (Matrix(Row,Col),Col=ColB,ColE)           ! Print row with label
      ELSE
         WRITE(Out_Unit,'(1x, 5x,"Row ", I3, 2x, 9(ES23.15))')Row, (Matrix(Row,Col),Col=ColB,ColE)    ! Print row without label
      ENDIF
   ENDDO
ENDDO

RETURN
ENDSUBROUTINE Matrix_Out
