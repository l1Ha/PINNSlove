MODULE MatrixLabels_Module
   LOGICAL LabelRows, LabelColumns          ! Set to true IF quantum labels are to be printed
   INTEGER NQOpen, NQLabels                 ! First dimension of QStates and number of quantum labels
   INTEGER, ALLOCATABLE:: QStates(:,:)      ! Quantum numbers
   CHARACTER(LEN=15) RowNames              ! Identification labels for Rows
   CHARACTER(LEN=15) ColumnNames           ! Identification labels for Columns
END MODULE MatrixLabels_Module
