MODULE POptAsy_Module
   IMPLICIT NONE
   SAVE
!-----------------------------------------------------------------------
!  subs  List of SUBROUTINE names with debugging options turned on.
!-----------------------------------------------------------------------
   INTEGER, PARAMETER:: NSubs=30
   LOGICAL AllNone/.False./, AllLittle/.False./, AllMedium/.False./, AllFull/.False./
   CHARACTER(LEN=21) Subs(NSubs)
   CHARACTER(LEN=6) Print_Level(NSubs)
   DATA  Subs /'ASymmetry',            &
               'Asymptotic',           &
               'Asymptotic_Analysis',  &
               'Asymptotic_Plotting',  &
               'AsymptoticTest',       &
               'AtomicWeights',        &
               'Eigenphase',           &
               'K_Gen_Test',           &
               'LgDiff',               &
               'Lifetimes',            &
               'Matrix_Out',           &
               'MxOut',                &
               'MxOut_Unit',           &
               'PoptAsy',              &
               'Print_Date_And_Time',  &
               'Print_QNumbr',         &
               'Probabilities',        &
               'Read_K1',              &
               'Read_K2',              &
               'SortEnergy',           &
               'S_Matrix',             &
               'Symmetrize',           &
               'Time_Delays',          &
               'T_Matrix',             &
               'Unit_Matrix',          &
               'Unitarity',            &
               'Write_Matrix2',        &
               'Write_Vector2',        &
               'WriteMat',             &
               'WriteVec'        /
  DATA Print_Level /'Medium',           &
                    'Medium',           &
                    'Medium',           &
                    'Medium',           &
                    'Medium',           &
                    'Medium',           &
                    'Medium',           &
                    'Medium',           &
                    'Medium',           &
                    'Medium',           &
                    'Medium',           &
                    'Medium',           &
                    'Medium',           &
                    'Medium',           &
                    'Medium',           &
                    'Medium',           &
                    'Medium',           &
                    'Medium',           &
                    'Medium',           &
                    'Medium',           &
                    'Medium',           &
                    'Medium',           &
                    'Medium',           &
                    'Medium',           &
                    'Medium',           &
                    'Medium',           &
                    'Medium',           &
                    'Medium',           &
                    'Medium',           &
                    'Medium'            /
!----------------------------------------------------------------------
!     NSubs ...... number of SUBROUTINEs which will have their debug
!                  WRITEs turned on
!
!  The following variables must be defined for each SUBROUTINE
!  for which you wish to turn on debug WRITEs:
!
!    Subs(i) ...... name of the SUBROUTINE to turn on debug WRITEs;
!                   must be a CHARACTER string of no more than
!                   eight CHARACTERs.
!    Print_Level... Default Print_level for each routine
!                   None         Nothing is printed
!                   Little       First row of matrices to be printed
!                   Medium       4-by-4 block of matrices to be printed
!                   Full         The entire matrix will be printed
!----------------------------------------------------------------------
END MODULE PoptAsy_Module
