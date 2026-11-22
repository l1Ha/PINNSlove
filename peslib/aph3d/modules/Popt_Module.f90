MODULE Popt_Module
   USE Numeric_Kinds_Module
   IMPLICIT NONE
   SAVE
!     NSubs ...... number of SUBROUTINEs which will have their debug
!                  WRITEs turned on
   INTEGER, PARAMETER:: NSubs=12
!    Subs(i) ...... name of the SUBROUTINE to turn on debug WRITEs;
!                   must be a CHARACTER string of no more than
!                   eight CHARACTERs.
   CHARACTER(LEN=21)  Subs(NSubs) &
    /'K_Gen_Test',     &
     'Unit_Matrix',    &
     'Asymptotic',     &
     'Eigenphase',     &
     'LgDiff',         &
     'Lifetimes',      &
     'Probabilities',  &
     'S_Matrix',        &
     'Symmetrize',     &
     'Time_Delays',    &
     'T_Matrix',        &
     'Unitarity'       /
!    Print_Level... Default Print_level for each routine
!                   NONE         Nothing is printed
!                   Little       First row of matrices to be printed
!                   Medium       4-by-4 block of matrices to be printed
!                   Full         The entire matrix will be printed
   CHARACTER(LEN=6) Print_Level(NSubs) &
    /'Medium',           &
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
   LOGICAL:: little=.False.        !
   LOGICAL:: medium=.FALSE.       !
   LOGICAL:: full=.FALSE.         !
   INTEGER(KIND=IW_Kind):: ithcll=0  !
   INTEGER(KIND=IW_Kind):: ithsub=0  !
   LOGICAL AllNone/.False./, AllLittle/.False./, AllMedium/.False./, AllFull/.False./
END MODULE Popt_Module
