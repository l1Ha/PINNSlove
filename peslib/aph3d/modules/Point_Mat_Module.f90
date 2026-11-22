MODULE Point_Mat_Module
   USE Numeric_Kinds_Module
   IMPLICIT NONE
   SAVE
   REAL(Kind=WP_Kind), ALLOCATABLE:: fbar(:, :, :)  ! fbar(nq, nquad, nquad)
   REAL(Kind=WP_Kind), ALLOCATABLE:: dfbar(:, :, :) ! dfbar(nq, nquad, nquad)
   REAL(Kind=WP_Kind), ALLOCATABLE:: hamil(:, :)    ! hamil(nq, nq)
   REAL(Kind=WP_Kind), ALLOCATABLE:: ff(:, :)       ! ff(numnp, nq)
END MODULE Point_Mat_Module
