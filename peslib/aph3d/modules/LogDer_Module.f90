MODULE LogDer_Module
   USE Numeric_Kinds_Module
   USE Ederiv_Module
   USE Boundary_Module
   USE Regins_Module
   USE DBS_Module
   IMPLICIT NONE
   SAVE
   CHARACTER(LEN=6)   :: iprtx      ='none  '        !
   CHARACTER(LEN=10)  :: method     ='new_logder'    !
   CHARACTER(LEN=12)  :: init_rmat  ='Initialize  '  !
   CHARACTER(LEN=100) :: rmat_file  ='R_matrix'      !
   NAMELIST/logder/ nsteps, iprtx, init_rmat, rmat_file, start, finish, method, Ederiv
   ! nsteps: Number of steps per asymptotic wavelength used in the Log-derivative method.
END MODULE LogDer_Module
