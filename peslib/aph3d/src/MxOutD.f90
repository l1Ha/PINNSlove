SUBROUTINE MxOutD(a, nrow, ncol, ms)
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
INTEGER nrow, ncol, irow, jcol, krow, lcol, index, ms
REAL(Kind=WP_Kind) a(nrow,ncol), b(10)

! If ms=0, the matrix is stored as a general matrix
! If ms=1, the matrix is stored as a symmetric matrix
! If ms=2, the matrix is stored as a diagonal matrix

IF(ms==0)THEN                                  ! General matrix storage
   WRITE(Out_Unit,230)(jcol,jcol=1,min(ncol,10))
   DO irow = 1, min(nrow,30)
       WRITE(Out_Unit,240)irow, (a(irow,jcol),jcol=1,min(ncol,10))
   ENDDO
ELSEIF(ms==1)THEN                              ! Symmetric matrix storage
   WRITE(Out_Unit,230)(jcol,jcol=1,min(ncol,10))
   DO irow = 1, min(nrow,30)
      DO jcol=1,min(ncol,10)
         IF(irow>=jcol)THEN
            index=irow*(irow-1)/2+jcol
         ELSE
            index=jcol*(jcol-1)/2+irow
         ENDIF
         lcol=1+index/nrow
         krow=index-lcol*nrow
         b(jcol)=a(krow,lcol)
      ENDDO
      WRITE(Out_Unit,240)irow, (b(jcol),jcol=1,min(ncol,10))
   ENDDO
ELSEIF(ms==2)THEN                              ! Diagonal matrix storage
   WRITE(Out_Unit,230)(jcol,jcol=1,min(ncol,10))
   DO irow = 1, min(nrow,30)
      b=0.d0
      b(irow)=a(irow,1)
      WRITE(Out_Unit,240)irow, (b(jcol),jcol=1,min(ncol,10))
   ENDDO
ELSE
   WRITE(Out_Unit,*)"Incorrect value of ms: ms=", ms
   WRITE(Out_Unit,*)"Stopping in MxOutd"
   STOP "Mxoutd"
ENDIF
 
230  FORMAT(//, 3x, "  row   ", 10(3x, i3, 6x))
240  FORMAT(1x, "col ", i3, 2x, 10ES12.4)
ENDSUBROUTINE MxOutd
