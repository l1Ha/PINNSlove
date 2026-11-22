SUBROUTINE read_pmat(wmat, n, nmodes, pmat_unit)
USE Numeric_Kinds_Module
USE FileUnits_Common_Module
IMPLICIT NONE
INTEGER pmat_unit, n, i, j, nmodes
REAL(Kind=WP_Kind) wmat(n,n)
!-------------------------------------------------------------------------
! wmat       Matrix to be read in from disk.
! n          Dimenison of the wmat matrix.
! nmodes     Maximum size of the wmat matrix as stored on disk.
! pmat_unit  Unit number to read from.
!-------------------------------------------------------------------------
!!!      Write(Msg_Unit,*)'read_pmat:nmodes',nmodes
IF(n>nmodes)THEN
   Write(Msg_Unit,*)'Error: n>nmodes', n,nmodes
   STOP 'read_pmat'
ENDIF

DO j=1,nmodes
   IF(j<=n)THEN
      READ(pmat_unit,err=100)(wmat(i,j),i=j,n)
      DO i = j,n
         wmat(j,i)=wmat(i,j)
      ENDDO
   ELSE
      READ(pmat_unit,err=200)
   ENDIF
ENDDO

RETURN
100   WRITE(0,*)'***Error reading wmat from the pmatrix file'
      WRITE(0,*)'     nmodes=',nmodes
      WRITE(0,*)'     n, j, i=',n, j, i
      STOP 'pmat_read'
200   STOP 'reading pmat_read'
  END
