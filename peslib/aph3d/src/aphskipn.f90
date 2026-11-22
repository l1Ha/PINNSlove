SUBROUTINE aphskipn(nrho, rhoval, lam, unit, id)
USE Numeric_Kinds_Module
USE CSBasis_Module
IMPLICIT NONE

!#include <Inc.csbasis>

LOGICAL mid_zero
INTEGER nrho, nmodes, i, j, lam, id, mega, megap, unit, nmodesp
REAL(Kind=WP_Kind) rhoval(nrho), rhocnt, rhofst

!-----------------------------------------------------------------------
! id = 0   skip pmatrix
!    = 1   skip diag_jtot, cor_unit
!-----------------------------------------------------------------------

IF(id==0)THEN
!        READ(pmat_unit)rhocnt,nmodes,nrho,(rhoval(i),i=1,nrho), mid_zero, r_cstst
!        READ(pmat_unit)(energy(i),i=1,n)
!        READ(pmat_unit,err=100)(wmat(i,j),i=j,n)

   READ(unit)rhocnt,nmodes,nrho,(rhoval(i),i=1,nrho),mid_zero
   READ(unit)

   DO j=1,nmodes
      READ(unit)
   ENDDO

   IF(.NOT.mid_zero)THEN
      IF(lam==1.AND..NOT.cstest)THEN
         DO j=1,nmodes
            READ(unit)
         ENDDO
      ENDIF
   ENDIF

   DO j=1,nmodes
      READ(unit)
   ENDDO
!-----------------------------------------------------------------------
! READ(jtot_unit)rhocent,rhofst,mega,megap,nmodes,nmodesp,ident
! READ(jtot_unit,err=100)(diag(i,j),i=j,n)
! or READ(cor_unit,err=100)(cor(i,j),i=1,n1)

ELSEIF(id==1)THEN
   READ(unit)rhocnt,rhofst,mega,megap,nmodes, nmodesp
   DO j=1,nmodesp
      READ(unit)
   ENDDO
ENDIF

RETURN
ENDSUBROUTINE aphskipn
