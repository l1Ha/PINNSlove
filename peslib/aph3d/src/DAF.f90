SUBROUTINE DAF(sigma, delta, mval, dafx, kmax, hermit, ndaf, nd_daf)
USE Numeric_Kinds_Module
USE FileUnits_Module
USE Numbers_Module
!-----------------------------------------------------------------------
! This routine was written by G. A. Parker
! IF you find an error or have an improvement please send a messge to
! parker@ou.edu
!-----------------------------------------------------------------------
IMPLICIT NONE
INTEGER lval, mval, kmax, mval2, kval, kmaxx, kdaf, ndaf, khalf, idaf, nd_daf
REAL(KIND=WP_Kind) sigma, delta, dafx(0:nd_daf,0:ndaf), hermit(0:mval+ndaf+1), harg, expfac, fac
REAL(KIND=WP_Kind), PARAMETER :: tol=1.d-50

WRITE(Out_Unit,*)
WRITE(Out_Unit,*)' Begin(DAF)'
!-----------------------------------------------------------------------
! calculate the daf deltas.
!-----------------------------------------------------------------------
mval2=mval/2
kmaxx=kmax
DO kval=0,kmax
   harg=(kval*delta)/(SQRT(Two)*sigma)
   expfac=EXP(-harg**2)
   CALL hermit_ex(harg, mval+ndaf, hermit, expfac)
   DO kdaf=0,ndaf
      khalf=kdaf/2
      dafx(kval,kdaf)=Zero
      fac=One
      DO idaf=1,khalf
         fac=fac*idaf
      ENDDO
      DO lval=0,mval2
         IF(lval/=0)fac=(lval+khalf)*fac/lval
         dafx(kval,kdaf)=dafx(kval,kdaf)+hermit(2*lval+kdaf)*fac
      ENDDO
      dafx(kval,kdaf)=(-One/(SQRT(Two)*sigma))**kdaf*dafx(kval,kdaf)*delta/(SQRT(twopi)*sigma)*(-Four)**khalf
   ENDDO

   IF(ABS(dafx(kval,ndaf))>tol)THEN
      kmaxx=kval
   ELSEIF(ABS(dafx(kval,ndaf))<tol)THEN
      GOTO 5
   ENDIF

ENDDO

WRITE(Msg_Unit,*)'Program must be modified to increase dimension'
WRITE(Msg_Unit,*)'on array dafx'
WRITE(Msg_Unit,*)'kmaxx=',kmaxx
WRITE(Msg_Unit,*)'mval=',mval
WRITE(Msg_Unit,*)'1  nd_daf=',nd_daf
STOP 'DAF'

5  CONTINUE
IF(kmaxx>nd_daf)THEN
   WRITE(Msg_Unit,*)'Program must be modified to increase dimension'
   WRITE(Msg_Unit,*)'on array dafx'
   WRITE(Msg_Unit,*)'kmaxx=',kmaxx
   WRITE(Msg_Unit,*)'mval=',mval
   WRITE(Msg_Unit,*)'2  nd_daf=',nd_daf
   STOP 'DAF'
ENDIF

WRITE(Out_Unit,*)'   kmaxx=',kmaxx
kmax=kmaxx

WRITE(Out_Unit,*)' END(DAF)'
WRITE(Out_Unit,*)

RETURN
ENDSUBROUTINE DAF
