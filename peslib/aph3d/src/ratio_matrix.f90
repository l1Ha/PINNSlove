SUBROUTINE ratio_matrix(Old2New,TNm1,TN,VecOld,VecNew,RN,NChanl,Nprop,W)
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
LOGICAL debug
CHARACTER(LEN=11) :: Blank
INTEGER i,j,info,nchanl,nprop
INTEGER ipiv(NProp)
REAL(Kind=WP_Kind) TNm1(NChanl),TN(NChanl)
REAL(Kind=WP_Kind) VecOld(NChanl,NProp)
REAL(Kind=WP_Kind) VecNew(NChanl,NProp)
REAL(Kind=WP_Kind) Old2New(NProp,NProp)
REAL(Kind=WP_Kind) RN(NProp,NProp)
REAL(Kind=WP_Kind) W(NChanl,NChanl)
REAL(Kind=WP_Kind) Work(NProp,NProp)
REAL(Kind=WP_Kind) Work2(NProp,NProp)

debug = .true.
Work = 0d0
DO i=1,NProp
    DO j=1,NProp
       Work2(i,j)=RN(i,j)*(1.d0-TNm1(j))
       Work(i,j)=Old2New(i,j)*(1.d0-TN(j))
    ENDDO
    !Work(i,i)=1.d0-TN(i)
ENDDO
IF(debug)THEN
    WRITE(Out_Unit,*)'RN matrix asis'
    CALL MatrixOut(RN,NChanl,nchanl,'RN_Matrix','RN',Blank, Blank,.False., Blank, .False.)
    WRITE(Out_Unit,*)'Work matrix asis'
    CALL MatrixOut(Work,NChanl,nchanl,'Work_Matrix','Work',Blank, Blank,.False., Blank, .False.)
ENDIF
! now W and Work are working matrices.  Calc W= Work**(-1)*W
CALL DGESV(NProp,NProp,Work,NProp,IPIV,Work2,NProp,INFO)
IF(debug)THEN
    WRITE(Out_Unit,*)'W matrix afterinv'
    CALL MatrixOut(Work2,NChanl,nprop,'W_Matrix','Work2',Blank, Blank,.False., Blank, .False.)
    IF(info>0) WRITE(Msg_Unit,*) 'info',info
ENDIF
!
! For now (00/11/08) we want the K matrix labelled by the
! primitive  basis set labels, so transform W before applying
! the boundary conditions.
!
W=MatMul(VecNew,MatMul(Work2,Transpose(VecOld)))
!
!   Transform to the diabatic basis near almost crossings
!
!CALL qbcmix(W,nchanl)
IF(debug)THEN
     WRITE(Out_Unit,*)'W matrix'
     CALL MatrixOut(W,NChanl,nchanl,'W_Matrix','W',Blank, Blank,.False., Blank, .False.)
ENDIF

RETURN
END
