SUBROUTINE trnsfm (t,w,a,nt,n1,n2,n3,n4,istop,isym)
USE Numeric_Kinds_Module
USE Popt_Module
USE FileUnits_Module
!
!-------------------------------------------------------------------
! written by g. a. parker and j. v. lill.
! this routine transforms the matrix w into a new
! basis set defined by t and stores the result
! into w again when istop=.false..  when istop=.true.
! the routine returns after the transpose of w
! is multiplied times t and the result store into a.
! when isym=.true. the resulting matrix is assumed
! symmetric.
! the matrix multiplication is parameterized as
! folllows:
!      a(j,i)=w(k,j)*t(k,i)
! where:
!        k=1,nt
!        j=n1,n2
!        i=n3,n4
! letting nd12=n2-n1+1 and nd34=n4-n3+1 the sizes
! of the arrays used are as follows: w is nt*nd12,
! t is nt*nd34, and a is thus nd12*nd34. provided
! that nd12=nt the similarity transform may be
! performed to yield a result of nd34*nd34.
! in general, w and t are assumed to share  the
! common dimension nt; the columns of w to be used
! are specified by n1 & n2, and those of t by n3 & n4.
! this determines how the arrays must be stored
! for the multiplication to be performed correctly.
!-------------------------------------------------------------------
IMPLICIT NONE
INTEGER nd12, n2, n1, nd34, n4, n3, nt, ks, ji, i, kj, j, ki, k, ij
CHARACTER(LEN=21), PARAMETER:: ProcName='trnsfm  '
LOGICAL istop, isym, cray
REAL(Kind=WP_Kind) sum
REAL(Kind=WP_Kind) t(1), w(1), a(1)
cray=.true.
full =.false.
!------------------------------------------------------------------
nd12=n2-n1+1
nd34=n4-n3+1
!-------------------------------------------------------------------
! multiply the transpose of the matrix w times t and
! store the result into matrix a.
!-------------------------------------------------------------------
IF(cray)THEN
CALL sgemm_junk(nd12,nt,nd34,w(nt*(n1-1)+1),nt,t(nt*(n3-1)+1), nt,a,nd12,4,1)
ELSE
ks=nt*(n3-2)
ji=0
DO i=1,nd34
   kj=nt*(n1-1)
   ks=ks+nt
   DO j=1,nd12
      sum=0.d0
      ki=ks
      DO k=1,nt
         kj=kj+1
         ki=ki+1
         sum=sum+w(kj)*t(ki)
      ENDDO
      ji=ji+1
      a(ji)=sum
   ENDDO
ENDDO
ENDIF
IF(full)THEN
WRITE(Out_Unit,*)'a=w(transpose)*t in trnsfm'
CALL MxOut(a,nd12,nd34)
ENDIF
IF(istop) RETURN
IF(isym) GOTO 50
!-------------------------------------------------------------------
! multiply the transpose of matrix a times matrix
! t and store the result into matrix w
!-------------------------------------------------------------------
IF(nd12/=nt) GOTO 80
IF(cray)THEN
CALL sgemm_junk(nd34,nt,nd34,a,nt,t(nt*(n3-1)+1),nt,w,nd34,4,1)
ELSE
ks=nt*(n3-2)
ji=0
DO i=1,nd34
   kj=0
   ks=ks+nt
   DO j=1,nd34
      sum=0.d0
      ki=ks
      DO k=1,nt
         kj=kj+1
         ki=ki+1
         sum=sum+a(kj)*t(ki)
      ENDDO
      ji=ji+1
      w(ji)=sum
   ENDDO
ENDDO
ENDIF
IF(full)THEN
WRITE(Out_Unit,*)'a=w(transpose)*t in trnsfm'
CALL MxOut(a,nd12,nd34)
ENDIF
RETURN
!-------------------------------------------------------------------
! this is reached only when the resulting matrix is
! symmetric hence only half the matrix need be computed
! and the other half stored by symmetry.
! multiply the transpose of matrix a times matrix
! t and store the result into matrix w
!-------------------------------------------------------------------
50   CONTINUE
IF(cray)THEN
CALL sgemm_junk(nd34,nt,nd34,a,nt,t(nt*(n3-1)+1),nt,w,nd34,4,1)
ELSE
ks=nt*(n3-2)
DO i=1,nd34
   kj=0
   ks=ks+nt
   DO j=1,i
      sum=0.d0
      ki=ks
      DO k=1,nt
         kj=kj+1
         ki=ki+1
         sum=sum+a(kj)*t(ki)
      ENDDO
      ji=j+(i-1)*nd34
      ij=i+(j-1)*nd34
      w(ji)=sum
      w(ij)=sum
   ENDDO
ENDDO
ENDIF
IF(full)THEN
WRITE(Out_Unit,*)'a=w(transpose)*t in trnsfm'
CALL MxOut(a,nd12,nd34)
ENDIF
RETURN
80 CONTINUE
WRITE(Out_Unit,90)
STOP 'trnsfm'
!----------------***end-trnsfm***---------------------------------------
!
90 FORMAT(1x,"matrix dimension error")
END
