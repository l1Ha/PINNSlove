SUBROUTINE TestXDAF(NX,kinetic,X,DafX,nd_daf)
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
INTEGER, PARAMETER :: nleg=8
INTEGER NX, IX, KX, jleg, nd_daf
REAL(Kind=WP_Kind) xval
REAL(Kind=WP_Kind) kinetic(NX,NX), X(NX)
REAL(Kind=WP_Kind) DafX(-nd_daf:nd_daf,0:2)
REAL(Kind=WP_Kind), ALLOCATABLE:: eig(:), Lmat(:,:)
REAL(Kind=WP_Kind), ALLOCATABLE:: pn(:,:), pnp(:,:), pnpp(:,:)
REAL(Kind=WP_Kind), ALLOCATABLE:: pmat(:,:), daf0(:,:), daf1(:,:), daf2(:,:)

ALLOCATE(eig(0:nleg), pn(NX,0:nleg), pnp(NX,0:nleg))
ALLOCATE(pmat(NX,0:nleg),Lmat(0:nleg,0:nleg))
ALLOCATE(daf0(NX,NX),daf1(NX,NX),daf2(NX,NX))

WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Begin(TestXDAF)'

WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'      X                    Legendre Polynomials'
DO IX=1,NX
   xval=X(IX)
   pn(IX,0)=1.d0
   pn(IX,1)=xval
   DO jleg=1,nleg-1
      pn(IX,jleg+1)=2.d0*xval*pn(IX,jleg)-pn(IX,jleg-1)-(xval*pn(IX,jleg)-pn(IX,jleg-1))/(jleg+1)
   ENDDO
   WRITE(Out_Unit,50)X(IX),(pn(IX,jleg),jleg=0,NLeg)
ENDDO

WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'      X                 D(Legendre Polynomials)'
DO IX=1,NX
   xval=X(IX)
   pnp(IX,0)=0.d0
   DO jleg=1,nleg
      pnp(IX,jleg)=-jleg*(xval*pn(IX,jleg)-pn(IX,jleg-1))/(1-xval**2)
   ENDDO
   WRITE(Out_Unit,50)X(IX),(pnp(IX,jleg),jleg=0,NLeg)
ENDDO

Lmat=0.d0
DO jleg=0,nleg
   Lmat(jleg,jleg)=jleg*(jleg+1)
ENDDO


  50  FORMAT(1x,11es12.4)

daf0=0.d0
daf1=0.d0
daf2=0.d0
DO IX=1,Min(NX,nd_daf)
   DO KX=1,Min(NX,nd_daf)
      daf0(IX,KX)=DafX(IX-KX,0)
      daf1(IX,KX)=DafX(IX-KX,1)
      daf2(IX,KX)=DafX(IX-KX,2)
   ENDDO
ENDDO

WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Daf(0)'
CALL MxOut(Daf0,NX,NX)

WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Daf(1)'
CALL MxOut(Daf1,NX,NX)

WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Daf(2)'
CALL MxOut(Daf2,NX,NX)

WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Kinetic'
CALL MxOut(Kinetic,NX,NX)

      !Test Pn-Daf(0)*Pn
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Test Pn-Daf(0)*Pn'
Pmat=Pn-Matmul(Daf0,Pn)
CALL MxOut(Pmat,NX,NLeg+1)

      !Test Pnp-Daf(1)*Pn
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Test Pnp-Daf(1)*Pn'
Pmat=Pnp-Matmul(Daf1,Pn)
CALL MxOut(Pmat,NX,NLeg+1)

      !Test L*(L+1)*Pn-Daf(1)*[1-x^2]*Daf(1)*Pn
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Test L*(L+1)*Pn-Daf(1)*[1-x^2]*Daf(1)*Pn'
Pmat=MatMul(Pn,Lmat)-Matmul(Daf1,Pn)
CALL MxOut(Pmat,NX,NLeg+1)

      !Test Number L*(L+1)*Pn-Kinetic*Pn
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Test L*(L+1)*Pn-Kinetic*Pn'
Pmat=Matmul(Pn,Lmat)-Matmul(Kinetic,Pn)
CALL MxOut(Pmat,NX,NLeg+1)

WRITE(Out_Unit,*)'END(TestXDAF)'
WRITE(Out_Unit,*)

RETURN

ENDSUBROUTINE TestXDAF
