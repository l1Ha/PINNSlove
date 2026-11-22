SUBROUTINE manolop(Y, n, left, right, nstpl, ksqmax, method)
USE Numeric_Kinds_Module
USE Numbers_Module
USE fileunits_Module
USE Masses_Module
USE popt_Module
USE Time_Module

!            P U R P O S E   O F    S U B R O U T I N E
!  This routine propagates a set of coupled second-order differential
!  equations using Manolopulos's Improved Log Derivative Method.
!
!  References: D.E. Manolopulos, J. Chem. Phys. 85, 6425 (1986)
!           D.E. Manolopulos, J. Chem. Phys. 105, 169 (1993)   
!
!     INPUT:
!  Y       On entering contains the Log Derivative matrix at left.
!  n       Number of coupled channels.
!  left   Starting distance.
!  right   Ending distance.
!  nstpl   Number of steps per asymptotic wavelength.
!
!  LOCAL:         
!  p2     p^2     
!  p       |sqrt(|p2|)|
!  p2, y1, y2, Qa, Qb, Qc   see the above references
!  UM      Unit matrix        
!  T       temporary matrix
!
!  Grids Scheme:  
!
!     left       h    h                     right  
!     |-----|-----|--+--|--------------------------|
!         0     1  ...                                nstep         
!     a  c  b     a  c  b     a  c  b   ...   a  c  b
!          a  c  b     a  c  b     a  ...    
!  
!  OUTPUT:
!  Y       On return Y contains the Wiger R-matrix at right.
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
IMPLICIT NONE
CHARACTER(LEN=10) hour
CHARACTER(LEN=10) method
CHARACTER(LEN=100) msg
CHARACTER(LEN=21), PARAMETER:: ProcName='manolop '
INTEGER n, nstpl, nstep, i, j
REAL(Kind=WP_Kind) left, right, wavel, h, ksqmax
REAL(Kind=WP_Kind) T(n,n), Wref(n,n), Qa(n,n),Qb(n,n),Qc(n,n),Y(n,n)
REAL(Kind=WP_Kind) p2(n),p(n) 
REAL(Kind=WP_Kind) UM(n,n), y1(n,n), y2(n,n), Q(n,n)
REAL(Kind=WP_Kind) a, b, c, TimeInfo
EXTERNAL smxinv, unit

CALL popt(procname, little, medium, full, ithcll, ithsub)

IF(little)THEN 
msg='manolop --Improved Log-Deriv. Propaga. with Invers.'
WRITE(Out_Unit,*) msg
WRITE(msg_unit,*) msg
ENDIF

t1 = TimeInfo('Manolop','Write')
WRITE(msg_unit,*)' In manolop, calling potmat: ',hour

wavel=2.d0*pi/sqrt(ksqmax)
IF(wavel>2.d0) wavel=2.d0    
nstep=nstpl*(right-left)/wavel   
IF(nstep < 1) nstep=1
h=(right-left)/nstep/2

WRITE(Out_Unit,*)'ksqmax.wavel,h=',ksqmax,wavel,h
WRITE(Out_Unit,*)'nstpl,nstep=',nstpl,nstep

CALL unit(n, UM)
DO i = 1, nstep                     !Propagation Loop
a= left+(i-1)*h*2
b= left+  i  *h*2
c= (a+b)/2

CALL potmat4( n, c, Q, Zero, 0)        
Wref=0d0
Do j=1, n
   p2(j)=-Q(j,j)              !P^2(j)=W(c)jj, W=-Q       
   p(j)=SQRT(ABS(p2(j)))
   Wref(j,j)=p2(j)   
ENDDO
T=UM + h*h/6*(Q+Wref)   
Call SmxInv(T, n)
Qc=4/h*(T-UM)                 !Q(c)=4/h*{[I+h*h/6*U(c)]^(-1)-I}

CALL potmat4( n, a, Q, Zero, 0)        
Qa=-h/3*(Q+Wref)              !Q(a)=h/3*U(a)          

CALL potmat4( n, b, Q, Zero, 0)        
Qb=-h/3*(Q+Wref)              !Q(b)=h/3*U(b)          

y1=0d0
y2=0d0
Do j=1, n
   IF(p2(j)>=0)THEN
      y1(j,j)=p(j)/tanh(p(j)*h) !coth=1/tanh
      y2(j,j)=p(j)/sinh(p(j)*h) !csch=1/sinh
   ELSE
      y1(j,j)=p(j)*1d0/tan(p(j)*h) 
      y2(j,j)=p(j)/sin(p(j)*h)  !csc=1/sin
   ENDIF
ENDDO

                           !(r', r") = (a, c)
T=Y+y1+Qa                     !Y + y1 + Q(r')
CALL SmxInv(T, n)             !invert [Y + y1 + Q(r')]
Y=y1+Qc-Matmul(y2, Matmul(T, y2) )  !propagate Y from a to c 

                           ! (r', r") = (c, b)
T=Y+y1+Qc                     !Y + y1 + Q(r')
CALL SmxInv(T, n)             !invert [Y + y1 + Q(r')]
Y=y1+Qb-Matmul(y2, Matmul(T, y2) )  !propagate Y from c to b 

ENDDO               

END
