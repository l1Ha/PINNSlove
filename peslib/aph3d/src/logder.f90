SUBROUTINE logder(Y, Yed, n, left, right, nstpl, ksqmax,method)
USE Numeric_Kinds_Module
USE Numbers_Module
USE fileunits_Module
USE Masses_Module
USE popt_Module
USE Time_Module
!            P U R P O S E   O F    S U B R O U T I N E
!  This routine propagates a set of coupled second-order differential
!  equations using Bob Johnson's Log Derivative method.

!     INPUT:
!  Y       On entering contains the Log Derivative matrix at left.
!  Yed   On entering contains the Derivative of R-Matrix at left.
!  n       Number of coupled channels.
!  left   Starting distance.
!  right   Ending distance.
!  nstpl   Number of steps per asymptotic wavelength.
!
!  LOCAL:
!  nstep   Number of grid points in interval (left, right]
!  UM      Unit matrix
!  Q       Interaction matrix Q=(2 u /hbar^2)[E-V(x)]
!  ZpI   Z+I and it also used as inverted ZpI, where Z=h*Y
!  IpQ     I+h*h/6*Q and it also used as inverted IpQ
!
!  OUTPUT:
!  Y       On return Y contains the Wiger R-matrix at right.
!  Yed   On return Yed contains the derivative of R-Matrix
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
IMPLICIT NONE
CHARACTER(LEN=10) hour, method
CHARACTER(LEN=100):: msg='none'
CHARACTER(LEN=21), PARAMETER:: ProcName='logder'
REAL(Kind=WP_Kind) timeinfo
INTEGER n, nstpl, nstep, i
REAL(Kind=WP_Kind) left, right, wavel, h, ksqmax
REAL(Kind=WP_Kind) Q(n,n), Y(n,n), Yed(n,n)
REAL(Kind=WP_Kind) UM(n,n), ZpI(n,n), IpQ(n,n)
LOGICAL deriv
LOGICAL DEBUG,firstcall
EXTERNAL smxinv, unit
DATA firstcall /.true./

IF(firstcall)THEN
    CALL get_debug(procname,debug)
    firstcall=.false.
ENDIF

IF(method /= 'logder  ' .AND. method /= 'dlogder ')THEN
    WRITE(Msg_Unit,*) 'Method = ',method
    CALL ErrorMsg(Msg_Unit,procname,'Method is not implemented!')
ENDIF
deriv= (method=='dlogder ')

IF(little.AND.Trim(msg)/='none')THEN
    IF(method=='logder ')THEN
        msg='METHOD -- Log-Deriv'
    ELSEIF(method=='logder ')THEN
        msg='METHOD -- Log-Deriv and D[Log-Deriv,E]'
    ENDIF
    WRITE(Out_Unit,*) msg
    WRITE(Msg_Unit,*) msg
ENDIF

wavel=2.d0*pi/sqrt(ksqmax)
IF(wavel>2.d0) wavel=2.d0
nstep=ceiling(nstpl*(right-left)/wavel/2.d0)*2 ! make nstep even
IF(nstep < 2) nstep=2
h=(right-left)/nstep

IF(DEBUG)THEN
   WRITE(Out_Unit,*)'ksqmax.wavel,h=',ksqmax,wavel,h
   WRITE(Out_Unit,*)'nstpl,nstep=',nstpl,nstep
ENDIF

CALL unit(n, UM)

SecndsNow =  TimeInfo('Logder','PotMat')
CALL potmat4( n, left, Q, Zero, 0)

ZpI=h*(Y-h/3*Q)+UM                              ! Compute starting (Z+I)
IF(deriv) Yed=Yed*h-h*h/3*usys2*UM              ! Compute starting Zed
   DO i = 1, nstep, 2                           ! Propagation Loop
      ! proagate an odd step ...
      CALL potmat4( n, left+h*i, Q, Zero, 0)           ! Compute inter. matrix
      IpQ=UM+h*h/6*Q                               ! IpQ= I + h^2/6*Q
      CALL smxinv (IpQ, n)                         ! invert IpQ
      CALL smxinv (ZpI, n)                         ! invert Z+Q
      IF(deriv) Yed=-MatMul(ZpI,MatMul(Yed,ZpI))+h*h*4/3*usys2*MatMul(IpQ,IpQ)
      ZpI=8*IpQ - ZpI - 6*UM                       ! Compute (Z+I)i+1
      ! proagate an even steps ...
      CALL potmat4( n, left+h*(i+1), Q, Zero, 0) ! Compute inter. matrix
      CALL smxinv (ZpI, n)                         ! invert Z+I
      IF(i+1 /= nstep)THEN                      ! if i+1=2,4,...,nstep-2
         ZpI=UM+UM -ZpI- h*h*2/3*Q                !     compute (Z+I)i+2
         IF(deriv)Yed=-Matmul(ZpI,Matmul(Yed,ZpI))-h*h*2/3*usys2*UM
      ELSE                                         ! ELSE, i.e. the last step
         Y  =( UM - ZpI - h*h/3*Q ) / h           !     compute Y
         IF(deriv) Yed=(-Matmul(ZpI,Matmul(Yed,ZpI))-h*h/3*usys2*UM ) / h
      ENDIF
   ENDDO
   
RETURN
ENDSUBROUTINE Logder
