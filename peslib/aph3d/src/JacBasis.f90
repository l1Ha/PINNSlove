SUBROUTINE JacBasis(nchanl, vecold, tstore, xsq, xk, finish, chinuj, nnuj, ntt)
!this routine is called by:
!         aph3d
!this routine calls
!         int
!-----------------------------------------------------------------------
USE Numeric_Kinds_Module
USE Narran_Module
USE FileUnits_Module
USE Integrat_Module
USE Narran_Module
USE Arrch_Module
USE region_Module
USE Parms_Module
USE Masses_Module, ONLY: AtomicSymbol
USE DiatomicPot_MODULE, ONLY: DiatomicPot
IMPLICIT NONE
INTEGER nchanl(narran), nnuj(narran), nst, ntt, ip1, ip2, ip3
REAL(Kind=WP_Kind) vecold(nvbrthrt)   !vecold(nvibrot**2)
REAL(Kind=WP_Kind) tstore(nvbrthrt)   !tstore(nvibrot**2)
REAL(Kind=WP_Kind) xsq(nvibrot)
REAL(Kind=WP_Kind) xk(nvibrot)
REAL(Kind=WP_Kind) chinuj(nvbrthrt)
REAL(Kind=WP_Kind) Finish
LOGICAL little
EXTERNAL input
little=.true.

!--------------------------------------------------------------------
! This routine calculates the asymptotic vibrational-rotational
! eigenfunctions and stores on unit=efun_unit=eigfuns.
!--------------------------------------------------------------------
iregion='jacobi '
nst=1
ntt=0

OPEN(Unit=efun_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'DiatomicOut/EigFuns.txt', form='formatted',status='unknown')
DO karran=1,narran
   nchanl(karran)=0
   nnuj(karran)=0
   IF(kArran==1)THEN
      Ip1=2
      Ip2=3
   ELSEIF(kArran==2)THEN
      Ip1=3
      Ip2=1
   ELSE
      Ip1=1
      Ip2=2
   ENDIF
   DiatomicPot=TRIM(AtomicSymbol(ip1))//TRIM(AtomicSymbol(ip2))
   IF(integrat(karran))THEN
      !--------------------------------------------------------------------
      !  CALL input which calls qlevel
      !--------------------------------------------------------------------
      Write(Out_Unit,*)'callinput',karran,nchanl(karran)
      CALL input (nchanl(karran), vecold, tstore, xsq, xk, finish, chinuj(nst), nnuj )
      nst=nst+nnuj(karran)
      !--------------------------------------------------------------------
      ! STOP execution if dimension size has been exceeded.
      !--------------------------------------------------------------------
      IF(nst>nvbrthrt)THEN
         WRITE(Out_Unit,*)'***error*** Execution terminated in APH3D'
         WRITE(Out_Unit,*) 'number of basis states', ' is greater than nvbrthrt'
         WRITE(Out_Unit,*)'nvbrthrt= ',nvbrthrt,'  nst= ',nst
         STOP 'JacBasis'
      ENDIF

      nst=nst-1
      !  ntt is the total number of vibration-rotation fcns in all channels.
      !  It includes the effect of symmetry in channel 1 but not 2 and 3.
      ntt=ntt+nchanl(karran)
      IF(little)WRITE(Out_Unit,*) 'nst=',nst,' ntt=',ntt
   ENDIF
ENDDO
!      WRITE(*,'(A,6e23.15)')"xsq2 =",xsq(1) ,xsq(36) ,xsq(71)
!      if(xsq(1)/=xsq(36).or.xsq(1)/=xsq(71))THEN
!         stop "error in routine jacbasis1"
!      endif

CLOSE(Unit=efun_unit,status='keep')
IF(ntt>nvibrot)THEN
   WRITE(Out_Unit,*) 'error: ntt is greater ', 'than nvibrot: ntt=',ntt, nvibrot
   STOP 'JacBasis'
ENDIF
RETURN
END
