SUBROUTINE Kmatrx_Open(nt, nopen, K_mat, KE_mat, wvec, nbasis,      &
                      kchan, elect, mvib, jrot, lorb, xksq,         &
                      energy, nchanl, narran, Out_Unit, calc_kmat,  &
                      EDeriv, phor, phoi)
!-----------------------------------------------------------------------
!
!     Author: Gregory A. Parker
!
!     This routine takes the K-matrix with both open and closed channel
!     contributions and deletes all of the open-closed, closed-open and
!     closed-closed elements. Hence, on RETURN from this routine K_mat
!     is the K-matrix with only open channels.
!
!  On Entering:
!     nt       Number of coupled channels (both open and closed).
!     K_mat    K_matrix with open and closed channels.
!     KE_mat   derivative of K_matrix with respect to energy.
!     kchan    Arrangement channel labels.
!     elect    Electronic state.
!     mvib     Vibrational state.
!     jrot     Rotational state.
!     lorb     Orbital angular momentum state.
!     xksq     Square of the wavenumber.
!     wvec     Wavenumber.
!     narran   Number of reactive arrangement channels.
!
!     phor     real part of photodissociation amplitude
!     phoi     imaginary part of photodissociation amplitude
!
!
!  On Exit:
!     nopen    Number of open channels.
!     K_mat    K_matrix with open channels only.
!     KE_mat   derivative of K_matrix with respect to energy.
!     nchanl   Number of states in each reactive arrangement channel.
!     kchan    Arrangement channel labels.
!     elect    Electronic state.
!     mvib     Vibrational state.
!     jrot     Rotational state.
!     lorb     Orbital angular momentum state.
!     xksq     Square of the wavenumber.
!     wvec     Wavenumber.
!
!-----------------------------------------------------------------------
USE dip_Module
USE NState_Module
IMPLICIT NONE
INTEGER nt, nopen, ithcll, ithsub, in, ij, j, i, iarran, narran, Out_Unit
INTEGER nbasis(nt), kchan(nt), elect(nt), mvib(nt), jrot(nt), lorb(nt), nchanl(narran)
LOGICAL little, medium, full, EDeriv, calc_kmat
REAL(Kind=WP_Kind) difi, difj, err
REAL(Kind=WP_Kind) K_mat(nt*nt), KE_mat(nt*nt), wvec(nt), xksq(nt), energy(nt)
REAL(Kind=WP_Kind) phor(nt), phoi(nt)
DATA ithcll/0/, ithsub/0/, little/.false./, medium/.false./, full/.false./

CALL popt('Kmatrx_Open', little, medium, full, ithcll, ithsub)

Little=.False.
IF(little)THEN
   WRITE(Out_Unit,*)'Quantum Numbers used in Kmatrx_Open'
   CALL Quant_Out(nt, kchan, elect, mvib, jrot, lorb, energy, xksq,'KMatrx_Open_1')
ELSEIF(medium)THEN
   WRITE(Out_Unit,*)'Quantum Numbers used in Kmatrx_Open'
   CALL Quant_Out(nt, kchan, elect, mvib, jrot, lorb, energy, xksq,'KMatrx_Open_2')
   IF(calc_kmat)THEN
      WRITE(Out_Unit,*)'K_matrix with open and closed channels'
      CALL MxOut(K_mat, nt, nt)
   ENDIF
   IF(EDeriv)THEN
      WRITE(Out_Unit,*)'Derivative of the K_matrix with respect',  &
                       'to Energy with both open and closed channels.'
      CALL MxOut(KE_mat, nt, nt)
   ENDIF
ELSEIF(full)THEN
   WRITE(Out_Unit,*)'Quantum Numbers used in Kmatrx_Open'
   CALL Quant_Out(nt, kchan, elect, mvib, jrot, lorb, energy, xksq,'KMatrx_Open_3')
   IF(calc_kmat)THEN
      WRITE(Out_Unit,*)'K_matrix with open and closed channels'
      CALL MxOutD(K_mat, nt, nt, 0)
   ENDIF
   IF(EDeriv)THEN
      WRITE(Out_Unit,*)'Derivative of the K_matrix with respect',  &
                      'to Energy with both open and closed channels.'
      CALL MxOutD(KE_mat, nt, nt, 0)
   ENDIF
ENDIF
!---------------------------------------------------------------------
! Keep the open-open portion of the K_matrix and delete the remaining
! portion.
!---------------------------------------------------------------------
in=0
nopen=0
Nopen_Max=0
ij=0
WRITE(Out_Unit,*)'Kmatrx_Open',nt
DO j=1,nt
   difj=xksq(j)
   wvec(j)=sqrt(ABS(difj))
   IF(difj > 0.d0)THEN
      nopen=nopen+1
      nbasis(nopen)=j
   ENDIF
   DO i=1,nt
      ij=ij+1
      difi=xksq(i)
      IF(difj>0.d0.AND.difi>0.d0)THEN
         in=in+1
         IF(lphoto)THEN
            IF(j==1)THEN
               phor(in)=phor(ij)
               phoi(in)=phoi(ij)
            ENDIF
         ENDIF
         K_mat(in)  = K_mat(ij)
         KE_mat(in) = KE_mat(ij)
      ENDIF
   ENDDO
ENDDO
IF(Nopen>Nopen_Max)Nopen_Max=Nopen

!--------------------------------------------------------------------
! Make sure that there are some open channels at this energy.
!--------------------------------------------------------------------
IF(nopen <= 0)THEN
   WRITE(*,*)'ERROR: nopen must be > 0; Energy is too small'
   WRITE(Out_Unit,*)'ERROR: nopen must be > 0; Energy is too small'
   STOP 'Stopping in Kmatrx_Open'
ENDIF
!--------------------------------------------------------------------
! Calculate the % of asymmetry.
!--------------------------------------------------------------------
CALL asym(K_mat,nopen,err)
!---------------------------------------------------------------------
! Reorder channel labels.
!---------------------------------------------------------------------
DO i=1,nopen
   mvib(i)=mvib(nbasis(i))
   jrot(i)=jrot(nbasis(i))
   lorb(i)=lorb(nbasis(i))
   kchan(i)=kchan(nbasis(i))
   xksq(i)=xksq(nbasis(i))
   wvec(i)=wvec(nbasis(i))
   elect(i)=elect(nbasis(i))
   energy(i)=energy(nbasis(i))
ENDDO

!----------------------------------------------------------------------
! Count the number of channels in each arrangement
!----------------------------------------------------------------------
DO iarran = 1, narran
    nchanl(iarran) = 0
ENDDO
DO i=1,nopen
   DO iarran=1,narran
      IF(kchan(i)==iarran)nchanl(iarran)=nchanl(iarran)+1
   ENDDO
ENDDO

WRITE(Out_Unit,*)'Number of Open Channels: nopen=',nopen
WRITE(Out_Unit,*)'Quantum Numbers used in Kmatrx_Open Open Channels Only'
CALL Quant_Out(nopen, kchan, elect, mvib, jrot, lorb, energy, xksq,'KMatrx_Open_4')
WRITE(Out_Unit,*)' Asymmetry of the open channels part of',   &
                 ' the K matrix: err=',err
IF(medium)THEN
   IF(calc_kmat)THEN
      WRITE(Out_Unit,*)'K_matrix with only open channels'
      CALL MxOut(K_mat, nopen, nopen)
   ENDIF
   IF(EDeriv)THEN
      WRITE(Out_Unit,*)'Derivative of the K_matrix with respect',   &
                       'to Energy with only open channels.'
      CALL MxOut(KE_mat, nopen, nopen)
   ENDIF
ELSEIF(full)THEN
   IF(calc_kmat)THEN
      WRITE(Out_Unit,*)'K_matrix with only open channels'
      CALL MxOutD(K_mat, nopen, nopen, 0)
   ENDIF
   IF(EDeriv)THEN
      WRITE(Out_Unit,*)'Derivative of the K_matrix with respect',   &
                       'to Energy with only open channels.'
      CALL MxOutD(KE_mat, nopen, nopen, 0)
   ENDIF
ENDIF

RETURN
ENDSUBROUTINE Kmatrx_Open
