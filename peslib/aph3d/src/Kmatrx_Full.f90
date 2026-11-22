SUBROUTINE Kmatrx_Full(n, R_mat, RE_mat, K_mat, KE_mat, wvec, &
                chanl, elect, nvib, jrot, lorb, jtot, etot,   &
                energy, usys2, Out_Unit, temp,                &
                A_mat,  B_mat,  E_mat,  F_mat,                &
                AE_mat, BE_mat, EE_mat, FE_mat,               &
                calc_kmat, EDeriv,njacobi,gamma,phor,         &
                phoi,temp1,temp2)

!
!     Author: Gregory A. Parker,
!     Modifications/additions: Zareh Darakjian
!     Revision: 1.0
!
! revision to include photodissociation option
! m. braunstein, april, 1991.
!
!
!     On entering:
!          n           Number of coupled channels.
!          R_mat       Contains the Wigner R-matrix.
!          chanl       Arrangement channel identification.
!          elect       Electronic state.
!          nvib        Vibrational state.
!          jrot        Rotational state.
!          lorb        Orbital angular momentum.
!          etot        Total scattering energy.
!          energy      Internal energy.
!          usys2       2 time the reduced mass in atomic units.
!          gamma       photodissociation inhomogeneity
!     On exit:
!          K_mat       Contains the K-matrix.
!          phor        real part of photodissociation amplitude
!          phoi        imaginary part of photodissociation amplitude
!
USE parms_Module
USE dip_Module
USE Numbers_MODULE
USE Convrsns_Module
IMPLICIT NONE
!#include <maxhermt.parms>
!#include <nvibrot.parms>
!#include <nvbrthrt.parms>
!
!#include <Inc.dip>
!
LOGICAL little, medium, full, calc_kmat, EDeriv, testing
INTEGER n, ithcll, ithsub, i, Out_Unit, njacobi, jtot,j
INTEGER chanl(n), elect(n), nvib(n), jrot(n), lorb(n)
REAL(Kind=WP_Kind) etot, usys2
REAL(Kind=WP_Kind) R_mat(n,n),   K_mat(n,n), RE_mat(n,n), KE_mat(n,n)
REAL(Kind=WP_Kind) A_mat(n,n),   B_mat(n,n),  E_mat(n,n),  F_mat(n,n)
REAL(Kind=WP_Kind) AE_mat(n,n), BE_mat(n,n), EE_mat(n,n), FE_mat(n,n)
REAL(Kind=WP_Kind) temp(n,n),    temp1(n,n),  temp2(n,n)
REAL(Kind=WP_Kind) oopsln(nvbrthrt), wvec(n), energy(n), xksq(n)
REAL(Kind=WP_Kind) gamma(n), phor(n), phoi(n)
DATA little/.false./,medium/.false./,full/.false./, ithcll/0/,ithsub/0/

CALL popt('kmatrx  ', little, medium, full, ithcll, ithsub)
  
IF(n<=0)THEN
   WRITE(Out_Unit,*)'n=',n
   STOP 'Stopping in Kmatrx_Full'
 ENDIF
!      little=.true.
!      medium=.true.
!      full=.true.
IF(little.or.medium.or.full)THEN
   WRITE(Out_Unit,*) 'JTOT=', jtot,'  Etot=', etot
ENDIF
!
!     Incomming matrix has space-fixed labels.
!
DO i=1,n
   xksq(i)=usys2*(etot-energy(i))
   wvec(i)=sqrt(ABS(xksq(i)))
ENDDO
IF(little)THEN
   WRITE(Out_Unit,*)'Quantum Numbers used in Kmatrx_Full'
   CALL Quant_Out (n, chanl, elect, nvib, jrot, lorb, energy, xksq,'KMatrx_Full')
ENDIF
!
!     print the incoming wigner r-matrix.
!
IF(full)THEN
   WRITE(Out_Unit,*)'Wigner R_Matrix'
   CALL MxOutL(R_mat, n, n, 0, 'space', 'space')
ENDIF
!
!     match onto the asymptotic states
!
WRITE(Out_Unit,*)'Calling Delvsf_New'
CALL Delvsf_New(n, A_mat, B_mat, E_mat, F_mat, AE_mat, BE_mat, EE_mat, FE_mat, wvec, oopsln)

IF(full)THEN
   WRITE(Out_Unit,*)'R_mat'
   CALL MxOut(R_mat,n,n)
   WRITE(Out_Unit,*)'A_mat'
   CALL MxOut(A_mat,n,n)
   WRITE(Out_Unit,*)'B_mat'
   CALL MxOut(B_mat,n,n)
   WRITE(Out_Unit,*)'E_mat'
   CALL MxOut(E_mat,n,n)
   WRITE(Out_Unit,*)'F_mat'
   CALL MxOut(F_mat,n,n)
   IF(EDeriv)THEN
      WRITE(Out_Unit,*)'RE_mat'
      CALL MxOut(RE_mat,n,n)
      WRITE(Out_Unit,*)'AE_mat'
      CALL MxOut(AE_mat,n,n)
      WRITE(Out_Unit,*)'BE_mat'
      CALL MxOut(BE_mat,n,n)
      WRITE(Out_Unit,*)'EE_mat'
      CALL MxOut(EE_mat,n,n)
      WRITE(Out_Unit,*)'FE_mat'
      CALL MxOut(FE_mat,n,n)
   ENDIF
ENDIF
!
!           form -[R*F - B] --> B_mat
!
B_Mat=B_Mat-MATMUL(R_Mat,F_Mat)
Temp= B_Mat
!
!           form -[R*E - A] --> A_mat
!
A_Mat=A_Mat-MATMUL(R_Mat,E_Mat)
K_Mat=A_Mat
!
! Now Solve for the K_Matrix (see APH_theory equation 122)
!     [R F - B] K = [R E - A]
!
B_Mat=-1.d0*B_Mat
CALL lineqd (K_mat, temp, n, n)
!
! calculate real and imaginary parts of the photodissociation
! transition amplitude
!
IF(lphoto)THEN
   B_Mat=-1.d0*B_Mat
   Temp=B_Mat
   CALL smxinv(temp,n)
   CALL fill_mat(temp,n)
   CALL sgemm_junk(n,n,n,temp,n,A_mat,n,temp1,n,0,1)
   Temp2=Temp1
   CALL sgemm_junk(n,n,n,A_mat,n,temp1,n,temp,n,0,1)
   CALL specadd(n,temp,B_mat)
   CALL smxinv(temp,n)
   CALL fill_mat(temp,n)
   CALL sgemm_junk(n,n,n,temp,n,R_mat,n,temp1,n,0,1)
   CALL dgemv('N',n,n,1.d0,temp1,n,gamma,1,0.d0,phor,1)
   CALL dgemv('N',n,n,1.d0,temp2,n,gamma,1,0.d0,phoi,1)
ENDIF
IF(EDeriv)THEN
   !
   !  form RE*F + R*FE - BE --> temp
   !
   Temp=MATMUL(RE_Mat,F_Mat)+MATMUL(R_Mat,FE_Mat)-BE_Mat
   !
   !  form RE*E + R*EE - AE - (RE*F + R*FE - BE)*K --> FE_mat
   !
   FE_Mat=MATMUL(RE_Mat,E_Mat) 
   FE_Mat=FE_Mat + MATMUL(R_Mat,EE_Mat) 
   FE_Mat=FE_Mat - AE_Mat 
   FE_Mat=FE_Mat - MATMUL(Temp,K_Mat)
   !
   !  solve the linear equation for KE
   !  (R*F-B)*KE = RE*E + R*EE - AE - (RE*F + R*FE - BE)*K
   !
   !CALL vcopy (n*n, KE_mat, FE_mat)
   KE_Mat=FE_Mat
   CALL lineqd (KE_mat, B_mat, n, n)
ENDIF
!
!     print the K-matrix.
!
IF(full.AND.calc_kmat)THEN
   WRITE(Out_Unit,*)'K_Matrix at E (au) and (ev)=', etot, etot*autoev
   CALL MxOutL(K_mat, n, n, 0, 'space', 'space')
ENDIF
!
!     print the KE-matrix
!
IF(full.AND.EDeriv)THEN
   WRITE(Out_Unit,*)'KE_Matrix at E (au) and (ev)=', etot, etot*autoev
   CALL MxOutL(KE_mat, n, n, 0, 'space', 'space')
ENDIF

RETURN
END
