SUBROUTINE kin_theta (ntheta, kinetic, theta, kmax, usys2, dafx, nd_daf, ndaf)
USE FileUnits_Module
USE DAFSigma_Module
USE InputFile_Module
IMPLICIT NONE
LOGICAL :: debug
INTEGER :: nd_daf, ndaf, k, dafopt, NTheta2
INTEGER :: i, j, ntheta, m, kmax, mloop, mp, lwork, info
REAL(Kind=WP_Kind) :: kinetic (ntheta, ntheta), usys2, dafx (0:nd_daf, 0: &
 ndaf), theta (ntheta), ctsq, xdaf, ddt, d1, d2, d3, d4, ww, &
 dtheta, d1daf, dik, dkj, eig_im (10000), norm
REAL(Kind=WP_Kind), ALLOCATABLE::daft (:, :), dafp (:, :), psi (:, :)
REAL(Kind=WP_Kind), ALLOCATABLE::eig (:)
NAMELIST/dafthetaopt/ dafopt
WRITE(Out_Unit, * )

WRITE(Out_Unit, * ) 'Begin(Kin_Theta)'
ALLOCATE (daft ( - nd_daf:nd_daf, 0:2), psi (ntheta, ntheta) )
ALLOCATE (dafp ( - nd_daf:nd_daf, 0:2), eig (ntheta) )
!-temp
OPEN(Unit=xdaf_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'GraphicsOut/xdafplt')
dtheta = theta (2) - theta (1)
DO i = 1, ntheta
   xdaf = d1daf (theta (1), theta (i), 16, sigma)
   WRITE(xdaf_Unit, * ) theta (i), dafx (i - 1, 1), xdaf * dtheta
ENDDO
dafopt = 314
OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile)
READ(In_Unit,NML=dafthetaopt)
WRITE(Out_Unit, * ) 'dafopt=', dafopt
CLOSE(unit=In_Unit)
WRITE(Out_Unit,NML=dafthetaopt)

debug = .true.
IF(ntheta<=0)THEN
   kinetic (1, 1) = 0.d0
   RETURN
ENDIF

kinetic = 0.d0
WRITE(Out_Unit, * ) 'Calculating Theta_Kinetic'
mloop = kmax / ntheta + 1

WRITE(Out_Unit, * ) 'mloop=', mloop

daft = 0.d0
DO i = 0, kmax
   daft (i, 0) = dafx (i, 0)
   daft (i, 1) = dafx (i, 1)
   daft (i, 2) = dafx (i, 2)
   daft ( - i, 0) = daft (i, 0)
   daft ( - i, 1) = - daft (i, 1)
   daft ( - i, 2) = daft (i, 2)
ENDDO
DO i = - nd_daf, nd_daf
   DO j = 0, 2
      dafp (i, j) = 0.d0
      DO m = - mloop, mloop
         k = i + 2 * m * (ntheta - 1)
         IF(iabs (k) <=kmax)THEN
            dafp (i, j) = dafp (i, j) + daft (k, j)
         ENDIF
      ENDDO
   ENDDO
ENDDO
!----
! NOTE: This produces a non-symmetric matrix
!----
IF(dafopt==1)THEN
   DO i = 1, ntheta
      DO j = 1, ntheta
         kinetic (i, j) = sin (2.d0 * theta (i) ) * (dafp (i - j, 2) &
          + dafp (i - (2 * ntheta - j), 2) ) + 2.d0 * cos (2.d0 * theta ( &
          i) ) * (dafp (i - j, 1) + dafp (i - (2 * ntheta - j), 1) )
      ENDDO
   ENDDO
   DO k = 1, ntheta
      kinetic (k, 1) = kinetic (k, 1) / sqrt (2.d0)
      kinetic (k, ntheta) = kinetic (k, ntheta) / sqrt (2.d0)
   ENDDO
ENDIF
!----
! NOTE: This produces a non-symmetric matrix
!----
IF(dafopt==411)THEN
   DO i = 1, ntheta
      IF(i==1.or.i==ntheta)THEN
         DO j = 1, ntheta
         kinetic (i, j) = 2.d0 * (dafp (i - j, 2) + dafp (i - (2 * ntheta - j), 2) )
         ENDDO
      ELSE
         DO j = 1, ntheta
         kinetic (i, j) = (dafp (i - j, 2) + dafp (i - (2 * ntheta - &
          j), 2) ) + 2.d0 * cos (2.d0 * theta (i) ) / sin (2.d0 * &
          theta (i) ) * (dafp (i - j, 1) + dafp (i - (2 * ntheta - j),1) )
         ENDDO
      ENDIF
   ENDDO
   DO i = 1, ntheta
   kinetic (i, 1) = kinetic (i, 1) / 2.d0
   kinetic (i, ntheta) = kinetic (i, ntheta) / 2.d0
   ENDDO
   IF(ntheta== - 1)THEN
!         CALL dgeev('N','V',ntheta,ksave,ntheta,eig,eig_im,0.d0,ntheta,psi,ntheta,work,lwork,info)
      DO i = 1, ntheta
      WRITE(Out_Unit, * ) 'eig_im=', eig_im (i) , - 4.d0 * eig (i) / usys2
      ENDDO
      CALL sorter (eig, eig_im, ntheta, psi)
      DO i = 1, ntheta
         DO j = 1, ntheta
            kinetic (i, j) = 0.d0
         ENDDO
      ENDDO
      DO k = 1, ntheta
         norm = 0.d0
         DO i = 1, ntheta
            norm = psi (i, k) * sin (2.d0 * theta (i) ) * psi (i, k) + norm
         ENDDO
         DO i = 1, ntheta
            psi (i, k) = sqrt (sin (2.d0 * theta (i) ) ) * psi (i, k) / sqrt (norm)
         ENDDO
      ENDDO
      DO i = 1, ntheta
         DO j = 1, ntheta
            DO k = 2, 2
               kinetic (i, j) = psi (i, k) * eig (k) * psi (j, k) +  kinetic (i, j)
            ENDDO
         ENDDO
      ENDDO
   ENDIF


ENDIF
!----
! NOTE:
!----
IF(dafopt==11)THEN
   DO i = 1, ntheta
      DO j = 1, ntheta
         kinetic (i, j) = 0.d0
         DO k = 1, ntheta
            dik = sqrt (sin (2.d0 * theta (k) ) ) * (dafp (k - i, 1)  + dafp (k - (2 * ntheta - i), 1) )
            dkj = sqrt (sin (2.d0 * theta (k) ) ) * (dafp (k - j, 1)  + dafp (k - (2 * ntheta - j), 1) )
            kinetic (i, j) = kinetic (i, j) - dik * dkj
         ENDDO
      ENDDO
   ENDDO
!         DO k=1,ntheta
!             kinetic(1,k)=kinetic(1,k)/sqrt(2.d0)
!             kinetic(k,1)=kinetic(k,1)/sqrt(2.d0)
!             kinetic(ntheta,k)=kinetic(ntheta,k)/sqrt(2.d0)
!             kinetic(k,ntheta)=kinetic(k,ntheta)/sqrt(2.d0)
!         ENDDO

ENDIF
!----
! NOTE:
!----
IF(dafopt==111)THEN
   ntheta2 = 2 * ntheta - 1
   OPEN(Unit=Ntheta2_Data_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//'ntheta2.data')
   READ(Ntheta2_Data_Unit, * ) ntheta2
   CLOSE(unit=Ntheta2_Data_Unit)
   ddt = 2.d0 * atan (1.d0) / (ntheta2 - 1)
   Write(Msg_Unit, * ) 'ddt=', ddt
   DO i = 1, ntheta
   DO j = 1, ntheta
   kinetic (i, j) = 0.d0
   DO k = 1, ntheta2
   d1 = 0.d0
   d2 = 0.d0
   d3 = 0.d0
   d4 = 0.d0
   DO m = - mloop - 1, mloop + 1
!c                  d1=d1+daft(k-i+2*m*(ntheta-1),1)
   d1 = d1 + d1daf ( (k - 1) * ddt, (i - 1 + 2 * m * (ntheta - 1) &
    ) * dtheta, 16, sigma) * ddt
!c                  d2=d2+daft(k-(2*ntheta-i)+2*m*(ntheta-1),1)
   d2 = d2 + d1daf ( (k - 1) * ddt, (2 * ntheta - i - 1 + 2 * m * &
    (ntheta - 1) ) * dtheta, 16, sigma) * ddt
!c                  d3=d3+daft(k-j+2*m*(ntheta-1),1)
   d3 = d3 + d1daf ( (k - 1) * ddt, (j - 1 + 2 * m * (ntheta - 1) &
    ) * dtheta, 16, sigma) * dtheta
!c                  d4=d4+daft(k-(2*ntheta-j)+2*m*(ntheta-1),1)
   d4 = d4 + d1daf ( (k - 1) * ddt, (2 * ntheta - j - 1 + 2 * m * &
    (ntheta - 1) ) * dtheta, 16, sigma) * dtheta
   ENDDO
   dik = sqrt (sin (2.d0 * (k - 1) * ddt) ) * (d1 + d2)
   dkj = sqrt (sin (2.d0 * (k - 1) * ddt) ) * (d3 + d4)
   kinetic (i, j) = kinetic (i, j) - dik * dkj
   ENDDO
   ENDDO
   ENDDO
   DO k = 1, ntheta
   kinetic (1, k) = kinetic (1, k) / sqrt (2.d0)
   kinetic (k, 1) = kinetic (k, 1) / sqrt (2.d0)
   kinetic (ntheta, k) = kinetic (ntheta, k) / sqrt (2.d0)
   kinetic (k, ntheta) = kinetic (k, ntheta) / sqrt (2.d0)
   ENDDO


ENDIF
!----
! NOTE:
!----
IF(dafopt==214)THEN
   ntheta2 = 2 * ntheta - 1
   OPEN(Unit=Ntheta2_Data_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//'ntheta2.data')
   READ(Ntheta2_Data_Unit, * ) ntheta2
   CLOSE(unit=Ntheta2_Data_Unit)
   ddt = 2.d0 * atan (1.d0) / (ntheta2 - 1)
   Write(Msg_Unit, * ) 'ddt=', ddt
   DO i = 1, ntheta
   DO j = 1, ntheta
   kinetic (i, j) = 0.d0
   DO m = - 1, 1
   DO mp = - 1, 1
   DO k = 1, ntheta2
   d1 = 0.d0
   d2 = 0.d0
   d3 = 0.d0
   d4 = 0.d0
   d1 = d1 + d1daf ( (i - 1 + 2 * mp * (ntheta - 1) ) * dtheta, &
    (k - 1) * ddt, 16, sigma) * ddt
   d2 = d2 + d1daf ( (2 * ntheta - i - 1 + 2 * mp * (ntheta - 1) ) &
    * dtheta, (k - 1) * ddt, 16, sigma) * ddt
   d3 = d3 + d1daf ( (k - 1) * ddt, (j - 1 + 2 * m * (ntheta - 1) &
    ) * dtheta, 16, sigma) * dtheta
   d4 = d4 + d1daf ( (k - 1) * ddt, (2 * ntheta - j - 1 + 2 * m * &
    (ntheta - 1) ) * dtheta, 16, sigma) * dtheta
   dik = sqrt (sin (2.d0 * (k - 1) * ddt) ) * (d1 + d2)
   dkj = sqrt (sin (2.d0 * (k - 1) * ddt) ) * (d3 + d4)
   kinetic (i, j) = kinetic (i, j) + dik * dkj
   ENDDO
   ENDDO
   ENDDO
   ENDDO
   ENDDO
   DO k = 1, ntheta
   kinetic (1, k) = kinetic (1, k) / sqrt (2.d0)
   kinetic (k, 1) = kinetic (k, 1) / sqrt (2.d0)
   kinetic (ntheta, k) = kinetic (ntheta, k) / sqrt (2.d0)
   kinetic (k, ntheta) = kinetic (k, ntheta) / sqrt (2.d0)
   ENDDO

ENDIF
!----
! NOTE:
!----
IF(dafopt==314)THEN
   DO i = 1, ntheta
      DO j = 1, ntheta
         kinetic (i, j) = 0.d0
         DO m = - 1, 1
            DO mp = - 1, 1
               DO k = 1, ntheta
                  IF(k==1)THEN
                     ww = 0.5d0
                  ELSEIF (k==ntheta)THEN
                     ww = 0.5d0
                  ELSE
                     ww = 1.d0
                  ENDIF
                  d1 = 0.d0
                  d2 = 0.d0
                  d3 = 0.d0
                  d4 = 0.d0
                  d1 = d1 + daft ( (i - 1 + 2 * mp * (ntheta - 1) ) - (k - 1),  1)
                  d2 = d2 + daft ( (2 * ntheta - i - 1 + 2 * mp * (ntheta - 1) ) - (k - 1), 1)
                  d3 = d3 + daft ( (k - 1) - (j - 1 + 2 * m * (ntheta - 1) ), 1)
                  d4 = d4 + daft ( (k - 1) - (2 * ntheta - j - 1 + 2 * m * (ntheta - 1) ), 1)
                  dik = sqrt (sin (2.d0 * (k - 1) * dtheta) ) * (d1 + d2)
                  dkj = sqrt (sin (2.d0 * (k - 1) * dtheta) ) * (d3 + d4)
                  kinetic (i, j) = kinetic (i, j) + dik * dkj * ww
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   DO k = 1, ntheta
      kinetic (1, k) = kinetic (1, k) / sqrt (2.d0)
      kinetic (k, 1) = kinetic (k, 1) / sqrt (2.d0)
      kinetic (ntheta, k) = kinetic (ntheta, k) / sqrt (2.d0)
      kinetic (k, ntheta) = kinetic (k, ntheta) / sqrt (2.d0)
   ENDDO


ENDIF
!----
! NOTE: There are two extra eigenvalues doing it this way
!----
IF(dafopt==12)THEN
   DO i = 1, ntheta
   DO j = 1, ntheta
   kinetic (i, j) = dafp (i - j, 2) - dafp (i + j - 2, 2)
   ENDDO
   ENDDO
   DO i = 1, ntheta
   ctsq = 0.d0
   IF(i/=1.AND.i/=ntheta)THEN
      ctsq = (cos (2.d0 * theta (i) ) / sin (2.d0 * theta (i) ) ) **2
   ELSE
      ctsq = - 2.d+0
   ENDIF
   kinetic (i, i) = kinetic (i, i) + 2.d0 + ctsq
   ENDDO
   DO i = 1, ntheta
   DO j = 1, ntheta
   kinetic (i, j) = sqrt (sin (2.d0 * theta (i) ) ) * kinetic (i, j) * sqrt (sin (2.d0 * theta (j) ) )
   ENDDO
   ENDDO


ENDIF
!----
! NOTE: There are two extra eigenvalues doing it this way
!----

IF(dafopt==2)THEN
   DO i = 1, ntheta
   DO j = 1, ntheta
   kinetic (i, j) = dafp (i - j, 2) - dafp (i + j - 2, 2)
   ENDDO

   ENDDO
   DO i = 1, ntheta
   kinetic (i, i) = kinetic (i, i) + 2.d0

   ENDDO
   DO i = 2, ntheta - 1
   ctsq = (cos (2.0d0 * theta (i) ) / sin (2.0d0 * theta (i) ) ) **2
   ctsq = 0.d0
   kinetic (i, i) = kinetic (i, i) + ctsq

   ENDDO
   DO j = 1, ntheta
   ctsq = (dafp (1 - j, 2) - dafp (1 + j - 2, 2) ) / 8.d0
   kinetic (1, j) = kinetic (1, j) + ctsq

   ENDDO
   DO j = 1, ntheta
   ctsq = (dafp (ntheta - j, 2) - dafp (ntheta + j - 2, 2) ) / 8.d0
   kinetic (ntheta, j) = kinetic (ntheta, j) + ctsq

   ENDDO
   DO k = 1, ntheta
   kinetic (1, k) = kinetic (1, k) / sqrt (2.0d0)
   kinetic (k, 1) = kinetic (k, 1) / sqrt (2.0d0)
   kinetic (ntheta, k) = kinetic (ntheta, k) / sqrt (2.0d0)
   kinetic (k, ntheta) = kinetic (k, ntheta) / sqrt (2.0d0)
   ENDDO


ENDIF
DO i = 1, ntheta
DO j = 1, ntheta
kinetic (i, j) = - 4.d0 * kinetic (i, j) / usys2
ENDDO

ENDDO
IF(dafopt==11) dafopt = 1
IF(dafopt==12) dafopt = 1
IF(dafopt==111) dafopt = 1
IF(dafopt==214) dafopt = 1
IF(dafopt==314) dafopt = 1
CALL TestThetaDAF (ntheta, kinetic, theta, Dafp, nd_daf)

CALL TestXDAF (ntheta, kinetic, theta, Dafp, nd_daf)
IF(debug)THEN
   WRITE(Out_Unit, * ) 'Theta_Kinetic Energy Matrix'
   CALL MxOut(kinetic, ntheta, ntheta)

ENDIF
WRITE(Out_Unit, * ) 'END(Kin_Theta)'

WRITE(Out_Unit, * )
RETURN
ENDSUBROUTINE kin_theta
