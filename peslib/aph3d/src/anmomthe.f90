SUBROUTINE anmomthe(ttheta,pttheta,ntheta,nthetmax,ltheta,angtot,thetaj,debug)
USE Numeric_Kinds_Module
!
!     SUBROUTINE angmom calculates the DVR angular momentum matrices
!     j_{theta}^2 and puts the result in angtot
!     for use in calculating the full 2-d hamiltonian
!     thetaj is a scratch array for intermediate results.
!     note that the full angular momentum array is calculated rather than
!     just calculating the elements for the truncated DVR. Some increase
!     in speed may be achieved by only calculating the elements of angtot
!     for those angles that have gaussians on them (as defined in npttheta).
!
!     ttheta=dvr transformation matrix from tdvf2b
!     pttheta=dvr points
!     ntheta=number of angular basis functions
!     nthetmax=maximum number of angular basis functions
!     angtot=transformed j^2 matrix
!     thetaj=REAL(Kind=WP_Kind) scratch array of dimension nthetmax
!
INTEGER i, j, jtheta, k, ltheta, ntheta, nthetmax
LOGICAL debug
EXTERNAL mprint
REAL(Kind=WP_Kind) Temp
REAL(Kind=WP_Kind) ttheta(nthetmax,ntheta)
REAL(Kind=WP_Kind) pttheta(ntheta),angtot(nthetmax,ntheta)
REAL(Kind=WP_Kind) thetaj(ntheta)
!
!     begin by calculating FBR j_{theta}^2 (which is diagonal) and placing
!     results in thetaj
!
DO j=1,ntheta
   jtheta=j-1+ltheta
   thetaj(j)=16.d0*(jtheta*(jtheta+1))
ENDDO
!
!     now transform the j^2 matrix using ttheta
!
DO i=1,ntheta
   DO j=1,ntheta
      temp=0.0d0
      DO k=1,ntheta
         temp=temp+ttheta(i,k)*thetaj(k)*ttheta(j,k)
      ENDDO
      angtot(i,j)=temp
   ENDDO
ENDDO
IF(debug) CALL printm(mprint,8,"angt",angtot  ,ntheta,ntheta,nthetmax)

RETURN
ENDSUBROUTINE anmomthe
