SUBROUTINE anchicp(tchi,ptchi,nchi,nchimax,antotchi,chij,debug)
USE Numeric_Kinds_Module
!
!     SUBROUTINE angmom calculates the DVR angular momentum matrices
!     j_{chi}^2 and puts the result in antotchi
!     for use in calculating the full hamiltonian for a(p)rime rep.
!     note that the full angular momentum array is calculated rather than
!     just calculating the elements for the truncated DVR.
!
!     tchi=chi dvr transformation matrix from tdvf2b
!     ptchi=chi dvr points
!     nchi=number of angular basis functions in chi
!     nchimax=maximum number of angular basis functions
!     antotchi=transformed j^2 matrix
!     chij=REAL(Kind=WP_Kind) scratch array of dimension nchimax
!
LOGICAL debug
EXTERNAL mprint
INTEGER nchi, nchimax, nchi2, j, jchi, i, k
REAL(Kind=WP_Kind) ptchi(nchi), antotchi(nchimax,nchi)
REAL(Kind=WP_Kind) tchi(nchimax,nchi)
REAL(Kind=WP_Kind) chij(nchi), temp
!
!     calculate FBR j_{chi}^2 (which is diagonal) and placing
!     results in chij for cos(2*j*chi)
!
nchi2=nchi/2
DO j=1,nchi2
!        jchi=2*(j-1)
   jchi=(j-1)
   chij(j)=jchi*jchi
ENDDO
!
!     calculate FBR j_{chi}^2 (which is diagonal) and placing
!     results in chij for sin(2*j*chi)
!
DO j=nchi2+1,nchi
!        jchi=(2*(j-nchi2))
   jchi=((j-nchi2))
   chij(j)=jchi*jchi
ENDDO
!
!     now transform the j^2 matrix using tchi
!
DO i=1,nchi
   DO j=1,nchi
      temp=0.0d0
      DO k=1,nchi
         temp=temp+tchi(i,k)*chij(k)*tchi(j,k)
      ENDDO
      antotchi(i,j)=temp
   ENDDO
ENDDO
IF(debug) CALL printm(mprint,8,"angt",antotchi,nchi  ,nchi  ,nchimax)
ENDSUBROUTINE anchicp
