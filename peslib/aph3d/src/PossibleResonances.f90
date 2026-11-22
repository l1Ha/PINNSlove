SUBROUTINE PossibleResonances(NCurves, NSectors, Curve, RhoVals)
USE Numeric_Kinds_Module
USE FileUnits_Module
USE Masses_Module
USE Convrsns_Module
IMPLICIT NONE
LOGICAL there
LOGICAL, PARAMETER:: dbug=.False.
CHARACTER(LEN=1), PARAMETER :: Jobz='N', Uplo='U'
CHARACTER(LEN=15) :: Identifier
CHARACTER(LEN=80) SubTitle
INTEGER NCurves, NSectors, i, j, kmax, Info, IOP(2), minindex(1), maxindex
INTEGER, PARAMETER:: nxval=2**6                  ! Number of distances
INTEGER, PARAMETER:: ndaf=2                      ! ndaf=2 for second order differential equations
INTEGER, PARAMETER:: nd_daf=40
INTEGER, PARAMETER:: lwork=1000
REAL(Kind=WP_Kind) hstep, xval(nxval), A(nxval), B(nxval), C(nxval), vmin, vmax
INTEGER, PARAMETER :: mval=32                    ! order of Hermite DAF expansion (always 32)
INTEGER, PARAMETER :: np=3                       ! degree of interpretation polynomial
REAL(kind=WP_Kind) dafx(0:nd_daf,0:ndaf), hermit(0:mval+ndaf+1), eig(nxval)
REAL(Kind=WP_Kind) Curve(NCurves,NSectors)
REAL(Kind=WP_Kind) RhoVals(NSectors), Work(lwork)
REAL(Kind=WP_Kind),  ALLOCATABLE:: Hamil(:,:)    ! Hamil(NSectors,NSectors)
REAL(Kind=WP_Kind),  ALLOCATABLE:: Toeplitz(:)   ! Toeplitz(NSectors)
REAL(Kind=WP_Kind),  ALLOCATABLE:: VPot(:)       ! VPot(nxval)
REAL(Kind=WP_Kind),  ALLOCATABLE:: YD(:)         ! YD(NSectors)

IOP(1)=5
IOP(2)=5

ALLOCATE(Hamil(nxval,nxval))    ! Allocate Hamil    
ALLOCATE(Toeplitz(nd_daf))      ! Allocate Toeplitz    
ALLOCATE(VPot(nxval))           ! Allocate VPot    
ALLOCATE(YD(NSectors))          ! Allocate YD


hstep=(Rhovals(NSectors)-Rhovals(1))/(nxval-1)
DO i=1,nxval
   xval(i)=Rhovals(1)+(i-1)*hstep
ENDDO


Identifier="_Resonances"
CALL get_daf(nxval, xval, hstep, mval, dafx, hermit, ndaf, nd_daf, toeplitz, UABC2, kmax, Identifier)

   
WRITE(9837,'(A)')"Possible resonances associated with diabatic curve="
WRITE(9837,'(2A11,6A17)')"Curve, ","Ith-Eig, ","Eig(Ha), ","Eig(eV), ","Vmin(eV), ","Vmax(eV), ","VAsy(eV), ","Classification"
DO j=1,NCurves
   DO i=1,NSectors
      YD(i)=Curve(j,i)
   ENDDO
   Hamil=0.D0
   CALL kin_radial(nxval, toeplitz, Hamil, kmax)
   CALL AKIMA (NP,NSectors,Rhovals,YD,nxval,xval, VPot)
   vpot=vpot/autoev
   DO i=1,nxval
      Hamil(i,i)=Hamil(i,i)+VPot(i)
   ENDDO
   CALL dsyev(Jobz, UpLo, nxval, Hamil, nxval, Eig, work, lwork, info)

   vmin=MINVAL(VPot)
   minindex=MINLOC(VPot)
   vmax=MAXVAL(VPot(Minindex(1):nxval))
   DO i=minindex(1),nxval
      IF(Vpot(i)>Curve(j,NSectors)/autoev)maxindex=i
   ENDDO
   DO i=1,nxval
      IF(Eig(i)>Curve(j,NSectors)/autoev.and.Eig(i)<vmax)THEN
         WRITE(9837,'(2(I10,","),5(ES16.7,","),A)')j,i,Eig(i),Eig(i)*autoev,vmin*autoev,vmax*autoev,Curve(j,NSectors),"Quasibound"
      ELSEIF(Eig(i)<Curve(j,NSectors)/autoev)THEN
         WRITE(9837,'(2(I10,","),5(ES16.7,","),A)')j,i,Eig(i),Eig(i)*autoev,vmin*autoev,vmax*autoev,Curve(j,NSectors),"Bound"
      ENDIF
   ENDDO
ENDDO


DEALLOCATE(Hamil)    ! Deallocate Hamil    
DEALLOCATE(Toeplitz) ! Deallocate Toeplitz    
DEALLOCATE(VPot)     ! Deallocate VPot    
DEALLOCATE(YD)       ! Deallocate YD


RETURN
ENDSUBROUTINE PossibleResonances