SUBROUTINE EigFun(nfuns,nhermt,elevls,chinuj,xpth,wpth,jang,karran,re,rx,alpha,pot,chi)
USE FileUnits_Module, Only: Efun_Unit              ! Unit used for printing eigfun.txt
USE Numeric_Kinds_Module, Only: WP_Kind, IW_Kind   ! Kind parameters for reals and integers
USE Convrsns_Module, ONLY: autoev                 ! Converts Ha to eV
IMPLICIT NONE
INTEGER(Kind=IW_Kind) NFuns, NHermt, jang, KArran, i, ifuns
REAL(Kind=WP_Kind)    PotBC, re, rx, alpha, rre
REAL(Kind=WP_Kind)    elevls(nfuns), chinuj(nhermt,nfuns), xpth(nhermt)
REAL(Kind=WP_Kind)    wpth(nhermt), pot(nhermt), chi(nhermt)
EXTERNAL PotBC  ! Diatomic potential
!-----------------------------------------------------------------------
! This routine prints eigenfunctions of isolated diatoms
! and diatomic potentials for subsequent plotting.
!-----------------------------------------------------------------------
WRITE(Efun_Unit,*)
WRITE(Efun_Unit,*)
WRITE(Efun_Unit,'(1x,"Arrangement channel =",i5)')karran
WRITE(Efun_Unit,'(1x,"No. of functions =",i5)')nfuns
WRITE(Efun_Unit,'(1x,"Angular momentum quantum number =",i5)')jang
WRITE(Efun_Unit,'(1x,"No. of quadrature points =",i5)')nhermt
WRITE(Efun_Unit,'("Distances in bohr")')
WRITE(Efun_Unit,'(1x,10e15.7)')(re*(xpth(i)/alpha+rx),i=1,nhermt)

! Calculate the isolated diatomic potential
DO i=1,nhermt
   rre=re*(xpth(i)/alpha+rx)
   pot(i)=PotBC( rre )
ENDDO

! Print the potential in Ha
WRITE(Efun_Unit,'("Potential in Ha")')
WRITE(Efun_Unit,'(1x,10es15.7)')(pot(i),i=1,nhermt)

! Print the potential in eV
WRITE(Efun_Unit,'("Potential in eV")')
WRITE(Efun_Unit,'(1x,10es15.7)')(pot(i)*autoev,i=1,nhermt)

! For each eigenstate eled with tau, nu, j, print the energy and the weighted eigenfunction
DO ifuns=1,nfuns
   WRITE(Efun_Unit,*)
   WRITE(Efun_Unit,'(1x,"tau, nu, j=",3i5," energy(Ha)=",ES18.10," energy(eV)=",ES18.10)')karran, ifuns-1, jang, elevls(ifuns), elevls(ifuns)*autoev
   DO i=1,nhermt
      chi(i)=chinuj(i,ifuns)*sqrt(wpth(i))  ! Weighted wigenfunction
   ENDDO
   WRITE(Efun_Unit,*)'Weighted wavefunctions=chinuj(i,ifuns)*sqrt(wpth(i))'
   WRITE(Efun_Unit,'(1x,10es15.7)')(chi(i),i=1,nhermt)
ENDDO

RETURN
ENDSUBROUTINE EigFun
