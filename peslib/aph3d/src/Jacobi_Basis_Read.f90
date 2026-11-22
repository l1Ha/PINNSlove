SUBROUTINE Jacobi_Basis_Read
USE FileUnits_Module, ONLY: Efun_Unit, Bas_Unit    ! Unit used for printing eigfun.txt and the binary jacobi basis unit
USE FileUnits_Common_Module, ONLY: Outdir          ! The output directory
USE Parms_Module, ONLY: Maxhermt                   ! Maximum number of Gauss-Hermite abscissas
USE Numeric_Kinds_Module, ONLY: WP_Kind, IW_Kind   ! Kind parameters for reals and integers
USE Convrsns_Module, ONLY: autoev                  ! Converts Ha to eV
USE Masses_Module, ONLY: AtomicSymbol              ! Atomic symbols
USE DiatomicPot_MODULE, ONLY: DiatomicPot          ! The diatomic molecular symbol
USE GaussB_Module, ONLY: NHermt                    ! Number of Gauss-Hermite abscissas in each arrangement channel
USE Quantb_module, ONLY: JminFull=>Jmin, JmaxFull=>Jmax, MinVib, MaxVib   !, MinVib=>MinVibFull, MaxVib=>MaxVibFull
USE GaussQuady_Module, ONLY: wpth, xpth
USE Region_Module, ONLY: iregion
IMPLICIT NONE
LOGICAL ioops
INTEGER(KIND=IW_Kind) iarran, i, jval, ip1, ip2, ip3, jmin, jmax, jmaxx, nmax
INTEGER(KIND=IW_Kind) Index, nt, nchanacc, ntmax, karran, nstat(3), nnuj(3), nnuj_nhermt, nhermt_iarran
REAL(KIND=WP_Kind) xktot, XsqHarmonic, XkHarmonic, Eint
INTEGER(KIND=IW_Kind), ALLOCATABLE:: j_Harmonic(:,:), l_Harmonic(:,:), m_Harmonic(:,:), k_Harmonic(:,:), n_Harmonic(:,:), nu_Harmonic(:,:)
REAL(KIND=WP_Kind),    ALLOCATABLE:: Psi_Harmonic(:,:,:,:), Eig_Harmonic(:,:,:), Xsq_Harmonic(:,:,:), Xk_Harmonic(:,:,:), x_Harmonic(:,:)

MaxHermt=MaxVal(NHermt)
OPEN(Unit=9945,File=TRIM(Outdir)//"DiatomicOut/Jacobi_Read.txt")
DO IARRAN=1,3
  IF(IArran==1)THEN
    Ip1=2
    Ip2=3
  ELSEIF(IArran==2)THEN
    Ip1=3
    Ip2=1
  ELSE
    Ip1=1
    Ip2=2
  ENDIF
  jmin=JminFull(0,Iarran,0)
  jmax=JmaxFull(0,Iarran)
  jmaxx=MAXVAL(JmaxFull)
  nmax=MAXVAL(MaxVib)
  DiatomicPot=TRIM(AtomicSymbol(ip1))//TRIM(AtomicSymbol(ip2))
   
! Harmonic oscillator method used in Jacobi
  IF(IARRAN==1)THEN
     OPEN(Unit=bas_unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Jac_Basis.bin', form = 'unformatted', status = 'Old')
     REWIND(Bas_Unit)
     DO karran=1,3
        READ(bas_unit)nstat(karran)
        READ(bas_unit)
        READ(bas_unit)
        READ(bas_unit)
        READ(bas_unit)
     ENDDO
     ntmax=SUM(nstat)
     REWIND(Bas_Unit)
  ENDIF
  READ(bas_unit)nt, nnuj(IARRAN), nchanacc, nnuj_nhermt, nhermt_iarran, ioops, iregion
  WRITE(9945,*)
  WRITE(9945,*)
  WRITE(9945,*)
  WRITE(9945,*)DiatomicPot, ioops, iregion
  WRITE(9945,'(7A15)')"IARRAN","nt", "nnuj(IARRAN)", "nchanacc", "nnuj_nhermt", "nhermt_iarran", "ntmax"
  WRITE(9945,'(7I15)')IARRAN, nt, nnuj(IARRAN), nchanacc, nnuj_nhermt, nhermt_iarran, ntmax
  IF(IARRAN==1)THEN
     ALLOCATE(nu_Harmonic(0:ntmax-1,1:3)) 
     ALLOCATE(j_Harmonic(0:ntmax-1,1:3))
     ALLOCATE(l_Harmonic(0:ntmax-1,1:3))
     ALLOCATE(m_Harmonic(0:ntmax-1,1:3))
     ALLOCATE(k_Harmonic(0:ntmax-1,1:3))
     ALLOCATE(x_Harmonic(0:ntmax-1,1:3))
     ALLOCATE(xsq_Harmonic(0:ntmax-1,0:jmaxx,1:3))
     ALLOCATE(xk_Harmonic(0:ntmax-1,0:jmaxx,1:3))
     ALLOCATE(n_Harmonic(0:ntmax-1,1:3))
     ALLOCATE(Psi_Harmonic(1:maxhermt,0:7,0:34,1:3))
     ALLOCATE(Eig_Harmonic(0:7,0:34,1:3))   !tmpmodGregParker
     j_Harmonic=0
     nu_Harmonic=0
     Psi_Harmonic=0.d0
     Eig_Harmonic=0.d0
  ENDIF
  READ(bas_unit)
  WRITE(9945,'(A6,ES15.7,A8,L5)')"xktot=",xktot,"  ioops=", ioops
  READ(bas_unit)(Nu_Harmonic(i,IARRAN), j_Harmonic(i,IARRAN), l_Harmonic(i,IARRAN), m_Harmonic(i,IARRAN), k_Harmonic(i,IARRAN), x_Harmonic(i,IARRAN), n_Harmonic(i,IARRAN), i=0,nt-1)
  IF(IARRAN==1)THEN
     Index=0
  ELSEIF(IARRAN==2)THEN
     Index=nnuj(1)
  ELSEIF(IARRAN==3)THEN
     Index=nnuj(1)+nnuj(2)
  ENDIF
  READ(bas_unit)(Psi_Harmonic(:, Nu_Harmonic(i,IARRAN), j_Harmonic(i,IARRAN), IARRAN), i=0,nt-1)
  WRITE(9945,'(A)')"Psi"
  READ(bas_unit)(Xsq_Harmonic(Nu_Harmonic(i,IARRAN), j_Harmonic(i,IARRAN), IARRAN), Xk_Harmonic(Nu_Harmonic(i,IARRAN), j_Harmonic(i,IARRAN),IARRAN), Eig_Harmonic(Nu_Harmonic(i,IARRAN), j_Harmonic(i,IARRAN),IARRAN), i=0,nt-1)
  DO i=0,nt-1
     Eint = Eig_Harmonic(Nu_Harmonic(i,IARRAN), j_Harmonic(i,IARRAN),IARRAN)
     WRITE(9945,'(4(A,I5),2(A12,ES23.15))')"Nbasis=", i, "  Iarran=", IARRAN, "  Nu=", Nu_Harmonic(i,IARRAN), "  j=", j_Harmonic(i,IARRAN), "  Eint(Ha)=", Eint,  "  Eint(eV)=", Eint*autoev
     WRITE(9945,'(10ES15.7)') Psi_Harmonic(:, Nu_Harmonic(i,IARRAN), j_Harmonic(i,IARRAN), IARRAN)*sqrt(wpth(:,IARRAN))
  ENDDO
  WRITE(9945,'(5A6,3A23,A7)')"Nu", "j", "l", "omega", "tau", "k^2", "k", "Eint", "nbasis"
  WRITE(9945,'(5I6,3ES23.15,I7)') (Nu_Harmonic(i,IARRAN), j_Harmonic(i,IARRAN), l_Harmonic(i,IARRAN), m_Harmonic(i,IARRAN), k_Harmonic(i,IARRAN), Xsq_Harmonic(Nu_Harmonic(i,IARRAN), j_Harmonic(i,IARRAN), IARRAN), Xk_Harmonic(Nu_Harmonic(i,IARRAN), j_Harmonic(i,IARRAN),IARRAN), Eig_Harmonic(Nu_Harmonic(i,IARRAN), j_Harmonic(i,IARRAN),IARRAN), i, i=0,nt-1)
ENDDO

CLOSE(UNIT=9945)
CLOSE(bas_Unit)

RETURN
ENDSUBROUTINE Jacobi_Basis_Read