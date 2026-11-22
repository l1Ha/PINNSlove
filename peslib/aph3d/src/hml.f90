SUBROUTINE hml (nrho, rho_val, ntheta, tht_val, nchi, nchir, &
 chi_val, t_rho, t_theta, t_chi, t_chir, t_chi1, t_chir1, &
 v_pot, nirrep, irrep, n1, nsize, jtot, parity, usys2, a, b, c, d, &
 rhom2, tblock)
USE Numeric_Kinds_Module
USE Convrsns_Module
USE FileUnits_Module
USE ONED_Module
USE DVR_Module
USE TWOD_Module
USE Spectro_Module
USE Arrch_Module
USE Oops_Module
USE parms_MODULE
IMPLICIT NONE
LOGICAL :: debug,once
INTEGER :: nirrep
INTEGER :: nrho, ntheta, nchi, nchir, n1 (nirrep), irho, itheta, &
 ichi, irrep, i, j, nsize, Indxr_New, jrho, jtheta, jchi, nbegin, &
 nlast, jrrep, jtot, lam, lammin, lammax, parity,jump,jumpmin,ibug

REAL(Kind=WP_Kind) :: t_rho (nrho, nrho), t_theta (ntheta, ntheta)
REAL(Kind=WP_Kind) :: t_chi (nchi, nchi), t_chir (nchir, nchir)
REAL(Kind=WP_Kind) :: v_pot (nrho, ntheta, nchi), rho_val (nrho)
REAL(Kind=WP_Kind) :: tht_val (ntheta), chi_val (nchi + 1)
REAL(Kind=WP_Kind) :: t_chi1 (nchi, nchi), t_chir1 (nchir, nchir)
REAL(Kind=WP_Kind) :: usys2, a (ntheta), b (ntheta), c (ntheta), d (ntheta)
REAL(Kind=WP_Kind) :: rhom2 (nrho), afac, afactor, cfac, cfactor, v_eff
REAL(Kind=WP_Kind) :: v0, v1, v2, v3, tblock (nchi, nchir)
REAL(Kind=WP_Kind),ALLOCATABLE :: hamil(:,:),sd(:),subd(:),vector(:,:)
INTEGER,ALLOCATABLE:: iwork(:),ifail(:)

INTEGER :: nevtot,nbound,jval,ngri,Npro,Npro1,Npro2,k,iupper,jupper,info,kindex,ntot,ii,m,mm,nbig2
INTEGER :: ntot2,jj
REAL(Kind=WP_Kind) :: potmax,pev2dmax,eiga1,prefac,cpu_begin, cpu_end,pev3dmax,temp,temp2
INTEGER :: nchan,indexi,indexj,nbig,kk,krho,kmax
INTEGER,PARAMETER::noscil=5
INTEGER,PARAMETER::nhermt=25
INTEGER,PARAMETER::nulow=0
INTEGER,PARAMETER::nuhigh=4
INTEGER,PARAMETER::scheme=1
INTEGER,PARAMETER::NBasisFac=1
REAL(Kind=WP_Kind) :: hamil_temp(noscil,noscil),vect(noscil,noscil)  !Several of these arrays are not used: tmpmod gregparker
REAL(Kind=WP_Kind) :: chinuj(nhermt,noscil),energy(noscil),alpha,rx
REAL(Kind=WP_Kind) :: xpth(nhermt),wpth(nhermt),re,e(noscil),h(noscil)
CHARACTER(LEN=1)lamlabel
CHARACTER(LEN=22)label
CHARACTER(LEN=4)Glabel
CHARACTER(LEN=19)Glabel2
WRITE(Out_Unit, * )
WRITE(Out_Unit, * ) 'Begin(hml)'
debug = .true.
!--------------------------------------------------------------------
! Determine Minimum and Maximum lambda projection.
!--------------------------------------------------------------------
IF(parity==0)THEN
   lammin = 0
ELSE
   lammin = 1
ENDIF

lammax = jtot
!--------------------------------------------------------------------
! WRITE input PARAMETERs
!--------------------------------------------------------------------
IF(debug)THEN
   WRITE(Out_Unit, * ) 'nrho=', nrho
   WRITE(Out_Unit, * ) 'ntheta=', ntheta
   WRITE(Out_Unit,  * ) 'nchi=', nchi, '   nchir=', nchir
   WRITE(Out_Unit, * ) 'nsize=', nsize
   WRITE(Out_Unit,  * ) 'nirrep=', nirrep, '   irrep=', irrep
   WRITE(Out_Unit, * ) 'n1=', n1
   WRITE(Out_Unit, * ) 'jtot=', jtot
   WRITE(Out_Unit, * ) 'parity=', parity
   WRITE(Out_Unit,  * ) 'lammin=', lammin, '   lammax=', lammax
ENDIF

!--------------------------------------------------------------------
! Store
!--------------------------------------------------------------------
nbegin = 1
DO jrrep = 1, irrep - 1
   nbegin = nbegin + n1 (jrrep)
ENDDO
nlast = nbegin + n1 (irrep) - 1

IF(debug)THEN
   WRITE(Out_Unit,  * ) 'nbegin', nbegin, '   nlast=', nlast
ENDIF

DO ichi = nbegin, nlast
   DO jchi = nbegin, nlast
      t_chir (ichi + 1 - nbegin, jchi + 1 - nbegin) = t_chi (ichi, jchi)
      t_chir1 (ichi + 1 - nbegin, jchi + 1 - nbegin) = t_chi1 (ichi, jchi)
   ENDDO
ENDDO

IF(debug)THEN
   WRITE(Out_Unit, * ) 't_chir'
   CALL MxOut(t_chir, nchir, nchir)
   WRITE(Out_Unit, * ) 't_chir1'
   CALL MxOut(t_chir1, nchir, nchir)
ENDIF

!--------------------------------------------------------------------
!New way to construct the hamiltonian
!--------------------------------------------------------------------
IF(.true.)THEN


!--------------------------------------------
!for debuging:
!1.one can set v=0.0d0 & t_rho=0.0d0 to the smith grand angular momentum
!2.one can set t_rho=0.0d0 to get the surface energies.
!--------------------------------------------
!v_pot=0.0d0
!t_rho=0.0d0
once=.true.

!-------------------------------------------
!Be careful with E_cut and V_cut, now we don't use it.
!-------------------------------------------
potmax=8.0d0
pev2dmax=5.0d0
pev3dmax=-1.0d-2
nstate=8
WRITE(Out_unit,*)'parameters for convergence study'
WRITE(Out_Unit,*)'Vcut=',potmax
WRITE(Out_Unit,*)'Ecut1=',pev2dmax
WRITE(Out_Unit,*)'Ecut2,Ecut3=',pev3dmax
WRITE(Out_unit,*)'nrho=',nrho
WRITE(Out_unit,*)'ntheta=',ntheta,'nchi=',nchi
WRITE(Out_unit,*)'nstate=',nstate

!--------------------------------------------
!ALLOCATE THE MATRIX
!--------------------------------------------
!ALLOCATE(e2d_new(nt2max,nrho),cc2dtot(ntheta*nchir,ntheta*nchir,nrho))
!note a new way to construct the code to save memory, it should be run twice to see the nbig
nbig2=100
WRITE(Msg_Unit,*)'newsize(cc2dtot)=',ntheta*nchir*nrho*nbig2*8/1024/1024,'MB'
WRITE(Out_Unit,*)'newsize(cc2dtot)=',ntheta*nchir*nrho*nbig2*8/1024/1024,'MB', " nbig2=",nbig2


ALLOCATE(e2d_new(nt2max,nrho),cc2dtot(ntheta*nchir,nbig2,nrho))
ALLOCATE(nc3tot(nrho),nev_3d(nrho))


IF(Jtot/=0)THEN
ALLOCATE(vectortemp(ntheta*nchir,ntheta*nchir))
ALLOCATE(cc2dtemp(ntheta*nchir,ntheta*nchir))
ALLOCATE(nev_3d1(nrho),nev_3d2(nrho))
OPEN(Unit=umat_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Hamil_Bound.bin',form='unformatted')
vectortemp=0.0d0
ENDIF

ngri=0
Npro=0
Npro1=0
Npro2=0
!-----------------------------------------------
!Start constructing the total hamiltonian
!-----------------------------------------------
DO lam=lammin,lammax
nbig=0
IF(Jtot/=0)THEN
WRITE(lamlabel,'(I1)')lam
Label='BinOut/transformation'//lamlabel
OPEN(Unit=Nzromg_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//Label,form='unformatted')
ENDIF

nev_3d=0
DO irho=1,nrho
nevtot=0
CALL cpu_time (cpu_begin)
DO itheta=1,ntheta
!-----------------------------------------------------------------------------------------------------------------
!Now we are doing a 1D surface calculation for chi
!-----------------------------------------------------------------------------------------------------------------
CALL NEWHam_chi(t_chir,nchir,itheta,b(itheta),rhom2(irho),potmax,v_pot,ntheta,nrho,irho)
   DO i=1,nc2tot(itheta)
      IF(e1d(i,itheta)*autoev<=pev2dmax)THEN
         nev(itheta)=nev(itheta)+1
      ENDIF
   ENDDO
   nevtot=nevtot+nev(itheta)
   IF(nevtot>nt2max .AND. once)THEN
      WRITE(Out_Unit,*)' ithetmax = ',itheta
      once = .false.
   ENDIF
ENDDO
CALL cpu_time (cpu_end)
   WRITE(Out_Unit, * ) 'NEWHAM_CHI: time=', cpu_end-cpu_begin
   WRITE(Msg_Unit,*)'NEWHAM_CHI: time=', cpu_end-cpu_begin

IF(ntheta/=1)THEN
CALL cpu_time (cpu_begin)
!-----------------------------------------------------------------------------------------------------------------
!Now we are doing a 2D surface calculation for theta and chi
!-----------------------------------------------------------------------------------------------------------------
CALL NEWHam_2d(rhom2(irho),ntheta,nchir,t_theta,irho,lam,jtot,a,b,c,parity,nbig2)
CALL cpu_time (cpu_end)
   WRITE(Out_Unit, * ) 'NEWHAM_2d: time=', cpu_end-cpu_begin
   WRITE(Msg_Unit,*)'NEWHAM_2d: time=', cpu_end-cpu_begin
IF(irrep==1) eiga1 = e2d_new (1,1)
   DO i=1,nc3tot(irho)
     IF(e2d_new(i,irho)*autoev<=pev3dmax)THEN
      nev_3d(irho)=nev_3d(irho)+1
     ENDIF
   ENDDO

IF(nbig<nev_3d(irho))nbig=nev_3d(irho)
WRITE(Out_unit,*)'nev_3d(irho)=',irho,nev_3d(irho),nc3tot(irho)
ELSE
IF(irrep==1) eiga1 = e1d (1,1)
ENDIF
ngri=ngri+nev_3d(irho)
ENDDO

IF(Jtot/=0)CLOSE(Nzromg_Unit)

!--------------------------------------------------------------------------
!USED to plot the surface curves
!--------------------------------------------------------------------------
!jumpmin = min(10,nchir*ntheta)
!DO jump = jumpmin, nchir*ntheta, 10
   DO irho = 1, nrho
      WRITE(e2d_new_Unit,9)rho_val(irho),e2d_new(1,irho)*autoev!  (e2d_new(i,irho)*autoev,i=max(jump-9,1),jump)
   ENDDO
!   IF(jump/=nchi*ntheta)WRITE(e2d_new_Unit,*)'&'
!ENDDO
    9 FORMAT(11e15.7)
    954 FORMAT(e15.7)
!STOP tmpmod GregParker
!--------------------------------------------------------------------------
!test for ground surface
!--------------------------------------------------------------------------
!DO irho=1,nrho
!IF(nev_3d(irho)>1)nev_3d(irho)=1
!ENDDO
!--------------------------------------------------------------------------
WRITE(Msg_Unit,*)'The number of total surfaces=',nbig
WRITE(Out_Unit,*)'The number of total surfaces=',nbig

ALLOCATE(e3d_pro(nrho,nbig),cc3dtot(nrho,nrho,nbig))
ALLOCATE(Genergy(nstate))
ALLOCATE(Gindex(nrho,nbig),Gindex2(nstate))
ALLOCATE(v_pro(nrho,nstate))

!WRITE(Msg_Unit,*)'size of V_pro=',size(v_pro)*8/1024/1024,'MB'

!--------------------------------------------------------------------------
!initialize some matrix
!--------------------------------------------------------------------------
Genergy=1.0d0
Gindex=0
Gindex2=0


IF(nrho/=1)THEN
CALL cpu_time (cpu_begin)
kindex=0
DO k=1,nbig
!-----------------------------------------------------------------------------------------------------------------
!Now we are doing a 1D calculation on ith surface to get the guess basis.
!Warning:!! Here we put ecut3=ecut2=pev3dmax, this is set in this code.
!-----------------------------------------------------------------------------------------------------------------
CALL NEWHam_3d(t_rho,rho_val,nrho,ntheta,nchir,lam,jtot,parity,a,b,d,t_chir1,npro,npro1,npro2,k,pev3dmax/autoev,kindex)
ENDDO
CALL cpu_time (cpu_end)
npro2=npro1
npro1=npro
npro=ngri
   WRITE(Out_Unit, * ) 'NEWHAM_3d: time=', cpu_end-cpu_begin
   WRITE(Msg_Unit,*)'NEWHAM_3d: time=', cpu_end-cpu_begin
ENDIF
DEALLOCATE(cc3dtot)
ENDDO

!--------------------------------------------------------------------------
!DeALLOCATE space
!--------------------------------------------------------------------------
DEALLOCATE(cc2dtot)

IF(Jtot/=0)THEN
DEALLOCATE(cc2dtemp)
DEALLOCATE(vectortemp)
ENDIF

DEALLOCATE(Gindex)
!------------------------------------------------------------
!Start to diagonalize the final hamiltonian from here
!------------------------------------------------------------
ngri=nrho

!-------------------------------------------------------------
!Our new way to do the diagonalization
!-------------------------------------------------------------
!1. let us read in all the basis and put it into Gbasis
!-------------------------------------------------------------
ALLOCATE(Gbasis(nrho*nchir*ntheta,nstate))


DO kk=1,nstate
k=Gindex2(kk)
IF(k<10)THEN
   WRITE(GLabel,'(4I1)')0,0,0,k
ELSEIF(k<100)THEN
   WRITE(GLabel,'(2I1,I2)')0,0,k
ELSEIF(k<1000)THEN
   WRITE(GLabel,'(I1,I3)')0,k
ELSEIF(k<10000)THEN
   WRITE(GLabel,'(I4)')k
ELSE
   WRITE(Out_Unit,*)'state number is greater than 9999:'
   STOP 'Diagonalization'
ENDIF
IF(k<=0)THEN
WRITE(Msg_Unit,*)'Warning! nstate=',nstate
WRITE(Msg_Unit,*)'Warning! only',kk-1,'basis'
STOP "Hml: k<=1"
ENDIF
GLabel2='BinOut/Gfunction'//Glabel
OPEN(Unit=GFunction_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//Glabel2,form='unformatted')
REWIND(GFunction_Unit)
i=0
918 READ(GFunction_Unit,END=919)temp
i=i+1
Gbasis(i,kk)=temp
GOTO 918

!-------------------------------------------------------------
!END OF READING
!-------------------------------------------------------------

919 CONTINUE
WRITE(OUT_UNIT,*)'2D surface=',kk,k
WRITE(OUT_UNIT,*)'surface energy in ev=',Genergy(kk)*autoev
CLOSE(GFunction_Unit)
ENDDO


!-------------------------------------------------------------
!2.let us put in the Hamiltonian elements now
!-------------------------------------------------------------
WRITE(Msg_Unit,*)'final size(Hamil)=',nstate*nstate/1024*8,'KB'

ntot=ntheta*nchir*nrho
ntot2=ntheta*nchir
ALLOCATE(Hamil(nstate,nstate))

!-----------------------------------
!2.1 Construct the second  part of Hamiltonian first
!------------------------------------
CALL cpu_time(cpu_begin)
DO j=1,nstate
     DO i=j,nstate
        temp=0.0d0
        DO k=1,ntot
        ii=(k-1)/ntot2+1
        temp=temp+Gbasis(k,i)*Gbasis(k,j)*v_pro(ii,j)
        ENDDO
     Hamil(i,j)=temp
!     Hamil(j,i)=Hamil(i,j)
     ENDDO
ENDDO

DEALLOCATE(v_pro)

CALL cpu_time(cpu_end)
WRITE(Msg_Unit,*)'transformation1 time=',cpu_end-cpu_begin
WRITE(OUT_UNIT,*)'transformation1 time=',cpu_end-cpu_begin
!-----------------------------------
!2.2 Construct the first part of Hamiltonian then
!------------------------------------
WRITE(Msg_Unit,*)'transmation2 to be allocated=',nrho*ntheta*nchir*nstate*8/1024/1024,'MB'

ALLOCATE(Gbasistemp(nrho*ntheta*nchir,nstate))
DO i=1,ntot
    irho=(i-1)/(ntheta*nchir)+1
    ii=i-(irho-1)*ntheta*nchir
    DO j=1,nstate
         temp=0.0d0
         DO k=1,nrho
         kk=ii+(k-1)*ntheta*nchir
         temp=temp+t_rho(irho,k)*Gbasis(kk,j)
        ENDDO
         Gbasistemp(i,j)=temp
    ENDDO
ENDDO
WRITE(Msg_Unit,*)'Done in first step in transform2'
DO j=1,nstate
   DO i=j,nstate
      temp=0.0d0
      DO k=1,ntot
      temp=temp+Gbasis(k,i)*Gbasistemp(k,j)
      ENDDO
   Hamil(i,j)=Hamil(i,j)+temp
   Hamil(j,i)=Hamil(i,j)
   ENDDO
ENDDO

DEALLOCATE(Gbasistemp)


DEALLOCATE(Gbasis)
CALL cpu_time(cpu_end)
WRITE(Msg_Unit,*)'transformation2 time=',cpu_end-cpu_begin
WRITE(OUT_UNIT,*)'transformation2 time=',cpu_end-cpu_begin

WRITE(Msg_Unit,*)'hamiltonian constructed'








!------------------------------------------------------
!2.3 chech for symmetry of the Hamiltonian
!-----------------------------------------------------
!DO i=1, nstate
!     DO j=1,i
!       IF(ABS(Hamil(i,j)-Hamil(j,i))>1.0d-8)THEN
!       WRITE(Msg_Unit,*)'symmetry of the Hamiltonian2 is wrong'
!       WRITE(Msg_Unit,*)'Hamil(',i,',',j,')=',hamil(i,j)
!       WRITE(Msg_Unit,*)'Hamil(',j,',',i,')=',hamil(j,i)
!       STOP
!       ENDIF
!     ENDDO
!ENDDO








!-------------------------------------------------------------
!3. DO the diagonalization
!-------------------------------------------------------------
GOTO  985


!------------------------------------------------------------
!READ in the final Hamiltonian
!------------------------------------------------------------
REWIND(GFunction_Unit)
WRITE(Msg_Unit,*)'size(Hamil)',ngri*ngri/1024/1024*8,'MB'
ALLOCATE(Hamil(ngri,ngri))

READ(Msg_Unit,*)ibug
WRITE(Msg_Unit,*)'ibug=',ibug

Hamil=0.0d0
READ(Msg_Unit,*)ibug
WRITE(Msg_Unit,*)'ibug=',ibug

111 READ(GFunction_Unit,END=999)indexi,indexj,temp
Hamil(indexi,indexj)=temp
Hamil(indexj,indexi)=Hamil(indexi,indexj)
GOTO 111

!END OF READING
999 CONTINUE
CLOSE(GFunction_Unit)



!------------------------------------------------------------
!3.DO the diagonalization
!------------------------------------------------------------
985 CONTINUE
ngri=nstate
!----------------------------------
!ALLOCATE SPACE TO DIAGONALIZE
!----------------------------------
ALLOCATE(sd(ngri*(ngri+1)/2))

k=0
DO  j=1,ngri
   DO  i=j,ngri
      k=k+1
      sd(k)=hamil(i,j)
ENDDO
ENDDO
DEALLOCATE(hamil)
!-------------------------------------------------
!paramenters
!-------------------------------------------------
iupper=min(4000,ngri)        !iupper is the upper limit of the output eigenvalue numbers
jupper=iupper
ALLOCATE(vector(1,1),e3d(ngri))
ALLOCATE(subd(8*ngri),iwork(5*ngri),ifail(ngri))
CALL dspevx('N','I','L',ngri,sd,0.d0,0.d0,1,iupper,0.0d-12,jupper,e3d,vector,ngri,subd,iwork,ifail,info)

DEALLOCATE(sd)
DEALLOCATE(subd,iwork,ifail,vector)


987 CONTINUE
IF(irrep==1.AND.lam==lammin) eiga1=e3d_pro(1,1)
nchan=min(800,ngri)
IF(nrho==1)THEN
         prefac = usys2 * rho_val (1) **2
         IF(ntheta==1)THEN
            prefac = usys2 * rho_val (1) **2 / (2 * b (1) )
         ENDIF
      ELSE
         prefac = usys2
ENDIF

!---------------------------------------------------------------------------
! here we compare to the VibFun_Old or VibFun_New code ground energy
! for debugiing. just comment out!
!---------------------------------------------------------------------------
! nrigid=.false.
! ioops=.false.
! karran=1
! CALL readaph
! CALL VibFun_New(NBasisFac*noscil,0,nhermt,hamil_temp,vect,chinuj,nulow,nuhigh,nuhigh-nulow+1,&
! energy,alpha,rx,xpth,wpth,re,e,h,scheme)
! energy=energy*autoev
!energy(1)=-3.850984142142656D-002   !tmpmodgregparker what in world is this for? ERROR>>>>
!---------------------------------------------------------------------------
!finish with output energy(1)
!---------------------------------------------------------------------------
      WRITE(Msg_Unit,*)'total channal number is',nchan
      WRITE(OUT_UNIT,*)'total channal number is',nchan
      nbound = nchan
      WRITE(Out_Unit, * ) 'v_pot(1),v_pot(npot),energy(1)', v_pot (1,1,1)*autoev , v_pot (nrho,ntheta,nchi)*autoev,energy(1),e2d_new

      IF(nrho/=1)THEN
      DO jval = 1, nchan
         IF(e3d(jval)>energy(1)/autoev)THEN
            nbound = jval - 1
            goto 13
         ENDIF
      ENDDO
      ELSEIF(ntheta/=1)THEN
      DO jval = 1, nchan
         IF(e2d_new(jval,1) >v_pot (1,1,1) .or.e2d_new (jval,1) >v_pot (nrho,ntheta,nchi) )THEN
            nbound = jval - 1
            goto 13
         ENDIF
      ENDDO
      ELSE
      DO jval = 1, nchan
         IF(e1d (jval,1) >v_pot (1,1,1) .or.e1d (jval,1) >v_pot (nrho,ntheta,nchi) )THEN
            nbound = jval - 1
            goto 13
         ENDIF
      ENDDO
      ENDIF
      !---------------------------------------------------------------------
      ! Print bound eigenvalues.
      !---------------------------------------------------------------------
      13    CONTINUE
      WRITE(Out_Unit, * ) 'prefac=', prefac
      WRITE(Out_Unit, * ) 'There are ', nbound, ' bound eigenvalues(hml)'
      IF(nbound/=0)THEN
         WRITE(Out_Unit, * ) 'Bound Eigenvalues in ev'
  IF(nrho/=1)THEN
         WRITE(Out_Unit, 15) 'e3d(jval)', 'e3d(jval)-e2d1', 'prefac*e3d(jval)', 'prefac*(e3d(jval)-eiga1)'
         DO jval = 1, nbound
     WRITE(Out_Unit, 30) e3d(jval)*autoev, (e3d(jval) - eiga1)*autoev, prefac * e3d(jval), prefac * (e3d(jval) - eiga1)
            WRITE(Msg_Unit,*)'The bound',jval,e3d(jval)*autoev
         ENDDO
         WRITE(Out_Unit, * )
  ELSEIF(ntheta/=1)THEN
         WRITE(Out_Unit, 15) 'e2d(jval)', 'e2d(jval)-e2d1', 'prefac*e2d(jval)', 'prefac*(e2d(jval)-eiga1)'
         DO jval = 1, nbound
     WRITE(Out_Unit, 30) e2d_new (jval,1)*autoev, (e2d_new (jval,1) - eiga1)*autoev, &
            prefac * e2d_new (jval,1), prefac * (e2d_new (jval,1)-eiga1)
         ENDDO
         WRITE(Out_Unit, * )
  ELSE
  WRITE(Out_Unit, 15) 'e1d(jval)', 'e1d(jval)-e2d1', 'prefac*e1d(jval)', 'prefac*(e1d(jval)-eiga1)'
         DO jval = 1, nbound
     WRITE(Out_Unit, 30) e1d (jval,1)*autoev, (e1d (jval,1) - eiga1)*autoev, prefac * e1d (jval,1), prefac * (e1d (jval,1) - eiga1)
         ENDDO
         WRITE(Out_Unit, * )
  ENDIF
      ENDIF
      !---------------------------------------------------------------------
      ! Print scattering eigenvalues.
      !---------------------------------------------------------------------
      IF(nchan - nbound>0)THEN
         WRITE(Out_Unit,  * ) 'There are ', nchan - nbound, ' scattering eigenvalues'
         WRITE(Out_Unit, * ) 'Scattering Eigenvalues in ev'
  IF(nrho/=1)THEN
         WRITE(Out_Unit, 15) 'e3d(jval)', 'e3d(jval)-e2d1', 'prefac*e3d(jval)', 'prefac*(e3d(jval)-eiga1)'
         DO jval = nbound+1, nchan
     WRITE(Out_Unit, 30) e3d(jval)*autoev, (e3d(jval) - eiga1)*autoev, prefac * e3d(jval), prefac * (e3d(jval) - eiga1)
         ENDDO
         WRITE(Out_Unit, * )
  ELSEIF(ntheta/=1)THEN
         WRITE(Out_Unit, 15) 'e2d(jval)', 'e2d(jval)-e2d1', 'prefac*e2d(jval)', 'prefac*(e2d(jval)-eiga1)'
         DO jval =nbound+1, nchan
     WRITE(Out_Unit, 30) e2d_new (jval,1)*autoev, (e2d_new (jval,1) - eiga1)*autoev, &
            prefac * e2d_new (jval,1), prefac * (e2d_new (jval,1)-eiga1)
         ENDDO
         WRITE(Out_Unit, * )
  ELSE
  WRITE(Out_Unit, 15) 'e1d(jval)', 'e1d(jval)-e2d1', 'prefac*e1d(jval)', 'prefac*(e1d(jval)-eiga1)'
         DO jval = nbound+1, nchan
     WRITE(Out_Unit, 30) e1d (jval,1)*autoev, (e1d (jval,1) - eiga1)*autoev, prefac * e1d (jval,1), prefac * (e1d (jval,1) - eiga1)
         ENDDO
         WRITE(Out_Unit, * )
  ENDIF
      ENDIF
       15 FORMAT(24x,A9,11x,A15,7x,A16,4x,A24)
       30 FORMAT(1x,'Eigenvalues: ',4ES23.12)
DEALLOCATE(nc3tot,nev_3d)
STOP 'oliver temp'


ELSE
!--------------------------------------------------------------------------------------------------------------------
!Only if we need the FULL Hamiltonian
!--------------------------------------------------------------------------------------------------------------------
ALLOCATE(hamil(nsize,nsize))
hamil=0.0d0
!--------------------------------------------------------------------
! Calculate diagonal Lambda blocks.
!--------------------------------------------------------------------
DO lam = lammin, lammax
   IF(debug)THEN
      WRITE(Out_Unit, * ) 'lam=', lam
   ENDIF
!--------------------------------------------------------------------
! Rho Kinetic energy term.
!--------------------------------------------------------------------
   DO irho = 1, nrho
      DO jrho = 1, nrho
         DO itheta = 1, ntheta
            DO ichi = 1, nchir
               i = Indxr_New (lam, lammin, nrho, ntheta, nchir, irho, itheta, ichi)
               j = Indxr_New (lam, lammin, nrho, ntheta, nchir, jrho, itheta, ichi)
IF(nrho/=1)  hamil (i, j) = t_rho (irho, jrho)
            ENDDO
         ENDDO
      ENDDO
   ENDDO

   IF(debug)THEN
      WRITE(Out_Unit, * ) 'After t_rho'
      CALL MxOut(hamil, nsize, nsize)
   ENDIF

!--------------------------------------------------------------------
! Chi Kinetic energy term.
!--------------------------------------------------------------------
   DO irho = 1, nrho
      DO itheta = 1, ntheta
         DO ichi = 1, nchir
            DO jchi = 1, nchir
               i = Indxr_New (lam, lammin, nrho, ntheta, nchir, irho, itheta, ichi)
               j = Indxr_New (lam, lammin, nrho, ntheta, nchir, irho, itheta, jchi)
               hamil (i, j) = t_chi (ichi, jchi) *  2.d0 * b (itheta) *rhom2 ( irho) + hamil (i, j)
            ENDDO
         ENDDO
      ENDDO
   ENDDO

   IF(debug)THEN
      WRITE(Out_Unit, * ) 'After t_chir'
      CALL MxOut(hamil, nsize, nsize)
   ENDIF

!---------------------------------------------------------------------
! Theta Kinetic energy term.
!---------------------------------------------------------------------
   DO irho = 1, nrho
      DO itheta = 1, ntheta
         DO jtheta = 1, ntheta
            DO ichi = 1, nchir
               i = Indxr_New (lam, lammin, nrho, ntheta, nchir, irho, itheta, ichi)
               j = Indxr_New (lam, lammin, nrho, ntheta, nchir, irho, jtheta, ichi)
               hamil (i, j) = hamil (i, j) + t_theta (itheta, jtheta) * rhom2 (irho)
            ENDDO
         ENDDO
      ENDDO
   ENDDO

IF(debug)THEN
   WRITE(Out_Unit, * ) 'After t_theta'
   CALL MxOut(hamil, nsize, nsize)
ENDIF

!---------------------------------------------------------------------
! Add Effective potential
!---------------------------------------------------------------------
   IF(debug)THEN
      WRITE(Out_Unit, * ) 'jtot=', jtot
      WRITE(Out_Unit, * ) 'lam=', lam
   ENDIF
   WRITE(Out_Unit,  * ) 'irrep=', irrep, '   nchir=', nchir, '   nbegin=', nbegin
   DO irho = 1, nrho
      DO itheta = 1, ntheta
         DO ichi = 1, nchir
            i = Indxr_New (lam, lammin, nrho, ntheta, nchir, irho, itheta, ichi)
            v0 = 0.d0
!            v0 = v_pot (irho, itheta, ichi)
            hamil (i, i) = hamil (i, i) + v0
            v1 = jtot * (jtot + 1) * (a (itheta) + b (itheta) ) / 2.d0 * rhom2 (irho) / usys2
!            v1 = 0.d0
            hamil (i, i) = hamil (i, i) + v1
!            v2 = 15.d0 / usys2 / 4.d0 * rhom2 (irho)
            v2 = 0.d0
            hamil (i, i) = hamil (i, i) + v2
            v3 = (2.d0 * c (itheta) - (a (itheta) + b (itheta) ) ) * lam**2 * rhom2 (irho) / usys2
            v3 = 0.d0
            hamil (i, i) = hamil (i, i) + v3
            v_eff = v0 + v1 + v2 + v3
         ENDDO
      ENDDO
   ENDDO

   IF(debug)THEN
      WRITE(Out_Unit, * ) 'After V_eff'
      CALL MxOut(hamil, nsize, nsize)
   ENDIF

ENDDO

!-------------------------------------------------------------------
! construct Coriolis coupling terms
!-------------------------------------------------------------------
DO lam = lammin, lammax - 1
   cfactor = cfac (jtot, lam, lam + 1, parity)
   DO irho = 1, nrho
      DO itheta = 1, ntheta
         DO ichi = 1, nchir
            DO jchi = 1, nchir
               i = Indxr_New (lam, lammin, nrho, ntheta, nchir, irho, itheta, ichi)
               j = Indxr_New (lam + 1, lammin, nrho, ntheta, nchir, irho, itheta, jchi)
!               hamil (i, j) = hamil (i, j) - rhom2 (irho) / usys2 * d (itheta) * t_chir1 (ichi, jchi) * cfactor
!               hamil (j, i) = hamil (i, j)
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

IF(debug)THEN
   WRITE(Out_Unit, * ) 'After Coriolis'
   CALL MxOut(hamil, nsize, nsize)
ENDIF

!-------------------------------------------------------------------
! construct Asymmetric Top coupling term for lambda=lambda'=1
! This is a diagonal contribution for only lambda=1.
!-------------------------------------------------------------------
IF(lammax>=1)THEN
   cfactor = cfac (jtot, 1, 1, parity)
   DO irho = 1, nrho
      DO itheta = 1, ntheta
         DO ichi = 1, nchir
            i = Indxr_New (1, lammin, nrho, ntheta, nchir, irho, itheta, ichi)
            hamil (i, i) = hamil (i, i) - rhom2 (irho) / usys2 * (a (itheta) - b (itheta) ) / 2.d0 * cfactor
         ENDDO
      ENDDO
   ENDDO
ENDIF

IF(debug)THEN
   WRITE(Out_Unit, * ) 'After Asymmetric Top'
   CALL MxOut(hamil, nsize, nsize)
ENDIF

!-------------------------------------------------------------------
! Construct Asymmetric Top Coupling terms.
!-------------------------------------------------------------------
DO lam = lammin, lammax - 2
   afactor = afac (jtot, lam, lam + 2, parity)
   DO irho = 1, nrho
      DO itheta = 1, ntheta
         DO ichi = 1, nchir
            i = Indxr_New (lam, lammin, nrho, ntheta, nchir, irho, itheta, ichi)
            j = Indxr_New (lam + 2, lammin, nrho, ntheta, nchir, irho, itheta, ichi)
!            hamil (i, j) = hamil (i, j) + (a (itheta) - b (itheta) ) / 2.d0 * rhom2 (irho) * afactor / usys2
!            hamil (j, i) = hamil (i, j)
         ENDDO
      ENDDO
   ENDDO
ENDDO

IF(debug)THEN
   WRITE(Out_Unit, * ) 'After Asymmetric Top (Lambda=1)'
   CALL MxOut(hamil, nsize, nsize)
ENDIF

DEALLOCATE(hamil)
ENDIF
!CALL sym_check (hamil, nsize)
WRITE(Out_Unit, * ) 'END(hml)'

WRITE(Out_Unit, * )


RETURN
ENDSUBROUTINE hml
