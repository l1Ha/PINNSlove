SUBROUTINE v_READ(nrho, ntheta, nchi, v_pot, rho_val, tht_val, chi_val, usys2, read_v, nsize, system, group)
USE Numeric_Kinds_Module
USE FileUnits_Module
USE InputFile_Module
USE Masses_Module, ONLY: Mass
IMPLICIT NONE
CHARACTER(LEN=3) system, group
CHARACTER(LEN=25) PES_Name
LOGICAL vmin_rho, vmin_theta, vmin_chi, groupc1, const_rho, const_theta, const_chi, read_v, groupc2
INTEGER nrho, nchi, ntheta, irho, itheta, ichi, mrho, mchi, mtheta, nsize
REAL(Kind=WP_Kind) Usys2, rho_min, rho_max, rho, theta, chi
REAL(Kind=WP_Kind) rho_val(nrho),tht_val(ntheta), chi_val(nchi+1), v_pot(nrho,ntheta,nchi+1)
NAMELIST/junk/ groupc1, groupc2

WRITE(Out_Unit,*)
WRITE(Out_Unit,*)'Begin(V_Read) READ in Potential'
groupc1=.false.
groupc2=.false.

OPEN(Unit=In_Unit,File=InputDIR(1:LEN(TRIM(InputDIR)))//InputFile,status='old')
WRITE(Out_Unit,*)'groupc1=',groupc1
WRITE(Out_Unit,*)'groupc2=',groupc2
CLOSE(Unit=In_Unit)

OPEN(Unit=Pot_Bin_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/Pot.bin',form='unformatted')
READ(Pot_Bin_Unit)PES_Name
READ(Pot_Bin_Unit)mass,usys2
READ(Pot_Bin_Unit)rho_min,rho_max
READ(Pot_Bin_Unit)mrho,const_rho,rho,vmin_rho
READ(Pot_Bin_Unit)mchi,const_chi,chi,vmin_chi
READ(Pot_Bin_Unit)mtheta,const_theta,theta,vmin_theta

IF(.NOT.read_v)THEN
   nrho=mrho
   ntheta=mtheta
   IF(mchi/=1)THEN
      nchi=mchi-1
   ELSE
      nchi=1
   ENDIF
ENDIF

IF(mass(1)==mass(2).AND.mass(2)==mass(3))THEN
   system='AAA'
   group='C6v'
   nsize=2*nrho*ntheta*nchi/12
ELSEIF(mass(2)==mass(3))THEN
   system='ABB'
   group='C2v'
   nsize=nrho*ntheta*nchi/4
ELSE
   system='ABC'
   group='C2 '
   nsize=nrho*ntheta*nchi/2
ENDIF

IF(groupc1)THEN
   group='C1'
   nsize=nrho*ntheta*nchi
ENDIF

IF(groupc2)THEN
   group='C2'
   nsize=nrho*ntheta*nchi/2
ENDIF

IF(nchi<=1)THEN
   group='C1 '
   nsize=nrho*ntheta
ENDIF

nsize=nrho*ntheta*nchi
WRITE(Out_Unit,*)'system=',system
WRITE(Out_Unit,*)'group=',group
WRITE(Out_Unit,*)'PES_Name=',PES_Name
WRITE(Out_Unit,'(1x,A7,5es14.7)')'masses=',mass(1),mass(2),mass(3),usys2
WRITE(Out_Unit,*)'rho min and max=',rho_min,rho_max
WRITE(Out_Unit,*)'nrho=',mrho,const_rho,rho,vmin_rho
WRITE(Out_Unit,*)'nchi=',mchi,const_chi,chi,vmin_chi
WRITE(Out_Unit,*)'ntheta=',mtheta,const_theta,theta,vmin_theta
WRITE(Out_Unit,*)'nsize=',nsize

IF(read_v)THEN
   IF(nrho/=mrho)THEN
      WRITE(Out_Unit,*)'Error: nrho, mrho=',nrho,mrho
      STOP 'v_read'
   ENDIF
   IF(mchi/=nchi+1.AND.mchi/=1)THEN
      WRITE(Out_Unit,*)'Error: nchi, mchi=',nchi,mchi
      STOP 'v_read'
   ENDIF
   IF(mtheta/=ntheta)THEN
      WRITE(Out_Unit,*)'Error: ntheta, mtheta=',ntheta,mtheta
      STOP 'v_read'
   ENDIF
   DO irho=1,nrho
      DO itheta=1,ntheta
         DO ichi=1,mchi
            READ(Pot_Bin_Unit)rho_val(irho),tht_val(itheta),chi_val(ichi),v_pot(irho,itheta,ichi)
         ENDDO
      ENDDO
   ENDDO
ENDIF

CLOSE(Unit=Pot_Bin_Unit)
WRITE(Out_Unit,*)'END(V_Read)'
WRITE(Out_Unit,*)

RETURN
ENDSUBROUTINE V_Read
