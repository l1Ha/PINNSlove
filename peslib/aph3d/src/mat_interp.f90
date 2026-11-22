SUBROUTINE mat_interp(rho, n, w, w0, w1, w2, xksq, energy, rho0, rho1, rho2)
!
USE TotalEng_Module
USE Masses_Module
USE fuzzy_MODULE
USE centerho_MODULE
USE FileUnits_MODULE
IMPLICIT NONE
LOGICAL little, medium, full
INTEGER n, nrho, ms, ithcll, ithsub, i, nmodes, j, ij, maxrho
REAL(Kind=WP_Kind) rhosave, rho, rho0, rho1, rho2, a0, a1, a2
REAL(Kind=WP_Kind) abs, energy(n), xksq(n), w(n*(n+1)/2), w0(n,n), w1(n,n)
REAL(Kind=WP_Kind) w2(n,n)

!----------------------------------------------------------------------
!  written by r. t pack. this routine transforms the wigner
DATA nrho/1/, ms/1/
!----------------------------------------------------------------------
! rho       hyperpsherical radius
! n         number of coupled equations
! w         potential matrix evaluated at rho. the potential matrix is
!           determined by a lagrange interpolation of w0, w1, w2.
! w0        stores potential matrix at rho-0
! w1        stores potential matrix at rho-1
! w2        stores potential matrix at rho-2
!-----------------------------------------------------------------------
!
! start of a new sector read in values of the potential.
!-----------------------------------------------------------------------
DATA ithcll/0/,ithsub/0/
DATA little/.false./,medium/.false./,full/.false./
CALL popt('mat_interp ',little,medium,full,ithcll,ithsub)

!----------------------------------------------------------------------
! set lagrange constants to determine w at rho.
!----------------------------------------------------------------------
IF(ABS(rho1-rho2)<fuzz)THEN
a0=(rho-rho1)/(rho0-rho1)
a1=(rho-rho0)/(rho1-rho0)
a2=0.0d0
ELSEIF(ABS(rho0-rho1)<fuzz)THEN
a0=0.0d0
a1=(rho-rho2)/(rho1-rho2)
a2=(rho-rho1)/(rho2-rho1)
ELSE
a0=(rho-rho1)*(rho-rho2)/(rho0-rho1)/(rho0-rho2)
a1=(rho-rho0)*(rho-rho2)/(rho1-rho0)/(rho1-rho2)
a2=(rho-rho0)*(rho-rho1)/(rho2-rho0)/(rho2-rho1)
ENDIF
!      WRITE(Msg_Unit,*)'a012=',a0,a1,a2
!      STOP
IF(full)THEN
    WRITE(Out_Unit,80)maxrho,nrho,rho,rho0,rho1,rho2
    WRITE(Out_Unit,85)a0, a1, a2
    WRITE(Out_Unit,86)energy
    WRITE(Out_Unit,81)rho0
    CALL MxOut(w0,n,n)
    WRITE(Out_Unit,82)rho1
    CALL MxOut(w1,n,n)
    WRITE(Out_Unit,83)rho2
    CALL MxOut(w2,n,n)
ENDIF
!-----------------------------------------------------------------------
! interpolate to find w.
!-----------------------------------------------------------------------
ij=0
DO j=1,n
   DO i=1,j
      ij=ij+1
      w(ij)=a0*w0(i,j)+a1*w1(i,j)+a2*w2(i,j)
   ENDDO
ENDDO
IF(full)THEN
    WRITE(Out_Unit,'(" w after the interpolation:")')
    CALL MxOutD(w,n,n,1)
ENDIF
DO i=1,n
   xksq(i)=usys2*(Etot-energy(i))
ENDDO

ij=0
DO j=1,n
   DO i=1,j
      ij=ij+1
      IF(i==j)THEN
         w(ij)=-usys2*(-(rhocnt/rho)**2*energy(i))+usys2*w(ij)
      ELSE
         w(ij)=usys2*w(ij)
      ENDIF
   ENDDO
ENDDO


IF( (rho < rho0 .AND. ABS(rho-rho0) > fuzz)   &
     .or. (rho > rho2 .AND. ABS(rho-rho2) > fuzz) )THEN
    WRITE(Out_Unit,*)'error in mat_interp'
    WRITE(Out_Unit,*)'rho=',rho
    WRITE(Out_Unit,*)'rho0,rho1,rho2 =', rho0,rho1,rho2
    WRITE(Msg_Unit,*)'error in mat_interp'
    WRITE(Msg_Unit,*)'rho=',rho
    WRITE(Msg_Unit,*)'rho0,rho1,rho2 =', rho0,rho1,rho2
    STOP 'Stopping in mat_interp'
ENDIF

IF(full)THEN
    WRITE(Out_Unit,84)rho
    CALL MxOutD(w,n,n,1)
ENDIF

RETURN
!
80  FORMAT(//,'maxrho,nrho,rho,rho0,rho1,rho2: ',2i5,4f10.5)
81  FORMAT(//,'rho0=',f10.5)
82  FORMAT(//,'rho1=',f10.5)
83  FORMAT(//,'rho2=',f10.5)
84  FORMAT(//,'rho=',f10.5)
85  FORMAT(//,'a0, a1, a2=',3f10.5)
86  FORMAT(//,'energy=',20f10.5)
END
