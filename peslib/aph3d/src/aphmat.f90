SUBROUTINE aphmat(rho, n, w, w0, w1, w2, jderiv)
USE FileUnits_MODULE
USE Narran_Module
USE nstate_MODULE
USE TotalEng_Module
USE Masses_Module
USE Qall_Module
USE fuzzy_MODULE
USE centerho_MODULE
IMPLICIT NONE
LOGICAL little, medium, full
INTEGER n, mx, ithcll, ithsub, jderiv, nmodes, i, j,ij, maxrho, ms
INTEGER nrho
REAL(Kind=WP_Kind) rholast, rhosave, rho, rho0, rho1, rho2, a0, a1, a2
REAL(Kind=WP_Kind) energy(nstate)
REAL(Kind=WP_Kind) rhoval(21)
REAL(Kind=WP_Kind) w(n,n),w0(n,n),w1(n,n),w2(n,n)

!----------------------------------------------------------------------
!  written by r. t pack. this routine transforms the wigner
DATA rhoval/21*-1.0d0/,nrho/1/,ms/1/
!----------------------------------------------------------------------
! rho       hyperpsherical radius
! n         number of coupled equations
! w         potential matrix evaluated at rho. the potential matrix is
!           determined by a lagrange interpolation of w0, w1, w2.
! jderiv     calculates the j-deruvatuve of w at rho.
! w0        stores potential matrix at rho-0
! w1        stores potential matrix at rho-1
! w2        stores potential matrix at rho-2
!-----------------------------------------------------------------------
!
! start of a new sector read in values of the potential.
!-----------------------------------------------------------------------
DATA rholast/-1.0d0/
DATA ithcll/0/,ithsub/0/
DATA little/.false./,medium/.false./,full/.false./
CALL popt('aphmat  ',little,medium,full,ithcll,ithsub)

IF(rholast==-1.d0)THEN
    rhosave = rho
ELSEIF(rho==rhosave)THEN
    rholast = -1.d0
ENDIF
!----------------------------------------------------------------------
! Read in potential Matrix elements.
!----------------------------------------------------------------------
!      WRITE(Out_Unit,*)'n,rho',n,rho,rho0,rho1,rho2
CALL aphREAD(rho, n, w, w0, w1, w2, jderiv, nrho, nmodes, rhoval, energy, rho0, rho1, rho2, rholast)
!      WRITE(Out_Unit,*)'n,rho',n,rho,rho0,rho1,rho2
!----------------------------------------------------------------------
! set lagrange constants to determine w at rho.
!----------------------------------------------------------------------
IF(jderiv<=0)THEN
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
ENDIF
!----------------------------------------------------------------------
! set lagrange constants to determine dw/drho at rho.
!----------------------------------------------------------------------
IF(jderiv==1)THEN
    a0=(2*rho-rho1-rho2)/(rho0-rho1)/(rho0-rho2)
    a1=(2*rho-rho0-rho2)/(rho1-rho0)/(rho1-rho2)
    a2=(2*rho-rho0-rho1)/(rho2-rho0)/(rho2-rho1)
ENDIF
!----------------------------------------------------------------------
! set lagrange constants to determine d2w/drho2 at rho.
!----------------------------------------------------------------------
IF(jderiv==2)THEN
    a0=2/(rho0-rho1)/(rho0-rho2)
    a1=2/(rho1-rho0)/(rho1-rho2)
    a2=2/(rho2-rho0)/(rho2-rho1)
ENDIF
!-----------------------------------------------------------------------
! interpolate to find the j-th derivative of w.
!-----------------------------------------------------------------------
IF(jderiv>0)THEN
    DO j=1,n
       DO i=1,n
          w(i,j)=a0*w0(i,j)+a1*w1(i,j)+a2*w2(i,j)
       ENDDO
    ENDDO
    DO i=1,n
       xksq3(i)=usys2*(Etot-energy(i))
       Energy3_APH(i)=energy(i)
       xksq3_APH(i)=xksq3(i) 
    ENDDO
    DO j=1,n
       DO i=1,n
          IF(i==j)THEN
             w(i,j)=usys2*(-(rhocnt/rho)**2*energy(i)) -usys2*w(i,j)
          ELSE
             w(i,j)=-usys2*w(i,j)
          ENDIF
      ENDDO
   ENDDO
ELSE
   ij=0
   DO j=1,n
      DO i=1,j
        ij=ij+1
        w(ij,1)=a0*w0(i,j)+a1*w1(i,j)+a2*w2(i,j)
      ENDDO
   ENDDO
   DO i=1,n
       xksq3(i)=usys2*(Etot-energy(i))
       Energy3_APH(i)=energy(i)
       xksq3_APH(i)=xksq3(i) 
   ENDDO
   ij=0
   DO j=1,n
      DO i=1,j
         ij=ij+1
         IF(i==j)THEN
            w(ij,1)=-usys2*(-(rhocnt/rho)**2*energy(i)) +usys2*w(ij,1)
         ELSE
            w(ij,1)=usys2*w(ij,1)
         ENDIF
      ENDDO
   ENDDO
ENDIF

rholast=rho

IF( (rho < rho0 .AND. ABS(rho-rho0) > fuzz)  &
    .or. (rho > rho2 .AND. ABS(rho-rho2) > fuzz) )THEN
    WRITE(Out_Unit,*)'error in aphmat'
    WRITE(Out_Unit,*)'rho,rho0,rho1,rho2 = ',rho,rho0,rho1,rho2
    STOP 'aphmat'
ENDIF

IF(full)THEN
    WRITE(Out_Unit,80)maxrho,nrho,rho,rho0,rho1,rho2
    WRITE(Out_Unit,81)rho0
    CALL MxOut(w0,n,n)
    WRITE(Out_Unit,82)rho1
    CALL MxOut(w1,n,n)
    WRITE(Out_Unit,83)rho2
    CALL MxOut(w2,n,n)
    WRITE(Out_Unit,84)rho
    IF(jderiv<=0)THEN
    CALL MxOut(w,n,n)
    ELSE
    CALL MxOut(w ,n,n)
    ENDIF
ENDIF

RETURN

80  FORMAT(//,'maxrho,nrho,rho,rho0,rho1,rho2: ',2i5,4f10.5)
81  FORMAT(//,'rho0=',f10.5)
82  FORMAT(//,'rho1=',f10.5)
83  FORMAT(//,'rho2=',f10.5)
84  FORMAT(//,'rho=',f10.5)
85 FORMAT(//,'nrho is less than 3')
ENDSUBROUTINE aphmat
