SUBROUTINE potmatn4(jrot,lorb,nvib,dt1,dt2,upsln,      &
               vleg,nchan,work,bfpotmat,sfpotmat,jmax, &
               xpt, weight, rho, vprime, pj, jtot,     &
               j_min, symmetry)
USE fileunits_Module
USE numbers_Module
USE int1d_Module
USE Masses_Module
USE ngrid4_Module
USE readWRITE_Module
USE Time_Module

! This routine calculates space-fixed and body-fixed matrix
! elements of the potential energy.

IMPLICIT NONE
CHARACTER(LEN=11) :: Blank
CHARACTER(LEN=21), PARAMETER:: ProcName='PotMatn4'
LOGICAL symmetry
INTEGER(KIND=IW_Kind) nchan, n, np, jmax, jtot, j_min
INTEGER(KIND=IW_Kind) jrot(nchan), lorb(nchan), nvib(nchan)
INTEGER(KIND=IW_Kind) ibtheta, iltheta
REAL(Kind=WP_Kind) dt1, dt2, upsln(nltheta+1,nchan), ltheta, btheta
REAL(Kind=WP_Kind) work(nltheta+1), bfpotmat(nchan,nchan,0:jmax), acos
REAL(Kind=WP_Kind) sfpotmat(nchan,nchan), vprime(nbtheta,nltheta+1)
REAL(Kind=WP_Kind) pj(nbtheta,0:jmax,0:jmax), vtemp
REAL(Kind=WP_Kind) xpt(nbtheta), weight(nbtheta), rho, potdelves4
REAL(Kind=WP_Kind) vleg(nltheta+1,0:jmax,0:jmax,0:jmax)
REAL(Kind=WP_Kind) TimeInfo
LOGICAL little, medium, full
INTEGER ithcll, ithsub
data ithcll/0/,ithsub/0/
data little/.true./, medium/.false./, full/.false./

! calculate the Gauss_Legendre points and weights. A Gauss_Legendre
! integration will be used for the big-theta integrals over the
! associated legendre functions.

t1 = TimeInfo(ProcName,'Glegen')
CALL glegen(nbtheta, xpt, weight, -one, one)
t2 = TimeInfo(procname,'Glegen')
CALL WRITEtime(procname,t2-t1)

! WRITE out points and weights.

IF(medium)THEN
   WRITE(Out_Unit,*)'Legendre points and weights' 
   DO ibtheta=1,nbtheta
      WRITE(Out_Unit,20)ibtheta,xpt(ibtheta),acos(xpt(ibtheta)),weight(ibtheta)
   ENDDO
ENDIF
IF(read_pot)THEN 
   WRITE(Msg_Unit,*) 'reading sfpotmat'
   READ(readpot)sfpotmat 
   CALL MatrixOut(sfpotmat, NChan, NChan,'SFPotMat', 'SFPotMat-read',   &
                       Blank, Blank,.False., Blank, .false.)
   RETURN
ENDIF 
! evaluate the v-prime (see rtp98/4/30-3 Eq. 20) potential
! energy at the quadrature points.
! This is a potential special to Delves' sector adiabatic CC eqs.
!
!      WRITE(msg_unit,*)' In potmatn, getting vprime'
DO iltheta=1,nltheta+1
!         IF(iltheta<=n1)THEN
   IF(iltheta<=n1+1)THEN         
      ltheta=lthetzro+(iltheta-1)*dt1
   ELSE
!            ltheta=lthetzro+(n1-1)*dt1+(iltheta-n1)*dt2
      ltheta=lthetzro+n1*dt1+(iltheta-n1-1)*dt2
   ENDIF
   DO ibtheta=1,nbtheta
      btheta=acos(xpt(ibtheta))
      IF(ABS(BTheta)>4.d0*Atan(1.d0))THEN
         WRITE(Out_Unit,*)'Btheta out of range',btheta
         STOP 'Btheta out of range'
      ENDIF
      vprime(ibtheta,iltheta)=potdelves4(rho,ltheta,btheta)*usys2
!  what follows are special cases for tests
!            vprime(ibtheta,iltheta)=(1.0d-02)*usys2/autoev
!  note:  vprime=0 is not the same as a zero interaction potential;
!  for that, see ne2hsur.
!            vprime(ibtheta,iltheta) = zero
!            vprime(ibtheta,iltheta) = 1.0d0*usys2
!            vtemp = vprime(ibtheta,iltheta)/usys2
   ENDDO
!         WRITE(Out_Unit,*)' ltheta, vprime=',ltheta,vtemp
ENDDO
! evaluate the associated legendre polynomials at the quadrature
! points.
t1 = TimeInfo(ProcName,'Legndr')
CALL legndr(jmax, nbtheta, pj, xpt)
t2 = TimeInfo(procname,'Legndr')
CALL WRITEtime(procname,t2-t1)

! integrate over the big-theta coordinate store results in vleg.
t1 = TimeInfo(ProcName,'VLeg_Mat')
CALL vleg_mat(pj,vprime,weight,vleg,jmax,jtot,j_min,symmetry)
t2 = TimeInfo(procname,'Vleg_Mat')
CALL WRITEtime(procname,t2-t1)
! WRITE it if desired.
!      WRITE(Out_Unit,*)'vleg'
!      WRITE(Out_Unit,*) (vleg(iltheta,0,0,0),iltheta=1,nltheta+1)      
! calculate the body-fixed matrix elements by integrating with
! respect to little-theta. Store results in bfpotmat.
t1 = TimeInfo(ProcName,'BFpmatrix')
CALL bfpmatrix(jrot,lorb,nvib,dt1,dt2,upsln,vleg,nchan,work,bfpotmat,jmax,jtot)
t2 = TimeInfo(procname,'BFpmatrix')
CALL WRITEtime(procname,t2-t1)

IF(full)THEN
   WRITE(Out_Unit,*)'bfpotmat'
   DO n=1,nchan
      WRITE(Out_Unit,10)(bfpotmat(np,n,0),np=1,nchan)
   ENDDO
10 FORMAT(1x,5e18.10)
ENDIF

! transform from the body-fixed potential energy matrix to the
! space-fixed potential energy matrix. store results in sfpotmat.
t1 = TimeInfo(ProcName,'SFpmatrix')
CALL sfpmatrix(jtot,jrot,lorb,nvib,nchan,bfpotmat,sfpotmat,jmax)
t2 = TimeInfo(procname,'SFpmatrix')
CALL WRITEtime(procname,t2-t1)

IF(full)THEN
   WRITE(Out_Unit,*)'sfpotmat'
   DO n=1,nchan
      WRITE(Out_Unit,10)(sfpotmat(np,n),np=1,nchan)
   ENDDO
ENDIF
! bfpotmat contains body-fixed potential energy matrix elements.
! sfpotmat contains space-fixed potential energy matrix elements.          

20 FORMAT(1x,i5,3e18.10)

WRITE(readpot)sfpotmat
IF(Full)THEN
   CALL MatrixOut(sfpotmat, NChan, NChan,'SFPotMat', 'SFPotMat-wr',   &
                       Blank, Blank,.False., Blank, .false.)
ENDIF


RETURN
ENDSUBROUTINE Potmatn4
