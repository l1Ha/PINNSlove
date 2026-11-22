module pes
use RKHS
implicit none
! The kernel type is needed to set up and evaluate a RKHS model
type(kernel), save :: pes11, pes12, pes22
logical, save :: krb2_debug = .false.
end module pes

!*************************************************************************
!program to get potential energy for a given geometry
! input parameter:
! q(1) --> R1  K1-Rb 
! q(2) --> R2  K2-Rb
! q(3) --> R3  Rb1-Rb2
! output:
!        diabatic matrix | Vd11   Vd12 |
!                        |-Vd12   Vd22 |
!       adiabatic potential  Va1, Va2
! The angle interpolation use x=(1-cos(theta))/2, theta is Rb1-K-Rb2 
!
      subroutine krb2_pes_eval(q,lambda,vd11,vd22,vd12,Va1,Va2)
      use pes
      implicit none
      double precision,dimension(:),intent(in) :: q(3)
      double precision,intent(in) :: lambda
      double precision,intent(out) :: vd11,vd22,vd12,Va1,Va2
      double precision,parameter :: R1f = 28.d0, R2f = 28.d0, R3f = 28.d0
      double precision, parameter :: au_2_cm = 219474.6313702d0
      real*8,parameter :: pi = acos(-1.0d0), piby180 = pi/180.0d0
      real*8,external :: pwpot_ab,SW
      !subroutine parameters
      real(kind(0d0)) :: R1,R2,R3,theta
      double precision :: Vip11,Vip22,Vip12

      R1 = q(1)
      R2 = q(2)
      R3 = q(3)
                  if (krb2_debug .and. (R1.GT.R1f .or. R2.GT.R2f .or. R3.GT.R3f)) then
                        write(*,*) 'DEBUG krb2_pes: out-of-range call R1=', R1, ' R2=', R2, ' R3=', R3
                  end if
      if (R1.le.0d0 .or. R2.le.0d0 .or. R3.le.0d0) then
        print*, "Abnormal exit: R1, R2, and/or R3 smaller than zero"
        stop
      end if      
      theta = (R1**2 + R2**2 - R3**2)/2.d0/R1/R2
      theta = min(1.d0,max(-1.d0,theta))
      theta = (1.d0-theta)/2.d0

                  Vip11 = 0.d0
                  Vip22 = 0.d0
                  Vip12 = 0.d0
                  Va1 = 0.d0
                  Va2 = 0.d0
                  if (R2.GE.R1) then
                        call pes11%evaluate_fast((/theta,R1,R2/),Vip11)
                        call pes22%evaluate_fast((/theta,R1,R2/),Vip22)
                        call pes12%evaluate_fast((/theta,R1,R2/),Vip12)
                  else
                        call pes11%evaluate_fast((/theta,R2,R1/),Vip11)
                        call pes22%evaluate_fast((/theta,R2,R1/),Vip22)
                        call pes12%evaluate_fast((/theta,R2,R1/),Vip12)
                        Vip12 = (-1.d0)*Vip12
                  end if
                  if (krb2_debug) then
                     write(*,*) 'DEBUG krb2_pes evaluate: Vip11=', Vip11, ' Vip22=', Vip22, ' Vip12=', Vip12
                  end if
      if (R2.eq.R1) Vip12 = 0.d0
     
      Vd11 = Vip11*SW(R1,R2,R3)*lambda + pwpot_ab(R1,R2,R3,11)
      Vd22 = Vip22*SW(R1,R2,R3)*lambda + pwpot_ab(R1,R2,R3,22) 
      Vd12 = Vip12*SW(R1,R2,R3)*lambda + pwpot_ab(R1,R2,R3,12)

      Va1 = (Vd11 + Vd22)/2 - dsqrt(0.25d0*(Vd11 - Vd22)**2 + Vd12**2)
      Va2 = (Vd11 + Vd22)/2 + dsqrt(0.25d0*(Vd11 - Vd22)**2 + Vd12**2)

      END subroutine krb2_pes_eval
!********************************************************************************
      subroutine rkhs_init
                        use pes
                        use mod_common
                        implicit none
                        logical, save :: ker1 = .false., ker2 = .false., ker3 = .false.
                  character (len=512), save :: datapath="./KRb2_PESs/peskernel/"
                  double precision, parameter :: R1f_def = 28.d0, R2f_def = 28.d0, R3f_def = 28.d0

                        ! file_exists will be true if the file exists and false otherwise
                        if (.not. ker1) then
                              inquire(file=trim(datapath)//"pes11.kernel", exist=ker1)
                        end if
                        if (.not. ker2) then
                              inquire(file=trim(datapath)//"pes12.kernel", exist=ker2)
                        end if
                        if (.not. ker3) then
                              inquire(file=trim(datapath)//"pes22.kernel", exist=ker3)
                        end if
                        
                        if (ker1 .and. ker2 .and. ker3 ) then
                              if (krb2_debug) write(*,*) 'DEBUG RKHS: Loading kernels from ', trim(datapath)
                              call pes11%load_from_file(trim(datapath)//"pes11.kernel")
                              call pes12%load_from_file(trim(datapath)//"pes12.kernel")
                              call pes22%load_from_file(trim(datapath)//"pes22.kernel")
                              ! set PES max distances so callers can pick thresholds
                              pes_rmax = (/ R1f_def, R2f_def, R3f_def /)
                        else
                              write(*,'(A,A)') "ERROR: Could not find PES kernels in ", trim(datapath)
                              stop
                        end if   
! finish the initialize RKHS
       end subroutine rkhs_init
! switching function, which can make the short-range region and long-range region
! connect smoothly.
! R1c, R2c, R3c are the switching parameters
! R1: K1-Rb
! R2: K2-Rb
! R3: Rb1-Rb2, or angle Rb1-K-Rb2
function SW(R1,R2,R3)
  implicit none
  double precision, parameter :: pi = dacos(-1.d0)
  double precision, parameter :: ALFA = -0.5d0
  double precision,parameter :: R1f = 28.d0, R2f = 28.d0, R3f = 28.d0  
  double precision, intent(in) :: R1,R2,R3
  double precision :: R11,R22,R33
  double precision :: SW

  R11 = R1
  R22 = R2
  R33 = R3
! if the R3 is the angle
!  R33 = dsqrt(R11**2+R22**2-2.d0*R11*R22*cos(R33*pi/180.d0))

  SW = 0.d0
  SW = (0.5D0+0.5D0*TANH(ALFA*(R11-R1f)))   &
       *(0.5D0+0.5D0*TANH(ALFA*(R22-R2f)))  &
       *(0.5D0+0.5D0*TANH(ALFA*(R33-R3f)))

  return
  end
! End of the interpolation subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine krb2_pes_init
      use PES      
      implicit none
      call rkhs_init     
      call syspot_init
    end subroutine krb2_pes_init  

!*******************************************************************
! initialize the pes
! read the kernels in the 'kernel' folder and prepare the coeffficients
! inilize the pairwise potential functions and read the dimer potential
! data in the 'dimer' folder
      module para_syspot
        integer :: ncheck
        save
      end module para_syspot

      subroutine syspot_init
      use para_syspot
      implicit real * 8 (a-h,o-z)
      logical :: lsetup
! The kernel type is needed to set up and evaluate a RKHS model

! initialize the syspot KRb, Rb2 functions
      lsetup = .true.
      NV=1
      aaa=SYS_SS_POT_KRb(R,V,NV,LSETUP)
      NV=2
      aaa=SYS_SS_POT_KRb(R,V,NV,LSETUP)
      NV=1
      aaa=SYS_SS_POT_Rb2(R,V,NV,LSETUP)
      NV=2
      aaa=SYS_SS_POT_Rb2(R,V,NV,LSETUP)
      ncheck=123456789
! end of the initilization of syspot part
     END
!
function pwpot_ab(r12,r13,r23,nd)
!          K
!        /~  \
!  r12  /angl \ r13
!      /       \
!     /         \
!   Rb1----------Rb2
!         r23
!
! angl-Rb1-K-Rb2
! r12 -K1Rb
! r13 -K2Rb
! r23 -Rb1Rb2
implicit real * 8 (a-h,o-z)
!pi= dacos(-1d0)
!  write(*,*) angl,r12,r13,nd
!  r23=r12**2 + r13**2 - 2*dcos(angl/180d0*pi)*r12*r13
!  r23=dsqrt(r23)
!  write(*,*) r12,r13,r23
  call KRbRbhij_ab(r12,r13,r23,h11,h22,h12,e1,e2)
!  write(*,*) h11,h22,h12
    pwpot_ab = 0.d0
    if(nd.eq.11) then
      pwpot_ab = h11
    else if(nd.eq.22) then
      pwpot_ab = h22
    else if(nd.eq.12) then
      pwpot_ab = h12
    else
      write(*,*)"Please give the correct nd!"
      stop
    end if
                                                                                        RETURN
   END 

!======================================================================
!     KRbRb pairwise potential subroutine, the dimer potentials are from
!          molpro ab initio calculations KRbKhij_ab
!======================================================================      
      subroutine KRbRbhij_ab(r12,r13,r23,h11,h22,h12,e1,e2)
        use para_syspot      
!          K 
!        /  \
!  r12  /    \ r13
!      /      \
!     /        \
!   Rb1---------Rb2
!         r23 
!
! r12 -K1Rb
! r13 -K2Rb
! r23 -Rb1Rb2

       implicit real * 8 (a-h,o-z)

       IF(ncheck.ne.123456789) then
       write(*,*) 'ncheck.ne.123456789'
       write(*,*) 'ncheck= ',ncheck
       write(*,*)'Please CALL SUBROUTINE INIT (only once)befor fit3dhij'
       STOP
       END IF

       NV=1
       aaa=SYS_SS_POT_KRb(r12,e12s,NV,.false.)
       NV=2
       aaa=SYS_SS_POT_KRb(r12,e12t,NV,.false.)

       NV=1
       aaa=SYS_SS_POT_KRb(r13,e13s,NV,.false.)
       NV=2
       aaa=SYS_SS_POT_KRb(r13,e13t,NV,.false.)

       NV=1
       aaa=SYS_SS_POT_Rb2(r23,e23s,NV,.false.)
       NV=2
       aaa=SYS_SS_POT_Rb2(r23,e23t,NV,.false.)

       h11 = 3d0/4d0*(e12s + e13s) + 1d0/4d0*(e12t + e13t) + e23t
       h22 = 1d0/4d0*(e12s + e13s) + 3d0/4d0*(e12t + e13t) + e23s
       h12 = dsqrt(3d0)/4d0*(e13s - e12s - e13t + e12t)

       e1 = (h11 + h22)/2d0 - dsqrt(0.25d0*(h11 - h22)**2 + h12**2)
       e2 = (h11 + h22)/2d0 + dsqrt(0.25d0*(h11 - h22)**2 + h12**2)
!       write(*,*) h11,h22,h12 

       END
!=========end of pairwise potential subroutine====================

!=======================================================================
!
!    Born-Oppenheimer potentials  in au
!
!          NV=1 KRb singlet  (X)
!          NV=2 KRb triplet  (a)
!
!=======================================================================
!
!     V (r) = - ( c6/r^6 + c8/r^8 + c10/r^10 )  + Vex r^beta
!     exp(-kappa*r)
!              for r to infinity   ( Note sign convention)
!
!=======================================================================
!=======================================================================
      DOUBLE PRECISION FUNCTION SYS_SS_POT_KRb(R,V,NV,LSETUP)
      IMPLICIT none
      double precision :: r,V,dddDe,pi
      integer :: nv
      logical :: lsetup
      integer, parameter :: Num_ABOS=2
!
!     name of file containing table with potential values
!
      character(LEN=40), dimension(Num_ABOS) :: datafile
      character(LEN=40) :: filename
!
!     long-range dispersion coefficient, exchange coefficients,
!     and transition point  in a.u.
!
      double precision, dimension(Num_ABOS), save :: lc6,lc8,lc10,    &
                                                     lVex,lbeta,lkappa
      double precision, dimension(Num_ABOS), save :: lfitpt_lr
      double precision :: c6,c8,c10, Vex,beta,kappa, fitpt_lr
      double precision, dimension(Num_ABOS), parameter ::       &
         c6def=(/ 4299.509938d0, 4299.509938d0 /),             &
         c8def=(/ 482340.1050d0, 482340.1050d0 /),             &
         c10def=(/ 68143258.45d0, 68143258.45d0  /),           &
         Vexdef=(/ -2.359438727471d-3, 2.359438727471d-3/),    &
         betadef=(/ 5.25669d0, 5.25669d0/),                    &
         kappadef=(/ 1.118918749d0, 1.118918749d0 /),          &
         fitpt_lrdef=(/25.0d0,25.0d0/)
      
!     short range control parameter in cm-1
!
      double precision, dimension(Num_ABOS), save :: ldDe,lRe,lRw
      double precision :: dDe, Re, Rw
      double precision, dimension(Num_ABOS), parameter ::      & 
         dDedef= (/ 0d0, 0d0/),    &     ! in cm-1                     
         Redef= (/ 12d0, 12d0/),   &     ! in a.u.
         Rwdef= (/ 3d0, 3d0/)            ! in a.u.
!
!     short range fit parameters and transition point 
!                   slope*(x-fitpt_sr) + b  in a.u.
!
      double precision,dimension(Num_ABOS),  save :: fitpt_sr
      double precision,dimension(Num_ABOS),  save :: slope,b
!
!     arrays with data points
!
!        n_sr and n_lr added to short and long range to mimimize 
!        oscillation in the interpolated potential.
!
      integer, parameter :: NMAX=300
      double precision,save, dimension(NMAX,Num_ABOS) :: xdata,ydata
      integer, dimension(Num_ABOS), save :: ntot
      integer, parameter :: n_sr=3,n_lr=5
!
!     akima control parameter 
!           Np=number of polys in spline interval 
!
      integer, parameter :: NP=3
!
!     local variables and arrays
!
!       n:  number of of points read from file
!
      integer :: n
!
      integer :: k,pot
      logical :: ios
      double precision :: rr,dr,xSmooth
      double precision :: smooth
      double precision :: longrangeXr6
      double precision :: xi(1), yi(1)
      logical,save :: debug=.false.
!
!     constants
!
      double precision,parameter :: ALPHA=137.035999139D0,  &
                                    BOHRA=0.52917721067d0,  &
                                    HARTREE=219474.6313702d0
      common /c6/ c6
!
      namelist / XRbK  / dDe,c6,c8,c10, Vex,beta,kappa,     & 
                            filename,fitpt_lr,debug,Re,Rw
      namelist / aRbK  / dDe,c6,c8,c10, Vex,beta,kappa,     &
                            filename,fitpt_lr,debug,Re,Rw
     
       pi=dacos(-1d0)
!
!=========Initial setup================================================
!

              
      IF(lsetup) THEN           ! for driverPD
!     IF(lsetup.and.nv.eq.1) THEN    ! for the close coupling code

!        do pot=1,Num_ABOs    !  for cc
          pot=nv              !  for driverPD

          c6=C6def(pot)
          c8=C8def(pot)
          c10=C10def(pot)
          Vex=Vexdef(pot)
          beta=betadef(pot)
          kappa=kappadef(pot)
          dDe=dDedef(pot)
          Re=Redef(pot)
          Rw=Rwdef(pot)
          fitpt_lr=fitpt_lrdef(pot)
                              select case(pot)
                              case(1)
                                    filename='./KRb2_PESs/dimer/XKRb.dat'
!            read(5,XRbK)
                              case(2)
                                     filename='./KRb2_PESs/dimer/aKRb.dat'      
!            read(5,aRbK)
          end select
          lC6(pot)=C6
          lC8(pot)=C8
          lC10(pot)=C10
          lVex(pot)=Vex
          lbeta(pot)=beta
          lkappa(pot)=kappa
          ldDe(pot)=dDe
          lRe(pot)=Re
          lRw(pot)=Rw
          lfitpt_lr(pot)=fitpt_lr
!
!          write(6,*) '# parameters for pot= ',pot
!          write(6,*) '#        dDe= ',dDe
!          write(6,*) '#        Vex= ',Vex
!          write(6,*) '#       beta= ',beta
!          write(6,*) '#      kappa= ',kappa
!          write(6,*) '#   fitpt_lr= ',fitpt_lr
!          write(6,*) '#         C6= ',C6
!          write(6,*) '#         C8= ',C8
!          write(6,*) '#        C10= ',C10
!          write(6,*) '#         Re= ',Re
!          write(6,*) '#         Rw= ',Rw
!
          datafile(pot)=filename
!
! Now I need to 1) Read in the potential from CCSD(T) data 
!                  The potentials are in bohr and cm-1.        !!!!!!
!               2) Fit the inner ab-initio potential and 
!                  outer analytic form in a continuous manner.
!                  This fitting occurs at fitpt
!               3) Extrapolation parameters for the really short range 
!
!
!           
!
!
!         read potentials for r< fitpt_lr, 
!         add shift dDe(cm-1), convert to atomic units, and 
!         multiply the potential by r**6.
!
!
          inquire(file=datafile(pot),EXIST=ios)
          if( .not.ios ) then
            write(6,*) 'Missing KRb2 dimer data file: ',trim(datafile(pot))
            write(6,*) 'Download the PCCP supplementary archive d6cp00282j1.zip'
            write(6,*) 'and place XKRb.dat/aKRb.dat/XRbRb.dat/aRbRb.dat under'
            write(6,*) './KRb2_PESs/dimer/ before running KRb2 calculations.'
            stop
          endif
          if (debug) then
            write(6,*) '----------------------------------'
            write(6,*) 'Initializing potential ',pot
            write(6,*) '----------------------------------'
          endif
          open(unit=20,status='old',file=datafile(pot))
          read(20,*) n
!
          k=1+n_sr
          read(20,*) xdata(k,pot),ydata(k,pot)
! simple shift
!          ydata(k,pot)=(ydata(k,pot)+ldDe(pot))*xdata(k,pot)**6/HARTREE
! arctan shift
!         dddDe=ldDe(pot)/pi*2d0
!    + *datan((lfitpt_lr(pot)-xdata(k,pot))
!    + /(lfitpt_lr(pot)-lRe(pot)))
! fermi dirac shift
          dddDe=dexp((xdata(k,pot)-lRe(pot))/lRw(pot))+1d0
            dddDe=ldDe(pot)/dddDe
          ydata(k,pot)=(ydata(k,pot)+dddDe)*xdata(k,pot)**6/HARTREE
!end pet          
          do while ( k.lt.n+n_sr .and. xdata(k,pot).lt.lfitpt_lr(pot) )
            k=k+1
            read(20,*) xdata(k,pot),ydata(k,pot)
! simple shift
!           ydata(k,pot)=(ydata(k,pot)+ldDe(pot))*xdata(k,pot)**6/
!    +                                                         HARTREE
! arctan shift
!         dddDe=ldDe(pot)/pi*2d0
!    + *datan((lfitpt_lr(pot)-xdata(k,pot))
!    + /(lfitpt_lr(pot)-lRe(pot)))
! fermi dirac shift     
            dddDe=dexp((xdata(k,pot)-lRe(pot))/lRw(pot))+1d0
            dddDe=ldDe(pot)/dddDe
            ydata(k,pot)=(ydata(k,pot)+dddDe)*xdata(k,pot)**6/HARTREE
!end pet
          end do
          close(20)
          n=k
          if ( xdata(k,pot).ge.lfitpt_lr(pot) ) n=n-1 

          ntot(pot)=n+n_lr
          if (n.le.2) then
          write(6,fmt="('Insufficient number of points in potential.' &
                       ' n=',I2,' but must larger than 2')") n
            stop
          endif
          if ( ntot(pot).gt.NMAX) then
            write(6,fmt="('Insufficient space to store the ground ', &
              'state potential. Increase NMAX to ',I4 )") ntot(pot)
            stop
          endif
!
          if(debug) then
            write(6,*) 'now it does not work properly'
            write(6,*) '# original data from ', trim(datafile(pot))
            do k=1+n_sr,n
              write(6,*) xdata(k,pot),                            &
                 ydata(k,pot)*HARTREE/xdata(k,pot)**6-ldDe(pot)
            end do
          endif
!
! DO 2) ->
!
!         replace or add to the ab initio potentials with the analytic
!         form for xdata > fitpt.  And make a smooth transition over 
!         the radial range xdata(n) and lfitpt_lr+(n_lr-1)*dr
!
          dr=0.5d0
          rr=lfitpt_lr(pot)
!ML          rr=lfitpt_lr(pot)+2*dr
          do k=n+1,ntot(pot)
             xdata(k,pot)=rr
             xSmooth=(rr-xdata(n,pot))/(lfitpt_lr(pot)-xdata(n,pot))
!            xSmooth=(rr-xdata(n,pot))/(dr*(n_lr-1))
             ydata(k,pot)=longrangeXr6(rr,lc6(pot),lc8(pot),lc10(pot), &
                                    lVex(pot),lbeta(pot), lkappa(pot))
             ydata(k,pot)=ydata(n,pot)*smooth(xSmooth) +               &
                          ydata(k,pot)*smooth(1d0-xSmooth)
             rr=rr+dr
          enddo
!
! DO 3) ->
!
!         define the short range ( r < xdata(n_sr+1) ) extrapolation
!         parameters
!                  V=slope*(r-fitpt_sr)+b
!
          fitpt_sr(pot)=xdata(n_sr+1,pot)
          slope(pot)=(ydata(n_sr+2,pot)-ydata(n_sr+1,pot))/           & 
                                (xdata(n_sr+2,pot)-xdata(n_sr+1,pot))
          b(pot)=ydata(n_sr+1,pot) 
!
          dr=0.1d0
          rr=xdata(n_sr+1,pot)-dr
          do k=n_sr,1,-1
             xdata(k,pot)=rr
             ydata(k,pot)=slope(pot)*(rr-fitpt_sr(pot))+b(pot)
             rr=rr-dr 
          enddo
!
          if(debug) then
          write(6,*)'# Potential that is actually used to interpolate'
            do k=1,ntot(pot)
              write(6,*) xdata(k,pot),                    &
                     ydata(k,pot)*HARTREE/xdata(k,pot)**6
            end do
          endif

!        enddo   ! pot loop    ! for driverPD

        if (debug) write(6,*) 'Syspot initialization complete'
!
! DONE with setup
!
        V=0d0
        SYS_SS_POT_KRb=V
!
        RETURN
      ENDIF
!
!=======================================================================
!
! --- SECTION FOR lsetup = .false.
!
      if(r.lt.fitpt_sr(NV)) then

       V=slope(NV)*(r-fitpt_sr(NV))+b(NV)

      else 
        dr=0.5d0
        if(r.LE.lfitpt_lr(NV)) THEN

          k=ntot(NV)
          xi(1) = r
          CALL akima(NP,k,xdata(1:k,NV),ydata(1:k,NV),1,xi,yi)
          V = yi(1)

        else

          V= longrangeXr6(r,lc6(nv),lc8(nv),lc10(nv),lVex(nv),  &
                             lbeta(nv),lkappa(nv))

        end if
      end if

      V=V/r**6        ! Note the r^6 factor 
      SYS_SS_POT_KRb=V

      RETURN
!
      END

!=======================================================================
!
!   RbRb singlet  (X) and triplet (a)  Born-Oppenheimer potentials  in
!   au
!
!          NV=1    singlet
!          NV=2    triplet      
!
!
!=======================================================================
!
!     V (r) = - ( c6/r^6 + c8/r^8 + c10/r^10 )  + Vex r^beta exp(-kappa*r)
!              for r to infinity   ( Note sign convention)
!
!=======================================================================
      DOUBLE PRECISION FUNCTION SYS_SS_POT_Rb2(R,V,NV,LSETUP)
      IMPLICIT none
      double precision :: r,V,dddDe,pi
      integer :: nv
      logical :: lsetup
      integer, parameter :: Num_ABOS=2
!
!     name of file containing table with potential values
!
      character(LEN=40), dimension(Num_ABOS) :: datafile
      character(LEN=40) :: filename
!
!     long-range dispersion coefficient, exchange coefficients,
!     and transition point  in a.u.
!
      double precision, dimension(Num_ABOS), save :: lc6,lc8,lc10,  & 
                                                     lVex,lbeta,lkappa
      double precision, dimension(Num_ABOS), save :: lfitpt_lr
      double precision :: c6,c8,c10, Vex,beta,kappa, fitpt_lr
      double precision, dimension(Num_ABOS), parameter ::    &
        c6def=(/ 4710.216338d0, 4710.216338d0 /),             &
        c8def=(/ 576696.4510d0, 576696.4510d0 /),             &
        c10def=(/ 75912808.97d0, 75912808.97d0 /),            &
        Vexdef=(/ -2.035435502016d-3, 2.035435502016d-3 /),   &
        betadef=(/ 5.317689d0, 5.317689d0 /),                 &
        kappadef=(/ 1.10799971d0,  1.10799971d0 /),           &
        fitpt_lrdef=(/25.0d0,25.0d0/)  
!
!!    short range control parameter in cm-1
!
      double precision, dimension(Num_ABOS), save :: ldDe,lRe,lRw
      double precision :: dDe, Re, Rw
      double precision, dimension(Num_ABOS), parameter ::   & 
         dDedef= (/ 0d0, 0d0/),     &  ! in cm-1
         Redef= (/ 12d0, 12d0/),    &  ! in a.u.
         Rwdef= (/  3d0, 3d0/)         ! in a.u.
!
!!    short range fit parameters and transition point 
!!                  slope*(x-fitpt_sr) + b  in a.u.
!
      double precision,dimension(Num_ABOS),  save :: fitpt_sr
      double precision,dimension(Num_ABOS),  save :: slope,b
!
!     arrays with data points
!
!        n_sr and n_lr added to short and long range to mimimize 
!        oscillation in the interpolated potential.
!
      integer, parameter :: NMAX=300
      double precision,save, dimension(NMAX,Num_ABOS) :: xdata,ydata
      integer, dimension(Num_ABOS), save :: ntot
      integer, parameter :: n_sr=3,n_lr=5
!
!     akima control parameter 
!           Np=number of polys in spline interval 
!
      integer, parameter :: NP=3
!
!     local variables and arrays
!
!       n:  number of of points read from file
!
      integer :: n
!
      integer :: k,pot
      logical :: ios
      double precision :: rr,dr,xSmooth
      double precision :: smooth
      double precision :: longrangeXr6
      double precision :: xi(1), yi(1)
      logical,save :: debug=.false.
!
!     constants
!
      double precision,parameter :: ALPHA=137.035999139D0,  &
                                    BOHRA=0.52917721067d0,  &
                                    HARTREE=219474.6313702d0
      common /c6/ c6
!
      namelist / XKK   / dDe,c6,c8,c10, Vex,beta,kappa,  & 
                           filename,fitpt_lr,debug,Re,Rw
      namelist / aKK   / dDe,c6,c8,c10, Vex,beta,kappa,  &
                           filename,fitpt_lr,debug,Re,Rw
     
       pi=dacos(-1d0)
!
!==========Initial setup================================================
!
              
      IF(lsetup) THEN           ! for driverPD
!     IF(lsetup.and.nv.eq.1) THEN    ! for the close coupling code

!        do pot=1,Num_ABOs    !  for cc
          pot=nv              !  for driverPD

          c6=C6def(pot)
          c8=C8def(pot)
          c10=C10def(pot)
          Vex=Vexdef(pot)
          beta=betadef(pot)
          kappa=kappadef(pot)
          dDe=dDedef(pot)
          Re=Redef(pot)
          Rw=Rwdef(pot)
          fitpt_lr=fitpt_lrdef(pot)
                              select case(pot)
                              case(1)
                                    filename='./KRb2_PESs/dimer/XRbRb.dat'
!            read(5,XKK)
                              case(2)
                                     filename='./KRb2_PESs/dimer/aRbRb.dat'      
!            read(5,aKK)
          end select
          lC6(pot)=C6
          lC8(pot)=C8
          lC10(pot)=C10
          lVex(pot)=Vex
          lbeta(pot)=beta
          lkappa(pot)=kappa
          ldDe(pot)=dDe
          lRe(pot)=Re
          lRw(pot)=Rw
          lfitpt_lr(pot)=fitpt_lr
!
!          write(6,*) '# parameters for pot= ',pot
!          write(6,*) '#        dDe= ',dDe
!          write(6,*) '#        Vex= ',Vex
!          write(6,*) '#       beta= ',beta
!          write(6,*) '#      kappa= ',kappa
!          write(6,*) '#   fitpt_lr= ',fitpt_lr
!          write(6,*) '#         C6= ',C6
!          write(6,*) '#         C8= ',C8
!          write(6,*) '#        C10= ',C10
!          write(6,*) '#         Re= ',Re
!          write(6,*) '#         Rw= ',Rw
!
          datafile(pot)=filename
!
! Now I need to 1) Read in the potential from CCSD(T) data 
!                  The potentials are in bohr and cm-1.        !!!!!!
!               2) Fit the inner ab-initio potential and 
!                  outer analytic form in a continuous manner.
!                  This fitting occurs at fitpt
!               3) Extrapolation parameters for the really short range 
!
!         read potentials for r< fitpt_lr, 
!         add shift dDe(cm-1), convert to atomic units, and 
!         multiply the potential by r**6.
!
          inquire(file=datafile(pot),EXIST=ios)
          if( .not.ios ) then
            write(6,*) 'Missing KRb2 dimer data file: ',trim(datafile(pot))
            write(6,*) 'Download the PCCP supplementary archive d6cp00282j1.zip'
            write(6,*) 'and place XKRb.dat/aKRb.dat/XRbRb.dat/aRbRb.dat under'
            write(6,*) './KRb2_PESs/dimer/ before running KRb2 calculations.'
            stop
          endif
          if (debug) then
            write(6,*) '----------------------------------'
            write(6,*) 'Initializing potential ',pot
            write(6,*) '----------------------------------'
          endif
          open(unit=20,status='old',file=datafile(pot))
          read(20,*) n
!
          k=1+n_sr
          read(20,*) xdata(k,pot),ydata(k,pot)
! simple shift
!          ydata(k,pot)=(ydata(k,pot)+ldDe(pot))*xdata(k,pot)**6/HARTREE
! arctan shift
!         dddDe=ldDe(pot)/pi*2d0
!    + *datan((lfitpt_lr(pot)-xdata(k,pot))
!    + /(lfitpt_lr(pot)-lRe(pot)))
! fermi dirac shift
            dddDe=dexp((xdata(k,pot)-lRe(pot))/lRw(pot))+1d0
            dddDe=ldDe(pot)/dddDe
          ydata(k,pot)=(ydata(k,pot)+dddDe)*xdata(k,pot)**6/HARTREE
!end pet          
          do while ( k.lt.n+n_sr .and. xdata(k,pot).lt.lfitpt_lr(pot) )
            k=k+1
            read(20,*) xdata(k,pot),ydata(k,pot)
! simple shift
!           ydata(k,pot)=(ydata(k,pot)+ldDe(pot))*xdata(k,pot)**6/
!    +                                                         HARTREE
! arctan shift
!         dddDe=ldDe(pot)/pi*2d0
!    + *datan((lfitpt_lr(pot)-xdata(k,pot))
!    + /(lfitpt_lr(pot)-lRe(pot)))
! fermi dirac shift     
          dddDe=dexp((xdata(k,pot)-lRe(pot))/lRw(pot))+1d0
            dddDe=ldDe(pot)/dddDe
            ydata(k,pot)=(ydata(k,pot)+dddDe)*xdata(k,pot)**6/HARTREE
!end pet
          end do
          close(20)
          n=k
          if ( xdata(k,pot).ge.lfitpt_lr(pot) ) n=n-1 

          ntot(pot)=n+n_lr
          if (n.le.2) then
          write(6,fmt="('Insufficient number of points in potential.' & 
                         ' n=',I2,' but must larger than 2')") n
            stop
          endif
          if ( ntot(pot).gt.NMAX) then
            write(6,fmt="('Insufficient space to store the ground ', &
              'state potential. Increase NMAX to ',I4 )") ntot(pot)
            stop
          endif
!
          if(debug) then
            write(6,*) 'now it does not work properly'
            write(6,*) '# original data from ', trim(datafile(pot))
            do k=1+n_sr,n
              write(6,*) xdata(k,pot),  &
                 ydata(k,pot)*HARTREE/xdata(k,pot)**6-ldDe(pot)
            end do
          endif
!
! DO 2) ->
!
!         replace or add to the ab initio potentials with the analytic
!         form for xdata > fitpt.  And make a smooth transition over 
!         the radial range xdata(n) and lfitpt_lr+(n_lr-1)*dr
!
          dr=0.5d0
          rr=lfitpt_lr(pot) 
!ML          rr=lfitpt_lr(pot)+2*dr
          do k=n+1,ntot(pot)
             xdata(k,pot)=rr
             xSmooth=(rr-xdata(n,pot))/(lfitpt_lr(pot)-xdata(n,pot))
!            xSmooth=(rr-xdata(n,pot))/(dr*(n_lr-1))
             ydata(k,pot)=longrangeXr6(rr,lc6(pot),lc8(pot),lc10(pot), &
                                     lVex(pot),lbeta(pot), lkappa(pot))
             ydata(k,pot)=ydata(n,pot)*smooth(xSmooth) +               &
                          ydata(k,pot)*smooth(1d0-xSmooth)
             rr=rr+dr
          enddo
!
!!DO 3) ->
!
!!        define the short range ( r < xdata(n_sr+1) ) extrapolation
!!        parameters
!!                 V=slope*(r-fitpt_sr)+b
!
          fitpt_sr(pot)=xdata(n_sr+1,pot)
          slope(pot)=(ydata(n_sr+2,pot)-ydata(n_sr+1,pot))/    & 
                      (xdata(n_sr+2,pot)-xdata(n_sr+1,pot))
          b(pot)=ydata(n_sr+1,pot) 
!
          dr=0.1d0
          rr=xdata(n_sr+1,pot)-dr
          do k=n_sr,1,-1
             xdata(k,pot)=rr
             ydata(k,pot)=slope(pot)*(rr-fitpt_sr(pot))+b(pot)
             rr=rr-dr 
          enddo
!
          if(debug) then
            write(6,*)'# Potential that is actually used to interpolate'
            do k=1,ntot(pot)
              write(6,*) xdata(k,pot),                      &
                     ydata(k,pot)*HARTREE/xdata(k,pot)**6
            end do
          endif

!        enddo   ! pot loop    ! for driverPD

        if (debug) write(6,*) 'Syspot initialization complete'
!
! DONE with setup
!
        V=0d0
        SYS_SS_POT_Rb2=V
!
        RETURN
      ENDIF
!
!=======================================================================
!
! --- SECTION FOR lsetup = .false.
!
      if(r.lt.fitpt_sr(NV)) then

       V=slope(NV)*(r-fitpt_sr(NV))+b(NV)

      else 
        dr=0.5d0
        if(r.LE.lfitpt_lr(NV)) THEN

          k=ntot(NV)
          xi(1) = r
          CALL akima(NP,k,xdata(1:k,NV),ydata(1:k,NV),1,xi,yi)
          V = yi(1)

        else

          V= longrangeXr6(r,lc6(nv),lc8(nv),lc10(nv),lVex(nv), &
                             lbeta(nv),lkappa(nv))

        end if
      end if

      V=V/r**6        ! Note the r^6 factor 
      SYS_SS_POT_Rb2=V

      RETURN
!
      END
!=======================================================================
      DOUBLE PRECISION FUNCTION longrangeXr6(r,c6,c8,c10,Vex,beta,kappa)
!
!       longrange= V_longrange*r^6   !!!!!!!!!!
!
      implicit none
      double precision, intent(IN) :: r,c6,c8,c10,Vex,beta,kappa
      double precision :: Ir2
      Ir2=1d0/(r*r)
      longrangeXr6=-(c6+c8*Ir2+c10*Ir2*Ir2)
      longrangeXr6=longrangeXr6 + Vex*r**(beta+6d0)*dexp(-kappa*r)
      end
!=======================================================================
      double precision function smooth(x)
      implicit none
!
!         This routine returns a monotonic function in x 
!         with the properties 
!                 \alpha(x)=1 for x<=0,  \alpha(x)=0 for x>=1, and
!                 \alpha'(x)=0 if x->0 and x->1.
!
      double precision :: x
      double precision, parameter ::  ONE=1d0,    &
                                     ZERO=0d0,    &
                                     HALF=0.5d0,  &
             PI= 3.14159265358979323846264330d0
!
      if ( x.gt.ONE ) then
        smooth=ZERO
      else
        if ( x.lt.ZERO ) then
          smooth=ONE
        else
          smooth=(COS( x*PI )+ONE)*HALF
        endif
      endif
      end
!================================================================================
!     End of dimer syspot functions!
!================================================================================


      MODULE Numeric_Kinds_Module
!USE Numbers_Module
!=========================================================================================
! This module defines types and constants to be used throughout other codes. Here is 
! defined the real type wp_kind (working precision) that defines the real type of all
! declared real variables. 
!
! D. A. Brue
!=========================================================================================
      IMPLICIT NONE

! Define compiler-dependent type constants as parameters
      INTEGER, PARAMETER :: i1_kind = SELECTED_INT_KIND(2)  ! -128, 127
      INTEGER, PARAMETER :: i2_kind = SELECTED_INT_KIND(4)  ! 
      INTEGER, PARAMETER :: i4_kind = SELECTED_INT_KIND(9)
      INTEGER, PARAMETER :: i8_kind = SELECTED_INT_KIND(12)
      
      INTEGER, PARAMETER :: sp_kind = SELECTED_REAL_KIND(p=6,r=30)    ! single precision
      INTEGER, PARAMETER :: dp_kind = SELECTED_REAL_KIND(p=14,r=300)  ! double precision
      INTEGER, PARAMETER :: qp_kind = SELECTED_REAL_KIND(p=30,r=4000) ! quadruple precision
! note above real kind types also work for complex declarations

! Select real kind to be used throughout the rest of code
      INTEGER, PARAMETER :: WP_Kind = dp_kind
      INTEGER, PARAMETER :: IW_Kind = i4_kind

      END MODULE Numeric_Kinds_Module
!=========================================================================

      SUBROUTINE AKIMA (NP, ND, XD, YD, NI, XI, YI)
      USE Numeric_Kinds_Module
      IMPLICIT NONE
!
! PREVIOUSLY THIS ROUTINE WAS CALLED UVIP3P AND WAS IN SINGLE PRECISION
! IT WAS DOWNLOADED FROM GAMS AT NIST ON OCT 21 1998
! IT PROVIDES SLIGHTLY BETTER FITTING THAN THE IMSL ROUTINE DCSAKM
! THIS IS A SUBJECTIVE COMMENT - IT DOES HOWEVER PREVENT
! OSCILLATIONS BETWEEN SPARCE DATA POINTS
!
! Univariate Interpolation (Improved Akima Method)
!
! Hiroshi Akima
! U.S. Department of Commerce, NTIA/ITS
! Version of 89/07/04
!
! This SUBROUTINE performs univariate interpolation.  It is based
! on the improved A method developed by Hiroshi Akima, 'A method
! of univariate interpolation that has the accuracy of a third-
! degree polynomial,' ACM TOMS, vol. xx, pp. xxx-xxx, 19xx.  (The
! equation numbers referred to in the comments below are those in
! the paper.)
!
! In this method, the interpolating function is a piecewise
! function composed of a set of polynomials applicable to
! successive intervals of the given data points.  This method
! uses third-degree polynomials as the default, but the user has
! an option to use higher-degree polynomial to reduce undulations
! in resulting curves.
!
! This method has the accuracy of a third-degree polynomial if
! the degree of the polynomials for the interpolating function is
! set to three.
!
! The input arguments are
!   NP = degree of the polynomials for the interpolating
!        function,
!   ND = number of input data points
!        (must be equal to 2 or greater),
!   XD = array of dimension ND, containing the abscissas of
!        the input data points
!        (must be in a monotonic increasing order),
!   YD = array of dimension ND, containing the ordinates of
!        the input data points,
!   NI = number of points for which interpolation is desired
!        (must be equal to 1 or greater),
!   XI = array of dimension NI, containing the abscissas of
!        the desired points.
!
! The output argument is
!   YI = array of dimension NI, where the ordinates of the
!        desired points are to be stored.
!
! If an INTEGER value smaller than 3 is given to the NP argument,
! this SUBROUTINE assumes NP = 3.
!
! The XI array elements need not be monotonic, but this
! SUBROUTINE interpolates faster if the XI array elements are
! given in a monotonic order.
!
! If the XI array element is less than XD(1) or greater than
! XD(ND), this SUBROUTINE linearly interpolates the YI value.
!
!
! Specification statement
      INTEGER :: ni, nd, idmn, idmx, idmd, iintpv, id, np
      INTEGER :: ID0, ID1, ID2, ID3, IEPT, II, IINT, IPE, NP0, NPM1
      REAL(KIND=DP_Kind) :: xd, yd, xi, yi, dx, dy, a0, a1, s2, a3, xii
      REAL(KIND=DP_Kind) :: x0,x1,x2,x3,y0,y1,y2,y3
      REAL(KIND=DP_Kind) :: A12, A13, A2, AA0, AA1, B0, B1, DLT, DNM,  &
                            DY0, DY1, DY2, DY3, EPSLN
      REAL(KIND=DP_Kind) :: PE, RENNM2, RENPM1, SMPEF, SMPEI, SMWTF,  &
                            SMWTI, SX, SXX, SXY, SY
      REAL(KIND=DP_Kind) :: T0, T1, U, UC, V, VOL, WT, XX, YP, YP0, YP1
      
      DIMENSION :: XD ( * ), YD ( * ), XI ( * ), YI ( * )
! Error check
      IF(ND.LE.1) GOTO 90
      IF(NI.LE.0) GOTO 91
      DO 11 ID = 2, ND
        IF(XD (ID) .LE.XD (ID-1) ) GOTO 92
   11 ENDDO
! Branches off special cases.
      IF(ND.LE.4) GOTO 50
! General case  --  Five data points of more
! Calculates some local variables.
      NP0 = MAX (3, NP)
      NPM1 = NP0 - 1
      RENPM1 = NPM1
      RENNM2 = NP0 * (NP0 - 2)
! Main calculation for the general case
! First (outermost) DO-loop with respect to the desired points
      DO 39 II = 1, NI
      IF(II.EQ.1) IINTPV = - 1
      XII = XI (II)
! Locates the interval that includes the desired point by binary
! search.
      IF(XII.LE.XD (1) ) THEN
        IINT = 0
      ELSEIF (XII.LT.XD (ND) ) THEN
        IDMN = 1
        IDMX = ND
        IDMD = (IDMN + IDMX) / 2
   31 IF(XII.GE.XD (IDMD) ) THEN
         IDMN = IDMD
      ELSE
         IDMX = IDMD
      ENDIF
      IDMD = (IDMN + IDMX) / 2
      IF(IDMD.GT.IDMN) GOTO 31
      IINT = IDMD
      ELSE
      IINT = ND
      ENDIF
! End of locating the interval of interest
! Interpolation or extrapolation in one of the three subcases
      IF(IINT.LE.0) THEN
! Subcase 1  --  Linear extrapolation when the abscissa of the
!                desired point is equal to that of the first data
!                point or less.
! Estimates the first derivative when the interval is not the
! same as the one for the previous desired point.  --
! cf. Equation (8)
      IF(IINT.NE.IINTPV) THEN
         IINTPV = IINT
         X0 = XD (1)
         X1 = XD (2) - X0
         X2 = XD (3) - X0
         X3 = XD (4) - X0
         Y0 = YD (1)
         Y1 = YD (2) - Y0
         Y2 = YD (3) - Y0
         Y3 = YD (4) - Y0
         DLT = X1 * X2 * X3 * (X2 - X1) * (X3 - X2) * (X3 - X1)
         A1 = ( ( (X2 * X3) **2) * (X3 - X2) * Y1 + ( (X3 * X1) ** &
          2) * (X1 - X3) * Y2 + ( (X1 * X2) **2) * (X2 - X1) &
          * Y3) / DLT
      ENDIF
! Evaluates the YI value.
      YI (II) = Y0 + A1 * (XII - X0)
! End of Subcase 1
      ELSEIF (IINT.GE.ND) THEN
! Subcase 2  --  Linear extrapolation when the abscissa of the
!                desired point is equal to that of the last data
!                point or greater.
! Estimates the first derivative when the interval is not the
! same as the one for the previous desired point.  --
! cf. Equation (8)
      IF(IINT.NE.IINTPV) THEN
         IINTPV = IINT
         X0 = XD (ND)
         X1 = XD (ND-1) - X0
         X2 = XD (ND-2) - X0
         X3 = XD (ND-3) - X0
         Y0 = YD (ND)
         Y1 = YD (ND-1) - Y0
         Y2 = YD (ND-2) - Y0
         Y3 = YD (ND-3) - Y0
         DLT = X1 * X2 * X3 * (X2 - X1) * (X3 - X2) * (X3 - X1)
         A1 = ( ( (X2 * X3) **2) * (X3 - X2) * Y1 + ( (X3 * X1) ** &
          2) * (X1 - X3) * Y2 + ( (X1 * X2) **2) * (X2 - X1) &
          * Y3) / DLT
      ENDIF
! Evaluates the YI value.
      YI (II) = Y0 + A1 * (XII - X0)
! End of Subcase 2
      ELSE
! Subcase 3  --  Interpolation when the abscissa of the desired
!                point is  between those of the first and last
!                data points.
! Calculates the coefficients of the third-degree polynomial (for
! NP.LE.3) or the factors for the higher-degree polynomials (for
! NP.GT.3), when the interval is not the same as the one for the
! previous desired point.
      IF(IINT.NE.IINTPV) THEN
         IINTPV = IINT
! The second DO-loop with respect to the two endpoints of the
! interval
         DO 37 IEPT = 1, 2
! Calculates the estimate of the first derivative at an endpoint.
! Initial setting for calculation
            ID0 = IINT + IEPT - 1
            X0 = XD (ID0)
            Y0 = YD (ID0)
            SMPEF = 0.0d0
            SMWTF = 0.0d0
            SMPEI = 0.0d0
            SMWTI = 0.0d0
! The third (innermost) DO-loop with respect to the four primary
! estimate of the first derivative
            DO 36 IPE = 1, 4
! Selects point numbers of four consecutive data points for
! calculating the primary estimate of the first derivative.
               IF(IPE.EQ.1) THEN
                  ID1 = ID0 - 3
                  ID2 = ID0 - 2
                  ID3 = ID0 - 1
               ELSEIF (IPE.EQ.2) THEN
                  ID1 = ID0 + 1
               ELSEIF (IPE.EQ.3) THEN
                  ID2 = ID0 + 2
               ELSE
                  ID3 = ID0 + 3
               ENDIF
! Checks if any point number falls outside the legitimate range
! (between 1 and ND).  Skips calculation of the primary estimate
! if any does.
      IF(ID1.LT.1.OR.ID2.LT.1.OR.ID3.LT.1.OR.ID1.GT.ND.OR.ID2.GT.ND &
                                             .OR.ID3.GT.ND) GOTO 36
! Calculates the primary estimate of the first derivative  --
! cf. Equation (8)
               X1 = XD (ID1) - X0
               X2 = XD (ID2) - X0
               X3 = XD (ID3) - X0
               Y1 = YD (ID1) - Y0
               Y2 = YD (ID2) - Y0
               Y3 = YD (ID3) - Y0
               DLT = X1 * X2 * X3 * (X2 - X1) * (X3 - X2) * &
                (X3 - X1)
               PE = ( ( (X2 * X3) **2) * (X3 - X2) * Y1 + ( (X3 * &
                X1) **2) * (X1 - X3) * Y2 + ( (X1 * X2) **2) &
                * (X2 - X1) * Y3) / DLT
! Calculates the volatility factor, VOL, and distance factor,
! SXX, for the primary estimate.  --  cf. Equations (9) and (11)
               SX = X1 + X2 + X3
               SY = Y1 + Y2 + Y3
               SXX = X1 * X1 + X2 * X2 + X3 * X3
               SXY = X1 * Y1 + X2 * Y2 + X3 * Y3
               DNM = 4.0d0 * SXX - SX * SX
               B0 = (SXX * SY - SX * SXY) / DNM
               B1 = (4.0d0 * SXY - SX * SY) / DNM
               DY0 = - B0
               DY1 = Y1 - (B0 + B1 * X1)
               DY2 = Y2 - (B0 + B1 * X2)
               DY3 = Y3 - (B0 + B1 * X3)
               VOL = DY0 * DY0 + DY1 * DY1 + DY2 * DY2 + DY3 * &
                DY3
! Calculates the EPSLN value, which is used to decide whether or
! not the volatility factor, VOL, is essentially zero.
               EPSLN = (YD (ID0) **2 + YD (ID1) **2 + YD (ID2) ** &
                2 + YD (ID3) **2) * 1.0D-12
! Accumulates the weighted primary estimates.  --
! cf. Equations (13) and (14)
               IF(VOL.GT.EPSLN) THEN
! - For finite weight.
                  WT = 1.0d0 / (VOL * SXX)
                  SMPEF = SMPEF + PE * WT
                  SMWTF = SMWTF + WT
               ELSE
! - For infinite weight.
                  SMPEI = SMPEI + PE
                  SMWTI = SMWTI + 1.0d0
               ENDIF
   36             ENDDO
! End of the third DO-loop
! Calculates the final estimate of the first derivative.  --
! cf. Equation (14)
            IF(SMWTI.LT.0.5d0) THEN
! - When no infinite weights exist.
               YP = SMPEF / SMWTF
            ELSE
! - When infinite weights exist.
               YP = SMPEI / SMWTI
            ENDIF
            IF(IEPT.EQ.1) THEN
               YP0 = YP
            ELSE
               YP1 = YP
            ENDIF
! End of the calculation of the estimate of the first derivative
! at an endpoint
   37          ENDDO
! End of the second DO-loop
         IF(NP0.LE.3) THEN
! Calculates the coefficients of the third-degree polynomial
! (when NP.LE.3).  --  cf. Equation (4)
            DX = XD (IINT + 1) - XD (IINT)
            DY = YD (IINT + 1) - YD (IINT)
            A0 = YD (IINT)
            A1 = YP0
            YP1 = YP1 - YP0
            YP0 = YP0 - DY / DX
            A2 = - (3.0d0 * YP0 + YP1) / DX
            A3 = (2.0d0 * YP0 + YP1) / (DX * DX)
         ELSE
! Calculates the factors for the higher-degree polynomials
! (when NP.GT.3).  --  cf. Equation (20)
            DX = XD (IINT + 1) - XD (IINT)
            DY = YD (IINT + 1) - YD (IINT)
            T0 = YP0 * DX - DY
            T1 = YP1 * DX - DY
            AA0 = (T0 + RENPM1 * T1) / RENNM2
            AA1 = - (RENPM1 * T0 + T1) / RENNM2
         ENDIF
      ENDIF
! End of the calculation of the coefficients of the third-degree
! polynomial (when NP.LE.3) or the factors for the higher-degree
! polynomials (when NP.GT.3), when the interval is not the same
! as the one for the previous desired point.
! Evaluates the YI value.
      IF(NP0.LE.3) THEN
! - With a third-degree polynomial.  --  cf. Equation (3)
         XX = XII - XD (IINT)
         YI (II) = A0 + XX * (A1 + XX * (A2 + XX * A3) )
      ELSE
! - With a higher-degree polynomial.  --  cf. Equation (19)
         U = (XII - XD (IINT) ) / DX
         UC = 1.0d0 - U
         V = AA0 * ( (U**NP0) - U) + AA1 * ( (UC**NP0) - UC)
         YI (II) = YD (IINT) + DY * U + V
      ENDIF
! End of Subcase 3
      ENDIF
   39 ENDDO
! End of the first DO-loop
! End of general case
      RETURN
! Special cases  --  Four data points or less
! Preliminary processing for the special cases
   50 X0 = XD (1)
      Y0 = YD (1)
      X1 = XD (2) - X0
      Y1 = YD (2) - Y0
      IF(ND.EQ.2) GOTO 60
      X2 = XD (3) - X0
      Y2 = YD (3) - Y0
      IF(ND.EQ.3) GOTO 70
      X3 = XD (4) - X0
      Y3 = YD (4) - Y0
      GOTO 80
! Special Case 1  --  Two data points
! (Linear interpolation and extrapolation)
   60 A1 = Y1 / X1
      DO 61 II = 1, NI
      YI (II) = Y0 + A1 * (XI (II) - X0)
   61 ENDDO
! End of Special Case 1
      RETURN
! Special Case 2  --  Three data points
! (Quadratic interpolation and linear extrapolation)
   70 DLT = X1 * X2 * (X2 - X1)
      A1 = (X2 * X2 * Y1 - X1 * X1 * Y2) / DLT
      A2 = (X1 * Y2 - X2 * Y1) / DLT
      A12 = 2.0d0 * A2 * X2 + A1
      DO 71 II = 1, NI
      XX = XI (II) - X0
      IF(XX.LE.0.0d0) THEN
      YI (II) = Y0 + A1 * XX
      ELSEIF (XX.LT.X2) THEN
      YI (II) = Y0 + XX * (A1 + XX * A2)
      ELSE
      YI (II) = Y0 + Y2 + A12 * (XX - X2)
      ENDIF
   71 ENDDO
! End of Special Case 2
      RETURN
! Special Case 3  --  Four data points
! (Cubic interpolation and linear extrapolation)
   80 DLT = X1 * X2 * X3 * (X2 - X1) * (X3 - X2) * (X3 - X1)
      A1 = ( ( (X2 * X3) **2) * (X3 - X2) * Y1 + ( (X3 * X1) **2) &
        * (X1 - X3) * Y2 + ( (X1 * X2) **2) * (X2 - X1) * Y3) / DLT
      A2 = (X2 * X3 * (X2 * X2 - X3 * X3) * Y1 + X3 * X1 * (X3 * X3 - &
          X1 * X1) * Y2 + X1 * X2 * (X1 * X1 - X2 * X2) * Y3) / DLT
      A3 = (X2 * X3 * (X3 - X2) * Y1 + X3 * X1 * (X1 - X3) * Y2 + X1 * &
            X2 * (X2 - X1) * Y3) / DLT
      A13 = (3.0d0 * A3 * X3 + 2.0d0 * A2) * X3 + A1
      DO 81 II = 1, NI
      XX = XI (II) - X0
      IF(XX.LE.0.0d0) THEN
      YI (II) = Y0 + A1 * XX
      ELSEIF (XX.LT.X3) THEN
      YI (II) = Y0 + XX * (A1 + XX * (A2 + XX * A3) )
      ELSE
      YI (II) = Y0 + Y3 + A13 * (XX - X3)
      ENDIF
   81 ENDDO
! End of Special Case 3
      RETURN
! Error exit
   90 WRITE ( *, 99090) ND
      GOTO 99
   91 WRITE ( *, 99091) NI
      GOTO 99
   92 WRITE ( *, 99092) ID, XD (ID-1), XD (ID)
   99 WRITE ( *, 99099)
      RETURN
! Format statements for error messages
99090 FORMAT(1X/ ' ***   Insufficient data points.', 7X,'ND =',I3)
99091 FORMAT(1X/ ' ***   No desired points.', 7X,'NI =',I3)
99092 FORMAT(1X/ ' ***   Two data points identical or out of ', &
            'sequence.'/ 7X,'ID, XD(ID-1), XD(ID) =',I5,2F10.3)
99099 FORMAT(' Error detected in the akima SUBROUTINE'/)
      END SUBROUTINE AKIMA
