!=====================================================================================
module mod_pot_HeStar_H2 
  implicit none 
  integer,parameter:: rk = selected_real_kind(15,307) 
  real(rk),parameter:: pi = 3.14159265358979323846_rk 
  private 
  public :: pot_HeStar_H2,cart2dist,jacob2cart,pes_gp,pot_HeStar_H2_rrr
contains
function pot_HeStar_H2_rrr(x,y,z) result(v)
  real(rk),intent(in) :: x,y,z 
  !x: R_{He--H2} (a.u.)
  !y: r_{H--H} (a.u.)
  !z: theta_{center(H-H)--He} (rad)
  !v: potential energy (a.u.)
  real(rk) :: v 

  real(rk):: carts(3,3),dist(3),cov(1),ys(1)

 ! call jacob2cart(x,y,z,carts)
!  call cart2dist(3,carts,dist)
dist(1) = x 
dist(2) = y 
dist(3) = z 
  call pes_gp(3,1,dist,ys,cov,1)  
  v = ys(1) 
end function pot_HeStar_H2_rrr

function pot_HeStar_H2(x,y,z) result(v)
  real(rk),intent(in) :: x,y,z 
  !x: R_{He--H2} (a.u.)
  !y: r_{H--H} (a.u.)
  !z: theta_{center(H-H)--He} (rad)
  !v: potential energy (a.u.)
  real(rk) :: v 

  real(rk):: carts(3,3),dist(3),cov(1),ys(1)

  call jacob2cart(x,y,z,carts)
  call cart2dist(3,carts,dist)
  call pes_gp(3,1,dist,ys,cov,1)  
  v = ys(1) 
end function pot_HeStar_H2 
!=========================================================================
subroutine cart2dist(natoms,x,r)
  implicit none
  integer, intent(in) :: natoms
  real*8, intent(in) :: x(3,natoms)
  real*8, intent(out) :: r(natoms*(natoms-1)/2)
  real*8 :: dx(3)
  integer :: i,j,k
  k=0
  do i=1,natoms-1
    do j=i+1,natoms
      k=k+1
      dx=x(:,i)-x(:,j)
      r(k)=sqrt(dot_product(dx,dx))
    end do
  end do
  return
end
!=========================================================================
subroutine jacob2cart(Z,r,theta,x)
  implicit none
  real*8, intent(in) :: Z,r,theta
  real*8, intent(out) :: x(3,3)
  x=0.d0
  x(1,1)=Z*cos(theta)
  x(2,1)=Z*sin(theta)
  x(1,2)=-0.5d0*r
  x(1,3)= 0.5d0*r
  return
end
!=========================================================================
subroutine pes_gp(n0,ns,xs0,ys,cov,icov)
  implicit none
  ! --- variables used for a PES using GP model --------------------------------
  integer,save :: ndim,np,hypsize,d
  real*8,allocatable,save :: x(:,:),hyp(:)
  real*8,save :: sn2
  character(kind=1,len=20),save :: covfunc,xhandle
  real*8,save :: mcoef
  real*8,allocatable,save :: alpha(:),y(:),Cinv(:,:),alpha0(:)
  integer,save :: init=0
  character(kind=1,len=200),parameter :: finput="GPinput.txt"
  integer :: info
  ! ----------------------------------------------------------------------------
  integer :: i,j,fid
  ! --- variables for unknown geometries, input xs, output ys -------------
  integer,intent(in) :: n0,ns,icov
  real*8,intent(in)  :: xs0(n0,ns)
  real*8 :: xs(n0,ns),norm(1)
  real*8,intent(out) :: ys(ns),cov(ns)
  real*8,allocatable :: Ks(:,:)

  ! --- loading variables for gaussian process at first call -------------------
  if(init.eq.0) then
    write(*,"('GP PES initialization:')")
    fid=101010
    open(fid,file=trim(finput),status='old',action='read')
    read(fid,*) ndim,np
    read(fid,*) covfunc
    if(trim(covfunc).eq.'covMaternard'.or.trim(covfunc).eq.'covSEard') then
       hypsize=ndim+1
    else
       hypsize=2
    end if
   
    allocate(hyp(hypsize))
    read(fid,*) d
    do i=1,hypsize
      read(fid,*) hyp(i)
    end do
    read(fid,*) sn2
    read(fid,*) xhandle
    read(fid,*) mcoef

    allocate(x(ndim,np),alpha(np),y(np),Cinv(np,np),alpha0(np))

    do i=1,np
        read(fid,*) x(1:ndim,i),alpha0(i),y(i)
    enddo

    close(fid)

    if(trim(covfunc).eq.'covMaternard') then
      call covMaternard(d,hyp,ndim,x,np,x,np,0,Cinv)
    elseif(trim(covfunc).eq.'covMaterniso') then
      call covMaterniso(d,hyp,ndim,x,np,x,np,0,Cinv)
    endif

    Cinv=Cinv/sn2
    forall(i=1:np) Cinv(i,i)=Cinv(i,i)+1.d0
    call dpotrf('U',np,Cinv,np,info)
    if(info.ne.0) stop 'Cholesky factorization failed in Cinv!'
    call dpotri('U',np,Cinv,np,info)
    if(info.ne.0) stop 'Matrix inverse failed!'
    !fill up
    Cinv=Cinv/sn2
    do i=2,np
      do j=1,i-1
            Cinv(i,j)=Cinv(j,i)
      end do
    end do

    call dsymv('U',np,1.d0,Cinv,np,y,1,0.d0,alpha,1)

    y=alpha0-alpha
    norm=sqrt(dot_product(y,y))

    write(*,"('Dimension:',1x,i5)") ndim
    write(*,"('No. of points:',1x,i8)") np
    write(*,*) covfunc
    write(*,"('Norm:',1x,E20.8)") norm

    init=1
  endif
  ! --- end of loading variables for gaussian process --------------------------

  ! --- do prediction using Gaussian Process -----------------------------------
  if(n0.ne.ndim) then
    print*,"input xs0 dimension wrong"
    stop
  endif

  !sort dist
  xs=xs0
  do i=1,ns
    if(xs0(1,i).lt.xs0(2,i)) then
      xs(1,i)=xs0(2,i)
      xs(2,i)=xs0(1,i)
    end if
  end do

  if(trim(xhandle).eq.'normal') then
      !do nothing
  elseif(trim(xhandle).eq.'reverse') then
      xs=1.d0/xs
  elseif(trim(xhandle).eq.'morse') then
      xs=exp(-mcoef*xs)
  endif

  !!calculate potential
  allocate(Ks(ns,np))

  if(trim(covfunc).eq.'covMaternard') then
    call covMaternard(d,hyp,ndim,xs,ns,x,np,0,Ks)
  elseif(trim(covfunc).eq.'covMaterniso') then
    call covMaterniso(d,hyp,ndim,xs,ns,x,np,0,Ks)
  endif

  call dgemv('N',ns,np,1.d0,Ks,ns,alpha,1,0.d0,ys,1)
  ys=ys/27.2114d0

  if(icov.gt.0) then
      do i=1,ns
            if(trim(covfunc).eq.'covMaternard') then
                  call covMaternard(d,hyp,ndim,xs(:,i),1,xs(:,i),1,0,norm)
            elseif(trim(covfunc).eq.'covMaterniso') then
                  call covMaterniso(d,hyp,ndim,xs(:,i),1,xs(:,i),1,0,norm)
            endif
            call dsymv('U',np,1.d0,Cinv,np,Ks(i,1:np),1,0.d0,y,1)
            cov(i)=norm(1)-dot_product(Ks(i,1:np),y(1:np))
      end do
  end if

  deallocate(Ks)
  return
end
!======================================================================================
subroutine sq_dist(D,a,n,b,m,C)
  ! sq_dist - a function to compute a matrix of all pairwise squared distances
  ! between two sets of vectors, stored in the columns of the two matrices, a
  ! (of size D by n) and b (of size D by m). If only a single argument is given
  ! or the second matrix is empty, the missing matrix is taken to be identical
  ! to the first.
  ! a is of size Dxn, b is of size Dxm, C is of size nxm.
  !
  ! Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2010-12-13.
  ! Fortran by Jun Chen (Email: NJUChenJun@gmail.com), 2015-11-26
  implicit none
  integer,intent(in) :: D,n,m
  real*8,intent(in)  :: a(D,n),b(D,m)
  real*8,intent(out) :: C(n,m)
  real*8 :: ax(D,n),bx(D,m),amean(D),bmean(D),mu(D),asum(n),bsum(m)
  integer :: i,j

  ! subtract mean
  do i=1,D
      amean(i)=sum(a(i,1:n))/(n*1.d0)
      bmean(i)=sum(b(i,1:m))/(m*1.d0)
  enddo
  mu = (1.d0*m/(n+m))*amean + (1.d0*n/(n+m))*bmean
  do i=1,n
      ax(1:D,i)=a(1:D,i)-mu(1:D)
  enddo
  do i=1,m
      bx(1:D,i)=b(1:D,i)-mu(1:D)
  enddo

  ! compute squared distances
  asum=0.d0
  do j=1,n
  do i=1,D
      asum(j)=asum(j)+ax(i,j)**2
  enddo
  enddo
  bsum=0.d0
  do j=1,m
  do i=1,D
      bsum(j)=bsum(j)+bx(i,j)**2
  enddo
  enddo

  call dgemm('t','n',n,m,D,-2.d0,ax,D,bx,D,0.d0,C,n)
  do j=1,m
  do i=1,n
      C(i,j)=C(i,j)+asum(i)+bsum(j)
  enddo
  enddo

  ! numerical noise can cause C to negative i.e. C > -1e-14
  C = max(C,0.d0)
  return
end subroutine sq_dist
!===================================================================================
subroutine covMaterniso(d, hyp, ndim, x, nx, z, nz, iderivatives, K)
  ! Matern covariance function with nu = d/2 and isotropic distance measure. For
  ! d=1 the function is also known as the exponential covariance function or the 
  ! Ornstein-Uhlenbeck covariance in 1d. The covariance function is:
  !
  !   k(x^p,x^q) = sf^2 * f( sqrt(d)*r ) * exp(-sqrt(d)*r)
  !
  ! with f(t)=1 for d=1, f(t)=1+t for d=3 and f(t)=1+t+t²/3 for d=5.
  ! Here r is the distance sqrt((x^p-x^q)'*inv(P)*(x^p-x^q)), P is ell times
  ! the unit matrix and sf2 is the signal variance. The hyperparameters are:
  !
  !
  ! Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2010-09-10.
  ! Fortran by Jun Chen (Email: NJUChenJun@gmail.com), 2015-11-26
  implicit none
  integer,intent(in) :: ndim,nx,nz,d,iderivatives
  real*8,intent(in)  :: hyp(2),x(ndim,nx),z(ndim,nz)
  real*8,intent(out) :: K(nx,nz)
  real*8 :: ell,sf2,kt(nx,nz)
  integer :: i,j 
  ell = dsqrt(hyp(1))
  sf2 = hyp(2)

  if(d.ne.1 .and. d.ne.3 .and. d.ne.5) then
    print*,'only 1, 3 and 5 allowed for d'
    stop
  endif
  if(iderivatives.ne.0 .and. iderivatives.ne.1 .and. iderivatives.ne.2) then
    print*,'Unknown hyperparameter'
    print*,'only 0, 1 and 2 allowed for iderivatives'
    stop
  endif

  call sq_dist(ndim,dsqrt(d*1.d0)/ell*x,nx,dsqrt(d*1.d0)/ell*z,nz,K)
  do i=1,nx
  do j=1,nz
      K(i,j)=dsqrt(K(i,j))
  enddo
  enddo

  if (iderivatives.eq.0) then     ! covariances
    call f_covMaterniso(K,kt,nx,nz,d,0)
    do i=1,nx
    do j=1,nz
      K(i,j) = sf2*kt(i,j)*exp(-K(i,j))
    enddo
    enddo
  elseif (iderivatives.eq.1) then ! derivatives
    call f_covMaterniso(K,kt,nx,nz,d,1)
    do i=1,nx
    do j=1,nz
      K(i,j) = sf2*kt(i,j)*exp(-K(i,j))*K(i,j)
    enddo
    enddo
  elseif (iderivatives.eq.2) then
    call f_covMaterniso(K,kt,nx,nz,d,0)
    do i=1,nx
    do j=1,nz
      K(i,j) = 2.d0*sf2*kt(i,j)*exp(-K(i,j))
    enddo
    enddo
  else
    print*,'Unknown hyperparameter'
    stop
  endif
  return
end subroutine covMaterniso
!==================================================================================
subroutine f_covMaterniso(t,f,na,nb,d,k)
  ! calculate f(t) or df(t), df(t)=f(t)-f'(t)
  ! if d=1: f = 1; df = 1;
  ! if d=3: f = 1+t; df = t;
  ! if d=5: f = 1 + t.*(1+t/3); df = t.*(1+t)/3;
  ! Fortran by Jun Chen (Email: NJUChenJun@gmail.com), 2015-11-26
  implicit none
  integer,intent(in) :: na,nb,d,k
  real*8,intent(in)  :: t(na,nb)
  real*8,intent(out) :: f(na,nb)
  integer :: i,j
  if(d.ne.1 .and. d.ne.3 .and. d.ne.5) then
    print*,'only 1, 3 and 5 allowed for d'
    stop
  endif
  if(k.ne.0 .and. k.ne.1) then
    print*,'only 0 and 1 allowed for k, to calculate f or df'
    stop
  endif
  if (d.eq.1) then
    f=1.d0
  elseif (d.eq.3) then
    f=t
    if(k.eq.0) f=f+1.d0
  elseif (d.eq.5) then
    if(k.eq.0) then
      do i=1,na
      do j=1,nb
        f(i,j)=1.d0+t(i,j)*(1.d0+t(i,j)/3.d0)
      enddo
      enddo
    elseif(k.eq.1) then
      do i=1,na
      do j=1,nb
        f(i,j)=t(i,j)*(1.d0+t(i,j))/3.d0
      enddo
      enddo
    endif
  endif
  return
end subroutine f_covMaterniso
!=====================================================================================
subroutine covMaternard(d, hyp, ndim, x, nx, z, nz, iderivatives, K)
  ! Matern covariance function with nu = d/2 and with Automatic Relevance
  ! Determination (ARD) distance measure. For d=1 the function is also known as
  ! the exponential covariance function or the Ornstein-Uhlenbeck covariance 
  ! in 1d. The covariance function is:
  !
  !   k(x^p,x^q) = sf^2 * f( sqrt(d)*r ) * exp(-sqrt(d)*r)
  !
  ! with f(t)=1 for d=1, f(t)=1+t for d=3 and f(t)=1+t+t²/3 for d=5.
  ! Here r is the distance sqrt((x^p-x^q)'*inv(P)*(x^p-x^q)), where the P matrix
  ! is diagonal with ARD parameters ell_1^2,...,ell_D^2, where D is the dimension
  ! of the input space and sf2 is the signal variance. The hyperparameters are:
  !
  ! Copyright (c) by Hannes Nickisch, 2013-10-13.
  ! Fortran by Jun Chen (Email: NJUChenJun@gmail.com), 2015-11-27
  implicit none
  integer,intent(in) :: ndim,nx,nz,d,iderivatives
  real*8,intent(in)  :: hyp(ndim+1),x(ndim,nx),z(ndim,nz)
  real*8,intent(out) :: K(nx,nz)
  real*8 :: ell(ndim),sf2,kt(nx,nz),xt(ndim,nx),zt(ndim,nz),ki(nx,nz)
  integer :: i,j

  ell = sqrt(hyp(1:ndim))
  ! print*,size(ell),size(hyp)
  ! do i=1,ndim 
  !   ell(i) = sqrt(hyp(i))
  ! enddo 

  sf2 = hyp(ndim+1)

  if(d.ne.1 .and. d.ne.3 .and. d.ne.5) then
    print*,'only 1, 3 and 5 allowed for d'
    stop
  endif
  if(iderivatives.lt.0 .or. iderivatives.gt.ndim+1) then
    print*,'only 0 to ndim+1 allowed for iderivatives'
    stop
  endif
  ! precompute distances
  do j=1,nx
  do i=1,ndim
    xt(i,j)=dsqrt(d*1.d0)/ell(i)*x(i,j)
  enddo
  enddo

  do j=1,nz
  do i=1,ndim
    zt(i,j)=dsqrt(d*1.d0)/ell(i)*z(i,j)
  enddo
  enddo

  call sq_dist(ndim,xt,nx,zt,nz,K)
  do i=1,nx
  do j=1,nz
      K(i,j)=dsqrt(K(i,j))
  enddo
  enddo

  !m = @(t,f) f(t).*exp(-t); dm = @(t,f) df(t).*exp(-t);
  if (iderivatives.eq.0) then ! covariances
    call f_covMaternard(K,kt,nx,nz,d,0)
    do i=1,nx
    do j=1,nz
      K(i,j)=sf2*kt(i,j)*exp(-K(i,j))
    enddo
    enddo
  else                        ! derivatives
    if (iderivatives.le.ndim) then              ! length scale parameter
      call sq_dist(1,dsqrt(d*1.d0)/ell(iderivatives)*x(iderivatives,1:nx),& 
      nx,dsqrt(d*1.d0)/ell(iderivatives)*z(iderivatives,1:nz),nz,ki);
      call f_covMaternard(K,kt,nx,nz,d,1)
      do i=1,nx
      do j=1,nz
        K(i,j) = sf2*kt(i,j)*exp(-K(i,j))*Ki(i,j)
        if(Ki(i,j).le.1.d-12) K(i,j)=0.d0
      enddo
      enddo
    elseif (iderivatives.eq.(ndim+1)) then      ! magnitude parameter
      call f_covMaternard(K,kt,nx,nz,d,0)
      do i=1,nx
      do j=1,nz
        K(i,j)=2.d0*sf2*kt(i,j)*exp(-K(i,j))
      enddo
      enddo
    else
      print*,'Unknown hyperparameter'
      stop
    endif
  endif

  return
end subroutine covMaternard
!=======================================================================================
subroutine f_covMaternard(t,f,na,nb,d,k)
  ! calculate f(t) or df(t), df(t)=(f(t)-f'(t))/t
  ! case 1, f = @(t) 1;               df = @(t) 1./t;
  ! case 3, f = @(t) 1 + t;           df = @(t) 1;
  ! case 5, f = @(t) 1 + t.*(1+t/3);  df = @(t) (1+t)/3;
  ! Fortran by Jun Chen (Email: NJUChenJun@gmail.com), 2015-11-27
  implicit none
  integer,intent(in) :: na,nb,d,k
  real*8,intent(in)  :: t(na,nb)
  real*8,intent(out) :: f(na,nb)
  integer :: i,j
  if(d.ne.1 .and. d.ne.3 .and. d.ne.5) then
    print*,'only 1, 3 and 5 allowed for d'
    stop
  endif
  if(k.ne.0 .and. k.ne.1) then
    print*,'only 0 and 1 allowed for k, to calculate f or df'
    stop
  endif

  if (d.eq.1 .and. k.eq.0) then
    f=1.d0
  elseif (d.eq.3 .and. k.eq.0) then
    f=t+1.d0
  elseif (d.eq.5 .and. k.eq.0) then
    do i=1,na
    do j=1,nb
      f(i,j)=1.d0+t(i,j)*(1.d0+t(i,j)/3.d0)
    enddo
    enddo
  elseif (d.eq.1 .and. k.eq.1) then
    do i=1,na
    do j=1,nb
      f(i,j)=1.d0/t(i,j)
    enddo
    enddo
  elseif (d.eq.3 .and. k.eq.1) then
    f=1.d0
  elseif (d.eq.5 .and. k.eq.1) then
    do i=1,na
    do j=1,nb
      f(i,j)=(1.d0+t(i,j))/3.d0
    enddo
    enddo
  endif
  return
end subroutine f_covMaternard
!========================================================================================
end module mod_pot_HeStar_H2 

! subroutine potface1(ipes,icoord,mass3,capr,smlr,theta,vtot)
!   !!!!input
!   !theta in radian
!   !R in a.u.
!   !r in a.u.
!   !!!output
!   !energy in hartree 
!   !zero is set to the energy of E(He)+E(H)+E(H+)
!   use amod_pot_HeStar_H2,only:pot_HeStar_H2
!   implicit none
!   real*8, intent (in) :: capr,smlr,theta,mass3(3)
!   integer,intent(in):: ipes,icoord
!   real*8, intent (out) :: vtot
!   ! print*,'NO.',ipes,'PES','Icoord',icoord,'Mass',mass3
!   vtot = pot_HeStar_H2(capr,smlr,theta)
!   return
!   ! print*,'Pot gave coord',capr,smlr,theta,vtot 
!   end subroutine potface1
!   subroutine potfacem(ipes,icoord,mass3,capr,smlr,theta,vtot)
!     !!!!input
!     !theta in radian
!     !R in a.u.
!     !r in a.u.
!     !!!output
!     !energy in hartree 
!     !zero is set to the energy of E(He)+E(H)+E(H+)
!     use amod_pot_HeStar_H2,only:pot_HeStar_H2
!     implicit none
!     real*8, intent (in) :: capr,smlr,theta,mass3(3)
!     integer,intent(in):: ipes,icoord
!     real*8, intent (out) :: vtot
!     ! print*,'NO.',ipes,'PES','Icoord',icoord,'Mass',mass3
!     vtot = pot_HeStar_H2(capr,smlr,theta)
!     return
!     ! print*,'Pot gave coord',capr,smlr,theta,vtot 
!     end subroutine potfacem
!   subroutine pes_init
!     implicit none

!     print*,'HE vs H2+ PES '

!     return
! end subroutine
