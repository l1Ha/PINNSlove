subroutine svdprop(Etot, nrchanl, nprop, start, finish, rn)
!-------------------------------------------------------------------
!
! package   : CID
!
! Language  : Fortran 90
!
! author    : F. Colavecchia  (flavioc@lanl.gov)
!             R. T Pack       (pack@lanl.gov)
!
! date      : Apr-2002       version: 0.0
! revision  :                version:
!
! purpose   :  Propagate the Delves basis with SVD-ERN method
!
! input     :    Etot -> Total energy
!              nrchanl -> number of functions in the basis
!               nprop -> number of functions actually propagating
!               start -> initial rho for the propagation
!              finish -> final rho for the propagation
!
! output    :    rn -> Ratio matrix of SVD-ERN
!
!
! modules   :     time -> Timing module
!                files -> files and units
!             boundary -> rho sectors
!             recovery -> recovery parameters
!                 nrho -> number of rho values
!               jacobi -> Jacobi BC
!               masses -> masses
! common    :
!
!
! notes     : Recovery code is obsolete.
!
!-------------------------------------------------------------------
USE Time_Module
use FileNames_Module
use boundary_Module
use recovery_Module
use nrho_Module
use jacobi_Module
use Masses_Module
use method_Module
USE BasisRegions_Module
!
!
IMPLICIT NONE

REAL(Kind=WP_Kind) Etot
CHARACTER(LEN=11) :: Blank
CHARACTER(LEN=21), parameter:: procname='svdprop'
logical                :: debug
logical                :: firstcall
logical                :: WRITE_angbf
integer                :: ienergy
integer                :: nenergy
integer                :: i, j, isector
integer                :: first_sector, last_sector
integer                :: nprop
integer, intent(in)    :: nrchanl
integer                :: info,lwork
integer,allocatable    :: ipiv(:)
REAL(Kind=WP_Kind), intent(in)     :: start, finish
REAL(Kind=WP_Kind)                 :: distance, delta,next_distance,prev_distance
REAL(Kind=WP_Kind)                 :: rhorec
REAL(Kind=WP_Kind)                 :: tnmin, tnmax
REAL(Kind=WP_Kind)                 :: sector_left, sector_right
REAL(Kind=WP_Kind) :: Secnds
REAL(Kind=WP_Kind), intent(out)    :: rn(nrchanl,nrchanl)
!REAL(Kind=WP_Kind), ALLOCATABLE    :: RN(:,:)
REAL(Kind=WP_Kind), ALLOCATABLE    :: Old2New(:,:)
REAL(Kind=WP_Kind), ALLOCATABLE    :: vecnew(:,:), vecold(:,:)
REAL(Kind=WP_Kind), ALLOCATABLE    :: w(:,:), work(:,:)
REAL(Kind=WP_Kind), ALLOCATABLE    :: Qn(:), Qnm1(:), tn(:), eig(:), tnm1(:), kh(:)
REAL(Kind=WP_Kind), ALLOCATABLE    :: eigen(:,:),omat(:,:,:),vmat(:,:,:)
!    REAL(Kind=WP_Kind), ALLOCATABLE    :: psn(:),psnm1(:)
!    REAL(Kind=WP_Kind), ALLOCATABLE    :: Old2Newm1(:,:),Old2Newm(:,:)
REAL(Kind=WP_Kind), ALLOCATABLE    :: RNm1(:,:)
REAL(Kind=WP_Kind), ALLOCATABLE    :: old2newp(:,:),workp(:,:),wp(:,:)
data firstcall /.true./
save

WRITE(Out_Unit,*) 'propagating with svd'
WRITE(msg_unit,*) 'propagating with svd'

lwork = 200*200
nprop = nrchanl
!ALLOCATE(RN(nprop,nprop))
ALLOCATE(Old2New(nrchanl,nrchanl) )
ALLOCATE(Old2Newp(nprop,nprop) )
ALLOCATE(vecnew(nrchanl,nprop), vecold(nrchanl,nprop) )
ALLOCATE(w(nrchanl,nrchanl), work(nrchanl,nrchanl) )
ALLOCATE(wp(nprop,nprop), workp(nprop,nprop) )
ALLOCATE(qn(nrchanl), tn(nrchanl), eig(nrchanl), tnm1(nrchanl) )
ALLOCATE(kh(nrchanl) )
ALLOCATE(ipiv(nrchanl) )
!
!
!
!ALLOCATE(pmat(NSectors,nrchanl,nrchanl))
ALLOCATE(vmat(NSectors,nrchanl,nrchanl))
ALLOCATE(omat(NSectors,nrchanl,nrchanl))
ALLOCATE(eigen(NSectors,nrchanl))
!
! First call
!
if(firstcall)THEN
   !
   ! get debugging options
   !
   call get_debug(procname, debug)
   firstcall = .false.
   !
   ! Determine sector information for all steps on this first call.
   ! note.  now, the mesh points used by the propagator are the centers of
   ! what used to be the sectors.  the sector boundaries are not used.
   ! SVD Propagator runs on a uniform grid.
   !
   if(start>=finish)THEN
      call errormsg(msg_unit, procname,' Not possible, start gt finish')
   ENDIF
   first_sector=1
   last_sector=NSectors
   do isector=2,NSectors
      if((basis_dist(isector-1)+basis_dist(isector))/2.d0<=start) first_sector=isector
      if((basis_dist(isector)+basis_dist(isector+1))/2.d0<finish) last_sector=isector
   ENDDO
   WRITE(Out_Unit,*)'first_sector=',first_sector
   WRITE(Out_Unit,*)' last_sector=',last_sector
   !
   !   On first call read all the matrices
   !
   call read_overlap(basis_dist(1), omat(1,1:nrchanl,1:nrchanl), nrchanl)
   call read_potmat(basis_dist(1),vmat(1,1:nrchanl,1:nrchanl),w,eig,nrchanl)

   if(debug) call MatrixOut(vmat(1,1:nrchanl,1:nrchanl),nrchanl,nrchanl,'w_read', 'potential',Blank, Blank,.false., Blank, .false.)
   do i=2,NSectors
      call read_overlap(basis_dist(i), omat(i,1:nrchanl,1:nrchanl), nrchanl)
      call read_potmat(basis_dist(i),vmat(i,1:nrchanl,1:nrchanl),w,eig,nrchanl)
      if(debug)THEN
            WRITE(Out_Unit,*) 'Matrices as read at distance',basis_dist(i)
            call MatrixOut(omat(i,1:nrchanl,1:nrchanl),nrchanl,nrchanl,'omat(i)', 'overlap',Blank, Blank,.false., Blank, .false.)
            call MatrixOut(vmat(i,1:nrchanl,1:nrchanl),nrchanl,nrchanl,'w_read', 'potential',Blank, Blank,.false., Blank, .false.)
!               WRITE(Out_Unit,*) eigen
      ENDIF
   ENDDO
   if (debug) call MatrixOut(eigen,NSectors,nrchanl,'eigen', 'eigenvalues',Blank, Blank,.false., Blank, .false.)
ENDIF
!
!  set step length
!
delta=(finish-start)/(last_sector-first_sector)
!prev_distance = start-delta
!
!if(recovery)THEN
!    call read_unformatted_matrix(done_dir//vecoldrec_name,rhorec,vecnew,nrchanl,nprop)
!    call read_unformatted_vector(done_dir//tnrec_name,rhorec,tn,nrchanl)
!    if(debug) WRITE(Msg_Unit,*) 'previous distance =',prev_distance, 'rhorec=',rhorec
!    if(abs(prev_distance-rhorec)>1d-4)THEN
!        WRITE(msg_unit,*) 'prev_distance  =',prev_distance
!        WRITE(msg_unit,*) 'rho_rec        =',rhorec
!        call errormsg(msg_unit,procname,'Recovery of Tn and VecOld failed prev_distance/=rhorec')
!    ENDIF
!else
   call unit(nprop,vecnew)
!ENDIF !
!  Get the method for boundary conditions
!  we will WRITE the angular functions if distance>= rho_boundary
!
WRITE_angbf = .false.      ! Delves BC
!WRITE_angbf = .true.      ! Jacobi BC
!
!
!  ******** propagate solutions from start to finish. *********
!  ******** start do loop on sectors (steps) ***********
!
do isector = first_sector, last_sector
   !
   ! Move all the info of the previous sector to the done directory
   ! At this point, working has the information of isector-1
   ! and            done    has the information of isector-2, if exists.
      WRITE(Out_Unit,*) '---------------------------------------------'
      WRITE(Out_Unit,*) 'SECTOR =',isector
      WRITE(Out_Unit,*) 'rho    =',basis_dist(isector)
      WRITE(Out_Unit,*) '---------------------------------------------'

   !call update
   !
   ! Timing
   !
   SecndsElapsed = SecndsNow
   SecndsNow=0.d0
   SecndsElapsed = SecndsNow-SecndsElapsed
   if (debug) WRITE(msg_unit,*) 'Time =',SecndsNow
   if (debug) WRITE(msg_unit,*) 'TimeE=',SecndsElapsed
   if (debug) WRITE(msg_unit,*) 'TimeF=',timeframe
   !
   ! Now, we STOP the code if we are near the time frame
   !
   if(abs(SecndsNow+SecndsElapsed)>timeframe) STOP
   !
   ! get sector information and basis functions.
   !
   distance=basis_dist(isector)
   if(distance>=rho_boundary.AND.jacobi) WRITE_angbf = .true.

   if(debug) WRITE(Out_Unit,*) 'distance      =',distance, 'distance+delta=',distance+delta, 'delta         =',delta
   !
   ! Read overlap matrix.
   !
   !if(isector>1.or.recovery)then
   !call read_overlap(distance, Old2New, nrchanl)
   Old2new = omat(isector,1:nrchanl,1:nrchanl)
   !ENDIF
   !
   ! Save previous transformation matrix
   !
   vecold=vecnew
   !
   ! Read the transformation matrix
   !
   !
   ! The diagonalization of the potential matrix (-w) is now done
   ! in the calculation of the Basis. There we change the sign of w
   ! and diagonalize, now we just read their eigenvalues and eigenvectors.
   !
   ! diagonalize the interaction matrix.
   !
   vecnew = vmat(isector,1:nrchanl,1:nprop)
   !call djacobi(w,nrchanl,nrchanl,Eig,work,info)
   !do i=1,nrchanl
   !    do j=1,nprop
   !        vecnew(i,j) = work(i,j)
   !    ENDDO
   !ENDDO
   !call deigsrt(Eig,VecNew,nrchanl,nrchanl)
   !
   ! w is now the interaction matrix chosen to have diagonal parts
   ! positive in repulsive (classically forbidden) region.
   !
   ! The transformation matrix we have is energy INDEPENDENT,
   ! so we need to add now the energy dependent term to the
   ! eigenvalues of the transformation matrix.
   !
   do i=1,nrchanl
      eig(i) = eig(i) - Etot*usys2
   ENDDO

   if(debug)then
      !WRITE(Out_Unit,*)'potential at distance  =',distance
      !call MatrixOut(w,nrchanl,nrchanl,'potential', 'potential',Blank, Blank,.false., Blank, .false.)
      WRITE(Out_Unit,*)'eigenvalues at distance=', distance
      WRITE(Out_Unit,'(5es17.8e3)') eig
      WRITE(Out_Unit,*)'eigenvectors at distance', distance
      call MatrixOut(vecnew, nrchanl, nprop,'vectors', 'eigenvectors',Blank, Blank,.false., Blank, .false.)
   ENDIF
   ! save previous tn then determine new one
   Tnm1=Tn
   Tn=Eig*(delta**2)/12.d0
   if(debug) WRITE(Out_Unit,*) 'Tn   =',Tn
   if(debug) WRITE(Out_Unit,*) 'Tnm1 =',Tnm1
   ! usual renormalized numerov.
   ! qn=(2.d0+10.d0*tn)/(1.d0-tn)
   ! new renormalized numerov to allow bigger stepsizes.
   do i=1,nrchanl
      kh(i) = delta*sqrt(abs(eig(i)))
      !  classically forbidden region is eig>0.
      if (eig(i)>=0.d0)then
            Qn(i) = exp(kh(i)) + exp(-kh(i))
      else
            Qn(i) = 2.d0*cos(kh(i))
      ENDIF
   ENDDO
   !
   !   end of new renormalized numerov determination of qn
   !   check step size
   !   see notes of rtp of 00/9/25
   !
   tnmin = minval(Tn)
   tnmax = maxval(Tn)
   ! debug
   if(debug)then
      if(tnmin<=-0.09d0.or.tnmax>=0.2d0)then
            WRITE(msg_unit,*)'tnmax= ', tnmax
            WRITE(msg_unit,*)'tnmin= ', tnmin
            call warning(msg_unit,procname, 'stepsize in renormalized numerov, tnmin, tnmax')
      ENDIF
      WRITE(Out_Unit,*)'tn'
      WRITE(Out_Unit,'(5es17.8e3)') tn
      WRITE(tn_Unit,'(26e18.6)') distance,(tn(i),i=1,25)
      WRITE(Out_Unit,*)'qn'
      WRITE(Out_Unit,'(5es17.8e3)') qn
   ENDIF
   !
   ! propagate solutions using svd algorithm
   ! this assumes that for isector=0, the wavefunction is zero.
   ! that is, it assumes all channels start out classically forbidden.
   ! rn is the renormalized numerov ratio matrix.
   !
   if(isector==1)THEN
      if(recovery)THEN
            call read_unformatted_matrix(done_dir//rnrec_name,rhorec,rn,nprop,nprop)
            if(debug) WRITE(Msg_Unit,*) 'previous distance =',prev_distance, 'rhorec=',rhorec
            if(abs(prev_distance-rhorec)>1d-4)THEN
               STOP 'svd: Recovery of Rn failed prev_distance/=rhorec'
            ENDIF
      else
            RN=0.d0
      !call unit(nrchanl,Old2New)
            call unit(nprop,Old2Newp)
            do i=1,nprop
               !  RN(i,i)=1.d+30
               !  usual renormalized numerov start
               !  RN(i,i) = qn(i)
               !  new renormalized numerov start
               RN(i,i) = exp(kh(i))
            ENDDO
      ENDIF
   ENDIF
   !
   ! Now, we need to propagate if isector>1 and it's a first run
   !                        or if we recover from previous calculations, for ALL isectors.
   !
   if(isector>1.or.recovery)THEN
      !
      ! transform Old2New to the fully adiabatic bases
      !
      call MatrixOut (Old2New,nrchanl,nrchanl,'Old2New_before', 'o2nb',Blank, Blank,.false., Blank, .false.)
      if(coord=='del')THEN
            Old2Newp=matmul(transpose(vecold),matmul(Old2New,vecnew))
      else
            Old2Newp = Old2New(1:nprop,1:nprop)
      ENDIF
      call MatrixOut (vecold,nrchanl,nprop,'vecold_matrix', 'vold',Blank, Blank,.false., Blank, .false.)
      call MatrixOut (Old2Newp,nprop,nprop,'Old2New_matrix', 'o2n',Blank, Blank,.false., Blank, .false.)
      !
      ! ********************start boundary conditions*********************
      ! we now have the necessary info to apply boundary conditions.
      ! here is where the bc must be applied before calculating
      ! the new rn matrix.  we are using r(n-1)u(n-1)=o(n-1,n)u(n).
      ! see notes of rtp of 00/10/20 to 00/11/08
      !
      if(distance>=rho_boundary)THEN
      call ratio_matrix(Old2Newp,tnm1,tn,vecold,vecnew,rn,nrchanl,nprop,w)
      !
      !    Output the ratio matrix if needed
      !

            !call cidbcond3(w, nrchanl, distance,delta,tn,tnm1)
      ENDIF
      ! *****************end boundary conditions ******************
      !
      ! go on to get new rn
      !
      call smxinv(RN,nprop)
      RN=matmul(transpose(Old2Newp),matmul(RN,Old2Newp))
      !  rn is now o(n,n-1)*r(n-1)**(-1)*o(n-1,n); i.e., the transformed
      !  inverse of the r at isector-1=n-1
      !  now get the new ratio matrix at isector, i.e., r(n).
      RN = -RN
      do i=1,nprop
         RN(i,i)=qn(i)+RN(i,i)
      ENDDO
   ENDIF
   !
   ! Save the data of the propagation
   !
   call save_prop(distance, RN, VecNew, Tn, nrchanl, nprop)
   !
   ! Save the data for recovery in *.rec
   ! In the case the run stops after getting the last_sector,
   ! the recovery run will STOP.
   next_distance = start+isector*delta
   !call WRITE_main_rec(mainrec_unit,next_distance,isector+1)
   call WRITE_main_rec(next_distance)
   !if(debug) call WRITE_main_rec(msg_unit,next_distance,isector+1)
   if(debug)then
      !
      ! WRITE overlap matrix o if desired.
      !
      call MatrixOut (Old2New, nrchanl, nrchanl,'Old2New', 'Old2New', Blank, Blank,.false., Blank, .false.)
      !
      ! WRITE renormalized numerov matrix
      !
      WRITE(Out_Unit,*)'renormalized_numerov_matrix at distance=', distance,'  for isector=',isector
      call MatrixOut (RN, nprop, nprop,'rn', 'renormalized numerov matrix',Blank, Blank,.false., Blank, .false.)
   ENDIF
   call flush(Out_Unit)
ENDDO
!
!  ********  end of do loop on sectors (steps)  ********************
!
!      RN=w             ! store ratio matrix in rn.
!
! DEALLOCATE storage
!
DEALLOCATE(Old2New )
DEALLOCATE(vecnew, vecold )
DEALLOCATE(w, work )
DEALLOCATE(qn, tn, eig, kh )
!
!   If the current energy is the last one, we can DEALLOCATE
!   some matrices.
!
if (ienergy==nenergy)THEN
   DEALLOCATE(vmat,omat,eigen)
ENDIF

!
! return to logprop
!
return
end
