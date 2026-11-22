SUBROUTINE dbsprop(nchanl, rmat, drmat, start, finish, nsteps, ksqmax, idname, method,frst_WRITE)
USE fileunits_Module
USE nrho_Module
USE boundary_Module
USE BasisRegions_Module
! METHOD1 -- Uses Matrix Inversions.
!
!
!             P U R P O S E    O F    S U B R O U T I N E
! Integrates second-order coupled differential equations across a region.
! written by G. A. Parker
!             I N P U T   A R G U M E N T S
!   nchanl    Number of coupled second-order differential equations.
!   rmat      LOG derivative matrix at the begining of the region.
!   drmat     Energy derivative of the R-matrix
!   w         Temporary storage array.
!   vecnew    Temporary storage array.
!   tstore    Temporary storage array.
!   nsteps    The number of steps in the log-derivative region.
!   start     Beginning distance for the log-derivative region.
!   finish    Ending distance for the log-derivative region.
!             O U T P U T     A R G U M E N T S
!   rmat    Wigner R-matrix at the end of the region.
!   drmat   Energy Derivative of the Wigner R-Matrix.
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
USE pmatdat_Module
IMPLICIT NONE
CHARACTER(LEN=10) method
CHARACTER(LEN=11) :: Blank
CHARACTER(LEN=21), PARAMETER:: ProcName='dbsprop'
CHARACTER(LEN=*) idname
LOGICAL frst_WRITE, calc_rmat
LOGICAL DEBUG,firstcall
LOGICAL WRITE_angbf
INTEGER nchanl, nsteps, isector, first_sector, last_sector, info
INTEGER isect,i
REAL(Kind=WP_Kind) rho_last,rho_now,rho_next
REAL(Kind=WP_Kind)  start, finish, sector_right, sector_left
REAL(Kind=WP_Kind)  ksqmax, left, right
INTEGER ipiv(nchanl)
REAL(Kind=WP_Kind) ui(nchanl,nchanl), work(nchanl)
REAL(Kind=WP_Kind), ALLOCATABLE :: omat(:,:,:),pm1b(:,:,:),pm2b(:,:,:)
REAL(Kind=WP_Kind), ALLOCATABLE :: e2db(:,:)
REAL(Kind=WP_Kind) w(nchanl,nchanl), rmat(nchanl,nchanl),vecnew(nchanl,nchanl)
REAL(Kind=WP_Kind) tstore(nchanl,nchanl), drmat(nchanl,nchanl)
data firstcall/.true./

SAVE

IF(firstcall)THEN
    CALL get_debug(procname,debug)
    IF(debug) WRITE(Out_Unit,*)sector_left, sector_right, nsteps
    firstcall = .false.
    WRITE(Out_Unit,*)'Sector Propagation with Matrix Inversions.'
    WRITE(msg_unit,*)'Sector Propagation with Matrix Inversions.'

    ALLOCATE(omat(NSectors,nchanl,nchanl))
    ALLOCATE(pm1b(NSectors,nchanl,nchanl))
    ALLOCATE(pm2b(NSectors,nchanl,nchanl))
    ALLOCATE(e2db(NSectors,nchanl))
    ALLOCATE (pm1(nchanl,nchanl),pm2(nchanl,nchanl))
    ALLOCATE(e2d(nchanl))
    DO i=1,NSectors
        CALL read_overlap(basis_dist(i), omat(i,:,:), nchanl)
        CALL read_potmat(basis_dist(i),pm1b(i,:,:),pm2b(i,:,:),e2db(i,:),nchanl)
    ENDDO

ENDIF
!
! Write out integration parameters.
!
WRITE(Out_Unit,'(A,I5,A,ES22.15,A,ES22.15)') ' nsteps=',nsteps, '   start=', start, '   finish=',finish
!
! Determine sector limits.
!
CALL bounds_of_sector(rho_last, rho_now, rho_next, sector_left, sector_right, start, finish,    &
                      isector, first_sector, last_sector)

WRITE(msg_unit,*)'method=',method
!  start loop over sectors
DO isector = first_sector, last_sector
    IF(right>=rho_boundary) WRITE_angbf = .true.
    CALL bounds_of_sector(rho_last, rho_now, rho_next, sector_left, sector_right, start, finish,    &
                          isector, first_sector, last_sector)
    !
    !   Interface for interpolation of potential matrix
    !
    pm1 = pm1b(isector,:,:)
    pm2 = pm2b(isector,:,:)
    e2d = e2db(isector,:)
    rhos(1) = sector_left
    rhos(2) = rho_now
    rhos(3) = sector_right
    !
    ! Calculate Overlaps of previous basis with the current basis.
    !
    IF(isector>first_sector)THEN
        vecnew = omat(isector,:,:)
        rmat=MatMul(Transpose(vecnew),MatMul(rmat,vecnew))
        IF(debug)THEN
            WRITE(Out_Unit,*)'logder_matrix at r=', sector_left
            WRITE(Out_Unit,*)'transformed with an overlp matrix calculated between functions at rho = ', &
                                rho_last,' and ',rho_now
            CALL MatrixOut(rmat, nchanl, nchanl,'LogDer','LogDer Matrix',Blank, Blank,.False., Blank, .False.)
        ENDIF
    ENDIF
    !
    !   integrate from sector_left to sector_right using the log-derivative integrator.
    !
    !  in each sector, propagate over the left half, then the right half
    !  with the following do loop.
    !
    IF(DEBUG)THEN
        WRITE(Out_Unit,*)' rho_last=',rho_last, ' rho_now  =',rho_now
        WRITE(Out_Unit,*)' sec_left=',sector_left, ' sec_right=',sector_right
        WRITE(Out_Unit,*)' start   =',start,       ' finish   =',finish
        WRITE(Out_Unit,*)' isector =',isector
        WRITE(Out_Unit,*)' fsector =',first_sector,' lsector  =',last_sector
        WRITE(Out_Unit,*)' nchanl  =',nchanl
    ENDIF
    DO isect=1,2
        IF(isect==1)THEN
            IF(isector==1.AND.sector_left>=rho_now) continue
            left=sector_left
            right=rho_now
            calc_rmat=.true.
        ELSE
            IF(isector==last_sector) exit
            calc_rmat=.false.
            left=rho_now
            right=sector_right
        ENDIF
        IF(debug) WRITE(Out_Unit,*)' left    =',left,        ' right    =',right
        !
        ! Propagate Coupled-Channel Equations.
        !
        IF(method=='vivsprop')THEN
            IF(isector==first_sector)rmat=0.0
            !IF(frst_WRITE)THEN
                !CALL vivs_WRITE(tstore, rmat, w, nchanl, left, right)
            !ELSE
            !    CALL vivs_READ(tstore,rmat,w,nchanl,left,right,nsteps,ksqmax,ui,work,ipiv)
            !ENDIF
        ELSEIF(method=='manolop ')THEN
            CALL manolop(rmat, nchanl, left,right, nsteps, ksqmax, method)
        ELSE
            CALL logder(rmat, drmat, nchanl, left, right, nsteps, ksqmax, method)
        ENDIF

        IF(debug)THEN
            IF(method=='vivsprop')THEN
                WRITE(Out_Unit,*) 'R_matrix at r=', right
                CALL MatrixOut(rmat, nchanl, nchanl,'VIVS R_Matrix=rmat', 'R_Matrix',Blank, Blank,.False., Blank, .False.)
            ELSE
                WRITE(Out_Unit,*) 'LogDer_Matrix at r=', right
                CALL MatrixOut(rmat, nchanl, nchanl,'LogDer=rmat', 'LogDer',Blank, Blank,.False., Blank, .False.)
                ui=rmat
                CALL SmxInv(ui,nchanl)
                CALL MatrixOut(ui, nchanl, nchanl,'R_Matrix=ui', 'R_Matrix',Blank, Blank,.False., Blank, .False.)
                IF(method=='dlogder ')THEN
                    WRITE(Out_Unit,*)'D[logder_matrix,E] at r=',right
                    CALL MatrixOut(drmat,nchanl,nchanl,'D[LogDer,E]', 'D[Logder_Matrix,E]', Blank, Blank,.False., Blank, .False.)
                ENDIF
            ENDIF
        ENDIF
        IF(calc_rmat)THEN
            ui=rmat
            IF(method/='vivsprop') CALL SmxInv(ui, nchanl)
            IF(debug)THEN
                WRITE(Out_Unit,*) 'R_matrix at r=', right
                CALL MatrixOut(ui, nchanl, nchanl,'R_Matrix','R_matrix',Blank, Blank,.False., Blank, .False.)
            ENDIF
!                CALL CIDBCond2(right, ui, nchanl)
            CALL flush(Out_Unit)
        ENDIF
222   continue
    ENDDO
ENDDO
!
! calculate R = Y**(-1)
!
333 Continue
IF(method=='logder  ')THEN
   CALL smxinv (rmat, nchanl)
ELSEIF(method=='manolop ')THEN
   CALL smxinv (rmat, nchanl)
ELSEIF(method=='dlogder ')THEN
   w=rmat
   CALL smxinv (rmat, nchanl)
   ui=Matmul(drmat,rmat)
   CALL dsysv ('l', nchanl, nchanl, w, nchanl, ipiv, ui,nchanl, work, nchanl, info)
   drmat=-ui
ENDIF

IF(debug)THEN
    WRITE(Out_Unit,*) 'R-Matrix at r=', right
    CALL MatrixOut(rmat,nchanl,nchanl,'R_Matrix','R_matrix',Blank, Blank,.False.,Blank,.False.)
    IF(method=='dlogder ')THEN
        WRITE(Out_Unit,*) 'D[R-Matrix,E] at r=', right
        CALL MatrixOut(drmat,nchanl,nchanl,'D[R_Matrix,E]','D[R_Matrix,E]',Blank, Blank,.False.,Blank,.False.)
    ENDIF
ENDIF

finish = sector_right
firstcall = .false.
RETURN
ENDSUBROUTINE dbsprop
