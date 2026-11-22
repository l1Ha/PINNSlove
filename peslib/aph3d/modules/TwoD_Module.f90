MODULE TwoD_Module
   USE Numeric_Kinds_Module
   USE Parms_Module
   IMPLICIT NONE
   SAVE
   REAL(Kind=WP_Kind), ALLOCATABLE:: csurf1(:,:)  ! csurf1(ndvptmax,nsfnmax)
   REAL(Kind=WP_Kind), ALLOCATABLE:: csurf2(:,:)  ! csurf2(ndvptmax,nsfnmax)
   REAL(Kind=WP_Kind), ALLOCATABLE:: ovlp(:,:)    ! ovlp(nsfnmax,nsfnmax)
! NEW FROM Xuan Li
   INTEGER nstate
   REAL(Kind=WP_Kind), ALLOCATABLE:: e2d_new(:,:),e3d_pro(:,:)
   REAL(Kind=WP_Kind), ALLOCATABLE:: cc2dtot(:,:,:),cc2dtot1(:,:,:),cc2dtot2(:,:,:),vectortemp(:,:),cc2dtemp(:,:)
   REAL(Kind=WP_Kind), ALLOCATABLE:: cc2d(:,:,:),cc3dtot(:,:,:)
   REAL(Kind=WP_Kind), ALLOCATABLE:: e3d(:)
   REAL(Kind=WP_Kind), ALLOCATABLE:: Genergy(:)
   INTEGER      , ALLOCATABLE:: Gindex(:,:),Gindex2(:)
   REAL(Kind=WP_Kind), ALLOCATABLE:: Gbasis(:,:),Gbasistemp(:,:),v_pro(:,:)
   INTEGER      , ALLOCATABLE:: nc3tot(:)
   INTEGER      , ALLOCATABLE:: nev_3d(:),nev_3d1(:),nev_3d2(:)
END MODULE TwoD_Module
