SUBROUTINE btrns (nrho, ntheta, nchi, nchanl, lammin, lammax, nchir, nsize_old, psi_old, nsize_new, psi_new, tblock)
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
LOGICAL :: debug
INTEGER :: Indxr_New, nrho, ntheta, nchi, irho, itheta, ichi, lam, &
 lammin, nchanl, lammax, nsize_old, nsize_new, i, k, jchanl, nchir, kchir

REAL(Kind=WP_Kind) :: psi_old (nsize_old, nchanl), psi_new (nsize_new, nchanl), tblock (nchi, nchir)

debug = .true.
IF(debug)THEN
   WRITE(Out_Unit, * ) 'In btrns'
   WRITE(Out_Unit, * ) 'nrho=', nrho
   WRITE(Out_Unit, * ) 'ntheta=', ntheta
   WRITE(Out_Unit, * ) 'nchi=', nchi
   WRITE(Out_Unit, * ) 'nchanl=', nchanl
   WRITE(Out_Unit, * ) 'lammin=', lammin
   WRITE(Out_Unit, * ) 'lammax=', lammax
   WRITE(Out_Unit, * ) 'nchir=', nchir
   WRITE(Out_Unit, * ) 'nsize_old=', nsize_old
   WRITE(Out_Unit, * ) 'nsize_new=', nsize_new
   WRITE(Out_Unit, * ) 'lam,jchanl,ichi,itheta,irho,i,kchir,k'
   WRITE(Out_Unit, * ) 'tblock'
!cc         DO ichi=1,nchi
!cc            WRITE(Out_Unit,*)ichi,(tblock(ichi,kchir),kchir=1,nchir)
!cc         ENDDO
!cc         WRITE(Out_Unit,*)'psi_old'
!cc         DO kchir=1,nchir
!cc            WRITE(Out_Unit,*)kchir,(psi_old(kchir,jchanl),jchanl=1,nchanl)
!cc         ENDDO


ENDIF
DO lam = lammin, lammax
DO jchanl = 1, nchanl
DO ichi = 1, nchi
DO itheta = 1, ntheta
DO irho = 1, nrho
i = Indxr_New (lam, lammin, nrho, ntheta, nchi, irho, itheta, ichi)
psi_new (i, jchanl) = 0.d0
DO kchir = 1, nchir
k = Indxr_New (lam, lammin, nrho, ntheta, nchir, irho, itheta, kchir)
psi_new (i, jchanl) = psi_new (i, jchanl) + tblock (ichi, kchir) * psi_old (k, jchanl)
ENDDO
ENDDO
ENDDO
ENDDO
ENDDO

ENDDO
IF(debug)THEN
   WRITE(Out_Unit, * ) 'psi_new'
!cc         DO ichi=1,nchi
!cc            WRITE(Out_Unit,*)ichi,(psi_new(ichi,jchanl),jchanl=1,nchanl)
!cc         ENDDO

ENDIF
RETURN

ENDSUBROUTINE btrns
