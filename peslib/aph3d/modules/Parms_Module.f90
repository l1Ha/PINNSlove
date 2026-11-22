MODULE Parms_Module
   USE Narran_Module
   USE DVR_Module
   IMPLICIT NONE
   SAVE
   INTEGER MxAng                       ! INTEGER, PARAMETER :: MxAng=  1920
   INTEGER MxBasis                     ! INTEGER, PARAMETER :: MxBasis= 215
   INTEGER MxBasiss                    !
   INTEGER MxGlegn                     !
   INTEGER MxHermt                     !
   INTEGER Mxl                         !
   INTEGER MxMega                      ! INTEGER, PARAMETER :: MxMega=    0
   INTEGER Mxn                         ! INTEGER, PARAMETER :: Mxn=       5
   INTEGER NVibRot                     !
   INTEGER NVbrthrt                    !
   INTEGER MaxOsc                      !
   INTEGER maxhermt                    !
   INTEGER mxnode
   INTEGER mxelmnt
   INTEGER, PARAMETER :: npoint=300    !
   INTEGER, PARAMETER :: NRhoSec=   3  !
   INTEGER, PARAMETER :: MaxDim=1000   ! PORT 2026-10-04: raised 600->1000 for lam_n=700
   INTEGER, PARAMETER :: MxVib=40      ! PORT: raised 12->40 for large ABM bases
END MODULE Parms_Module
