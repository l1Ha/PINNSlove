PROGRAM Program_KRb2
!-----------------------------------------------------------------------
! PORT (2026-10-02): driver for the K + Rb2 / Rb + KRb scattering
! calculation with the KRb2 2x2 ESI potential (surface.f90 dispatcher
! selects surface_krb2 because PES_Name contains 'KRb2').
!-----------------------------------------------------------------------
USE FileUnits_Module
CHARACTER(LEN=100), PARAMETER:: SysDir='APH3D_KRb2/'

CALL APH3D_InitRoots

OutDir=TRIM(OutputRoot)//TRIM(SysDir)
InputDir=TRIM(InputRoot)//TRIM(SysDir)

CALL BasisForScattering

STOP 'Program_KRb2'
END PROGRAM Program_KRb2
