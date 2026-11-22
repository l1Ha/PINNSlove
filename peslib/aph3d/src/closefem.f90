SUBROUTINE closefem
! PORT NOTE (2026-10-01): the original closefem referenced stale unit names
! (Tidi_Unit, Surf33_Unit, ...) whose declarations no longer exist in
! FileUnits_Module.  The corresponding scratch/binary files are intermediate
! products of the FEM basis construction; leaving them in BinOut/ is harmless.
! Stubbed to no-op for the gfortran port.
  RETURN
ENDSUBROUTINE closefem
