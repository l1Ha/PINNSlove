MODULE PhiLis_Module
   ! PORT NOTE (2026-10-01): duplicates of Mshdat_Module variables; re-export
   ! to avoid ambiguous use-association in consumers.
   USE Mshdat_Module, ONLY: tolfn1, tolfn2, xkmesh
   USE Numeric_Kinds_Module
   IMPLICIT NONE
   SAVE
   REAL(Kind=WP_Kind) tolfn   ! local (no Mshdat counterpart)
END MODULE PhiLis_Module
