MODULE QDiv_Module
   USE Numeric_Kinds_Module
   USE Numbers_Module
   IMPLICIT NONE
   SAVE
   ! PORT NOTE (2026-10-01, gfortran port): mxdiv raised 7 -> 24.  The shipped
   ! input decks (e.g. APH3D_H3_DMBE/DMBE_ABM.data) supply ~21 entries for
   ! nqth/nqch/thdiv/chdiv, which overflowed the 2017 arrays and made the
   ! qdivp namelist read fail ("Cannot match namelist object name 0").
   INTEGER, PARAMETER    :: mxdiv=  24   !
   INTEGER               :: nthdiv=  2                   !
   INTEGER(KIND=IW_Kind) :: nchdiv=  7        !
   INTEGER(KIND=IW_Kind) :: nqth(mxdiv)  = [20,20,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]   !
   INTEGER(KIND=IW_Kind) :: nqch(mxdiv)  = [13,10,20,10,20,10,13, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] !
   REAL(Kind=WP_Kind)    :: thdiv(mxdiv) = [0.D0, 0.75D0, 1.5707963267949D0, 0.D0, 0.d0, 0.d0, 0.d0, &
                                            0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,   &
                                            0.d0,0.d0,0.d0,0.d0,0.d0,0.d0] !
   REAL(Kind=WP_Kind)    :: chdiv(mxdiv+1) = [-HalfPi, -1.2472D0, -0.8472D0,  -0.2D0, 0.2D0, 0.8472D0, 1.2472D0, HalfPi, &
                                              0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,                        &
                                              0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0] !
   NAMELIST/qdivp/ nthdiv, nqth, thdiv, nchdiv, nqch, chdiv
END MODULE QDiv_Module
