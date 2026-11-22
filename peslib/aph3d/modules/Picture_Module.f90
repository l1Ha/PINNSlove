MODULE Picture_Module
   USE Numeric_Kinds_Module
   USE Dipole_Module, ONLY: NSurf
   IMPLICIT NONE
   SAVE
   LOGICAL logpot                   !
   LOGICAL irdsurf                  !
   LOGICAL denplot                  !
   INTEGER(KIND=IW_Kind) kolb(2)       !
   INTEGER(KIND=IW_Kind) nthmode(100)  !
   INTEGER(KIND=IW_Kind) nl            !
   INTEGER(KIND=IW_Kind) nil           !
   INTEGER(KIND=IW_Kind):: iproj=2       !
   INTEGER(KIND=IW_Kind) klrplt        !
   INTEGER(KIND=IW_Kind) kgloss        !
   INTEGER(KIND=IW_Kind) kase          !
   INTEGER(KIND=IW_Kind) npltskp       !
   INTEGER(KIND=IW_Kind) iplt          !
   REAL(Kind=WP_Kind) cutoffl       !
   REAL(Kind=WP_Kind) cutoffh       !
   REAL(Kind=WP_Kind) radii         !
   REAL(Kind=WP_Kind) chij          !
   REAL(Kind=WP_Kind) range         !
   REAL(Kind=WP_Kind) view(3)       !
   REAL(Kind=WP_Kind) xsn(4)        !
   REAL(Kind=WP_Kind) cent(3)       !
   REAL(Kind=WP_Kind) win(4)        !
   REAL(Kind=WP_Kind) shad(5)       !
   REAL(Kind=WP_Kind) hue(11)       !
   REAL(Kind=WP_Kind) sat(11)       !
   REAL(Kind=WP_Kind) cint(11)      !
END MODULE Picture_Module
