    SUBROUTINE Bisect(r, v, rleft, vleft, rright, vright, delta)
    USE Numeric_Kinds_Module
USE Numbers_Module
    USE FileUnits_OneDim_Module
    IMPLICIT NONE
    REAL(KIND=WP_Kind) r, v, rleft, vleft, rright, vright, delta
    LOGICAL, PARAMETER:: idbug=.false.

    IF(idbug)THEN
       WRITE(Dbug_Unit,*)
       WRITE(Dbug_Unit,*)"On entering Bisect:"
       WRITE(Dbug_Unit,'(2(a,es23.15))')"     r=",r,     "      v=",v
       WRITE(Dbug_Unit,'(2(a,es23.15))')" rleft=",rleft, "  vleft=",vleft
       WRITE(Dbug_Unit,'(2(a,es23.15))')"rright=",rright," vright=",vright
    ENDIF
    IF(vleft*v.lt.Zero)THEN
       rright=r
       vright=v
    ELSE IF(vleft*v.gt.Zero)THEN
       rleft=r
       vleft=v
    ELSE
      rleft=r
      vleft=v
      rright=r
      vright=v
    ENDIF
    r=(rleft+rright)/Two
    delta=ABS(rright-rleft)/r

    IF(idbug)THEN
       WRITE(Dbug_Unit,*)"Bisect results:"
       WRITE(Dbug_Unit,'(3(a,es23.15))')"     r=",r,     "      v=",v," delta=",delta
       WRITE(Dbug_Unit,'(2(a,es23.15))')" rleft=",rleft, "  vleft=",vleft
       WRITE(Dbug_Unit,'(2(a,es23.15))')"rright=",rright," vright=",vright
    ENDIF

    RETURN
    END SUBROUTINE Bisect
