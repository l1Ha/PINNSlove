!*==FORMELEM.f90  processed by SPAG 6.55Dc at 15:05 on 21 Aug 2007
SUBROUTINE FormElem(BeginVec,EndVec,NTheta,NChi,valmax)
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
!*--FORMELEM6
 
! P U R P O S E    O F   S U B R O U T I N E
! This routine WRITEs coordinates of the Elements and the connectivitiy
! array to disk for subsequent plotting.
 
! I N P U T   A R G U M E N T S
!  BeginVec Coordinates for the start of each vector.
!  EndVec   Coordinates for the termination of each vector.
!  NTheta      Number of APH Theta coordinates.
!  NChi        Number of APH Chi coordinates.
INTEGER ntheta, nchi, i, ichi, itheta
REAL(KIND=WP_Kind) valmax
REAL(KIND=WP_Kind) BeginVec(3,NTheta-1,NChi-1), EndVec(3,NTheta-1,NChi-1)
REAL(KIND=WP_Kind) CoordMax(6)
DATA CoordMax/ - 1.0, 1.0, -1.0, 1.0, -1.0, 1.0/

OPEN(Unit=e_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'PlotVee/FormElem.txt')
WRITE(e_unit,99001)NChi*(NTheta-1) + 1, (NTheta-1)*(NChi-1)             !      0
 
99001 FORMAT(1x,3I10)
WRITE(e_unit,99004)(valmax*CoordMax(i),i=1,6)
WRITE(e_unit,99004)BeginVec(1,1,1), BeginVec(2,1,1), BeginVec(3,1,1)
 
DO IChi = 1, NChi - 1
   DO ITheta = 1, NTheta - 1                                            !      0
       WRITE(e_unit,99004)EndVec(1,Itheta,IChi), EndVec(2,Itheta,IChi),         &
                        & EndVec(3,ITheta,IChi)
   ENDDO
ENDDO
 
DO ITheta = 1, NTheta - 1                                               !      0
   WRITE(e_unit,99004)EndVec(1,ITheta,1), EndVEc(2,ITheta,1), EndVec(3,ITheta,1)
ENDDO
 
DO IChi = 1, NChi - 1                                                   !      0
   DO ITheta = 1, NTheta - 1                                            !      0
       IF(ITheta==1)THEN                                                !      0
           WRITE(e_unit,99002)ITheta, 2 + (IChi-1)*(NTheta-1),                  &
                            & 2 + IChi*(NTheta-1)                       !      0
99002      FORMAT(1x,'3  ',3I10)
       ELSE
           WRITE(e_unit,99003)ITheta + (IChi-1)*(NTheta-1),                     &
                            & ITheta + 1 + (IChi-1)*(NTheta-1),                 &
                            & ITheta + 1 + IChi*(NTheta-1),                     &
                            & ITheta + IChi*(NTheta-1)
99003      FORMAT(1x,'4  ',4I10)
       ENDIF
   ENDDO
ENDDO
99004 FORMAT(1x,6F10.5)
 
ENDSUBROUTINE FORMELEM
