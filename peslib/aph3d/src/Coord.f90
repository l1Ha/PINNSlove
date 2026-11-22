SUBROUTINE Coord(x,y,z,Theta,Chi,r)
USE FileUnits_Module
USE Numeric_Kinds_Module
IMPLICIT NONE
REAL(KIND=WP_Kind) x, y, z, Theta, Chi, r
 
x = r*sin(Theta)*cos(Chi)
y = r*sin(Theta)*sin(Chi)
z = r*cos(Theta)
 
ENDSUBROUTINE COORD
