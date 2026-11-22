SUBROUTINE funct2 (r, s, h, p, xj, det, xx, iintp)
USE FileUnits_Module
USE ToDim_Module
USE Var_Module
USE Numbers_Module

!            P U R P O S E   O F    S U B R O U T I N E
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   p r o g r a m                                                   .
! .                                                                   .
! .      to find interpolation functions ( h )                        .
! .      and derivatives ( p ) corresponding to the nodal points      .
! .      of a 4- to 8-node isoparametric quadrilateral                .
! .                                                                   .
! .      to find jacobian ( xj ) and its determinant ( det )          .
!
! .                                                                   .
! .                    node numbering convention                      .
!
! .                                                                   .
! .               2                 5                 1               .
! .                                                                   .
! .                 o . . . . . . . o . . . . . . . o                 .
! .                 .                               .                 .
! .                 .                               .                 .
! .                 .               s               .                 .
! .                 .               .               .                 .
! .                 .               .               .                 .
! .               6 o             9 o . . r         o 8               .
! .                 .                               .                 .
! .                 .                               .                 .
! .                 .                               .                 .
! .                 .                               .                 .
! .                 .                               .                 .
! .                 o . . . . . . . o . . . . . . . o                 .
! .                                                                   .
! .               3                 7                 4               .
! .                                                                   .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .


IMPLICIT NONE
INTEGER iintp, iperm, i, nn, ih, in, i2, i1, j, k, jend
REAL(Kind=WP_Kind) r, s, h, p, xj, det, xx, rp, rm, spvar, r2, s2, sm, dum
DIMENSION h(*), p(2,*), iperm(4), xj(2,2), xx(2,*)
DATA iperm /2,3,4,1/
!
rp=one+r
spvar=one+s
rm=one-r
sm=one-s
r2=one-r*r
s2=one-s*s
!
!  interpolation functions and their derivatives
!
!  4-node element
!
h(1)=quarter*rp*spvar
h(2)=quarter*rm*spvar
h(3)=quarter*rm*sm
h(4)=quarter*rp*sm
p(1,1)=quarter*spvar
p(1,2)=-p(1,1)
p(1,3)=-quarter*sm
p(1,4)=-p(1,3)
p(2,1)=quarter*rp
p(2,2)=quarter*rm
p(2,3)=-p(2,2)
p(2,4)=-p(2,1)
i=0
 10   i=i+1
IF(i > 5) GOTO 70
nn=i
IF(nn==1)THEN
!
   h(5)=half*r2*spvar
   p(1,5)=-r*spvar
   p(2,5)=half*r2
   GOTO 10
ELSEIF(nn==2)THEN
   h(6)=half*rm*s2
   p(1,6)=-half*s2
   p(2,6)=-rm*s
   GOTO 10
ELSEIF(nn==3)THEN
   h(7)=half*r2*sm
   p(1,7)=-r*sm
   p(2,7)=-half*r2
   GOTO 10
ELSEIF(nn==4)THEN
   h(8)=half*rp*s2
   p(1,8)=half*s2
   p(2,8)=-rp*s
   GOTO 10
ELSEIF(nn==5)THEN
   h(9)=r2*s2
   p(1,9)=-two*r*s2
   p(2,9)=-two*s*r2
   GOTO 10
ENDIF
!
!  correct functions and derivatives IF 5 or more nodes are
!  used to describe the element
!
 70   ih=0
 80   ih=ih+1
IF(ih > 5) GOTO 130
in=ih+4
IF(in == 9) GOTO 100
i1=in-4
i2=iperm(i1)
h(i1)=h(i1)-half*h(in)
h(i2)=h(i2)-half*h(in)
h(ih+4)=h(in)
DO j=1, 2
   p(j,i1)=p(j,i1)-half*p(j,in)
   p(j,i2)=p(j,i2)-half*p(j,in)
   p(j,ih+4)=p(j,in)
ENDDO
GOTO 80
100  DO j=1, 4
   h(j)=h(j)+quarter*h(9)
   p(1,j)=p(1,j)+quarter*p(1,9)
   p(2,j)=p(2,j)+quarter*p(2,9)
ENDDO
jend=4
DO j=1, jend
   in=j+4
   h(in)=h(in)-half*h(9)
   p(1,in)=p(1,in)-half*p(1,9)
   p(2,in)=p(2,in)-half*p(2,9)
ENDDO
GOTO 80
!
!  evaluate the jacobian matrix at point (r,s)
!
 130  IF(iintp > 0) RETURN
DO i=1, 2
   DO j=1, 2
      dum=zero
      DO k=1, 9
         dum=dum+p(i,k)*xx(j,k)
      ENDDO
      xj(i,j)=dum
   ENDDO
ENDDO
!
!  compute the determinant of the jacobian matrix at point (r,s)
!
det=xj(1,1)*xj(2,2)-xj(2,1)*xj(1,2)
IF(det<.00000001d0)THEN
   WRITE(Out_Unit,*) 'error: execution stopped in routine funct2'
   WRITE(Out_Unit,*) 'zero jacobian determinant for element'
   STOP 'funct2'
ENDIF
!
RETURN
ENDSUBROUTINE funct2
