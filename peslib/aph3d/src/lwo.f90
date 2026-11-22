SUBROUTINE lwo (npnts, vertice, nsrfs, polslist, flname)
USE FileUnits_Module
USE Numeric_Kinds_Module
!-------------------------------------------------------------------------------
!
!                        LIGHTWAVE OBJECT (lwo_)
!
! This program creates a binary file of an object that can be viewed through
! Lightwave 3D.  By default, the object has a gray color*, is double sided**,
! is smoothed**, and has 100% diffusion**.
!
!   npnts == total number of vertices in the object
!
!   vertice[] == array of all the x, y, and z coordinates of the vertices
!                in the object (ex.  {0,1,2,3,4,5}  correspondes to the
!                zeroth vertice with the coordinate (0,1,2) and the first
!                vertice with the coordinate (3,4,5))
!
!   nsrfs == total number of individual polygons that make up the object
!
!   polslist[] == connectivity array of the polygons.  Each polygon is
!                 represented by the number of vertices in that polygon
!                 followed by that many integers refferring to the vertices
!                 from the vertice[] array.  It is only necessary to list each
!                 once for a given polygon, but each vertice must be listed in
!                 clockwise order with respect to the side of the plygon that
!                 will be seen (ex. in {3,0,2,1,4,2,3,1,4}, the first polygon
!                 has three vertices: the zeroth vertice, the second vertice,
!                 and the first vertice, connected in that order (Lightwave
!                 automatically contects the first vertice to the zeroth).
!                 The second polygon has four vertices, etc.)
!
!   flname[] == filename string (".lwo" is added automatically.
!               "keith.sub" is changed to "keith.lwo" and "keith" is changed
!               to "k.lwo".)
!
!   * "low_" checks to see if the filename starts with p and ends in a three
!     digit number (ex. "pkeith045.obj"), in which case the object is given a
!     color corresponding to this number.  This can be turned off by setting
!     'Default' to 1.  The default color can be changed in the section
!     marked "Default Color".
!
!   ** Currently doubled sided, smoothing, and diffusion cannot be changed.
!---------------------------------------------------------------------------------------

IMPLICIT NONE
LOGICAL, PARAMETER :: debug=.false.
CHARACTER(len=1) surfname(0:7), thousands, tens, ones
CHARACTER(len=15) flname
INTEGER(KIND=I1_Kind) position, colrlist(0:3), strchr
INTEGER npolygon
INTEGER b, c, i, j, l, m, n, npols, default
INTEGER npntsbytes, nsrfsbytes, npolsbytes, nformbytes, nsurfbytes
INTEGER(Kind=IW_Kind) ncolrbytes, nflagbytes, lwoflag, ndiffbytes, diff, npnts, nsrfs
INTEGER polslist(0:5*npnts)
INTEGER(Kind=IW_Kind) temp, one
REAL(KIND=WP_Kind) vertice(3,0:npnts-1)

default=0
npntsbytes=12*npnts

l=0
DO n=0,nsrfs-1
   l=l+polslist(l)+1
ENDDO
npols=l+nsrfs

npolsbytes =2*npols
nsrfsbytes=8
ncolrbytes=0
nflagbytes=0
ndiffbytes=0
nsurfbytes=(ncolrbytes+14+nflagbytes+6+ndiffbytes+6)
nsurfbytes=0

flname(13:15)='lwo'
WRITE(0,*)flname
!flname(0)='d'
!flname(1)='e'
!flname(2)='f'
!flname(3)='a'
!flname(4)='u'
!flname(5)='l'
!flname(6)='t'

nformbytes=(28+npntsbytes+nsrfsbytes+npolsbytes+nsurfbytes)

IF(debug)THEN
   WRITE(Out_Unit,*)'Writing LWO file'
ENDIF
!!OPEN(Unit=LWO_Unit,File=flname,form=OutDIR(1:LEN(TRIM(OutDIR)))//'binary',convert='BIG_ENDIAN') !tempmod

IF(debug)THEN
   WRITE(Out_Unit,*)"FORM"
   WRITE(Out_Unit,*)
ENDIF
WRITE(LWO_Unit)"FORM"
WRITE(LWO_Unit)nformbytes

IF(debug)THEN
   WRITE(Out_Unit,*)"LWOB"
ENDIF
WRITE(LWO_Unit)"LWOB"

IF(debug)THEN
   WRITE(Out_Unit,*)"SRFS"
   WRITE(Out_Unit,*)nsrfsbytes
ENDIF
WRITE(LWO_Unit)"SRFS"
WRITE(LWO_Unit)nsrfsbytes
surfname(0)='D'
surfname(1)='e'
surfname(2)='f'
surfname(3)='a'
surfname(4)='u'
surfname(5)='l'
surfname(6)='t'
!surfname(7)=' "

IF(flname(1:1)=='p')THEN
   surfname(0)='W'
   surfname(1)='a'
   surfname(2)='v'
   surfname(3)='e'
ENDIF

IF(default==1)THEN
   surfname(4)='g'
   surfname(5)='o'
   surfname(6)='n'
!ELSE
!   surfname(4)=thousands
!   surfname(5)=tens
!   surfname(6)=ones
ENDIF

IF(debug)THEN
   WRITE(Out_Unit,*)surfname
ENDIF
WRITE(LWO_Unit)surfname

IF(debug)THEN
   WRITE(Out_Unit,*)"PNTS"
   WRITE(Out_Unit,*)npntsbytes
ENDIF
WRITE(LWO_Unit)"PNTS"
WRITE(LWO_Unit)npntsbytes
DO j=0,npnts-1
   WRITE(LWO_Unit)(vertice(i,j),i=1,3)
   IF(debug)THEN
      WRITE(Out_Unit,*)"j=",j,"   vertice(i,j)=",(vertice(i,j),i=1,3)
   ENDIF
ENDDO


IF(debug)THEN
   WRITE(Out_Unit,*)"POLS"
   WRITE(Out_Unit,*)npolsbytes
ENDIF
WRITE(LWO_Unit)"POLS"
WRITE(LWO_Unit)npolsbytes
one=1
l=0
DO n=0,nsrfs-1
   b=l
   l=l+polslist(l)+1
   WRITE(LWO_Unit)(polslist(m),m=b,l-1),one
   IF(debug)THEN
      WRITE(Out_Unit,*)"b=",b,'  polslist(b)=',(polslist(m),m=b,l-1),one
   ENDIF
ENDDO

IF(debug)THEN
   !WRITE(Out_Unit,*)"SURF"
   !WRITE(Out_Unit,*)nsurfbytes
   !WRITE(Out_Unit,*)surfname
ENDIF
!WRITE(LWO_Unit)"SURF"
!WRITE(LWO_Unit)nsurfbytes
!WRITE(LWO_Unit)surfname

!colrlist(0)=0
!colrlist(1)=0
!colrlist(2)=0
!colrlist(3)=0
IF(debug)THEN
   !WRITE(Out_Unit,*)"COLR"
   !WRITE(Out_Unit,*)ncolrbytes
   !WRITE(Out_Unit,*)colrlist
ENDIF
!WRITE(LWO_Unit)"COLR"
!WRITE(LWO_Unit)ncolrbytes
!WRITE(LWO_Unit)colrlist

!lwoflag=256+4        ! Double Sided and Smoothing
IF(debug)THEN
   !WRITE(Out_Unit,*)"FLAG"
   !WRITE(Out_Unit,*)lwoflag
ENDIF
!WRITE(LWO_Unit)"FLAG"
!WRITE(LWO_Unit)lwoflag

!diff=256             ! 100% Diffuse
IF(debug)THEN
   !WRITE(Out_Unit,*)"DIFF"
   !WRITE(Out_Unit,*)diff
ENDIF
!WRITE(LWO_Unit)"DIFF"
!WRITE(LWO_Unit)diff

IF(debug)THEN
   WRITE(Out_Unit,*)'LWO file has been written'
ENDIF

CLOSE(Unit=LWO_Unit)
RETURN
END
