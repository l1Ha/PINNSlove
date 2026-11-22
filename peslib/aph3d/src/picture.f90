SUBROUTINE picture(id, th, ch, phifun, phirms, igloss, kolr,  center, &
         rad, nelem, numnp, nelemx, ntheta,  nchi, nidm, chivals, thetaval, naph)
USE Numeric_Kinds_Module
USE PES_MODULE
Use FileUnits_Module
USE Converge_Module
USE point_Module
USE Narran_Module
USE opts_Module
USE STST_Module
USE main_Module
USE convrsns_Module
USE rhosur_Module
USE totj_Module
USE chiang_Module
USE pcture_Module
USE Numbers_Module
USE delang_Module
USE mdata_Module
USE Masses_Module
USE PlotNL_Module
!
!   $RCSfile: picture.f,v $ $Revision: 1.15 $
!   $Date: 89/11/16 11:03:18 $
!   $State: Stable $
!
!----------------------------------------------------------------
!  id       list of node numbers.
!  th       list of theta values.
!  ch       list of chi values.
!  nelem    number of elements.
!  numnp    number of nodal points.
!  nelemx   For symmetry=false nelemx=nelem.
!           For symmetry=true, nelemx=2*nelem.
!  phifun   surface functions.
!  phi      scratch array.
!  igloss   scratch array.
!  kolr     scratch array.
!  center   scratch array.
!  rad      scratch array.
!  phirms
!  ptheta
!  nchi
!  nidm
!  chivals
!  thetaval
!  naph

IMPLICIT NONE
LOGICAL ifirst
INTEGER irho, nthmod, nelem, numnp, nelemx, ntheta, nchi, nidm, naph, icnt
INTEGER id(9,nelem), th(numnp), ch(numnp), igloss(nelemx), kolr(nelemx)
REAL(Kind=WP_Kind) Etot
REAL(Kind=WP_Kind) center(3,nelemx), rad(nelemx), chivals(nchi), thetaval(ntheta)
REAL(Kind=WP_Kind) phifun(numnp,naph), phirms(ntheta,nchi)

!----------------------------------------------------------------
! set default colors
! hue=0.0  black
! hue=0.5  magenta
! hue=1.0  red
! hue=1.5  yellow
! hue=2.0  green
! hue=2.5  cyan
! hue=3.0  blue
! hue=3.5  magenta
! hue=4.0  white
!
!DATA hue/0.0d0, 0.5d0, 1.0d0, 1.5d0, 2.0d0, 2.5d0, 3.0d0, 3.5d0, 4.0d0, 0.5d0, 1.5d0/
!----------------------------------------------------------------
! nil gives the number of intensity levels for shaded surface pictures.
!     nil = 30 gives good quality pictures. nil = 10 is satisfactory for
!     many purposes.  more computer time is required for larger nil.
!     maximum nil is 100.
!DATA nil/10/
!----------------------------------------------------------------
! view is a vector of size 3 giving the x, y, and z coordinates
!     of the view point.
!DATA view/0.d0,-50.d0,50.d0/
!----------------------------------------------------------------
! win is a vector of size 4 defining the x, y, and z coordinates of
!     the center of the viewing window.  win(4) defines the width of the
!     window.
!DATA win/0.d0,0.d0,0.d0,8.d0/
!----------------------------------------------------------------
! sun=xsn is a 4 by nxsn array giving the x, y, and z coordinates of
!     the illuminating light positions.  a good position for the light
!     is at the view position.  xsn(4,n) gives the light intensity of
!     the nth light source (values between 0.0 and 1.0).  this should
!     normally be set to 1.0, but it should be less for multiple light
!      sources.
!DATA xsn/-50.d0,0.d0,200.d0,8.d0/
!----------------------------------------------------------------
! nxsn gives the number of light sources.
!DATA nxsn/1/
!----------------------------------------------------------------
! kolr is a vector of size nob that points to the colors in a color
!     palette.  example, IF kolr(3) = 6, THEN object 3 will be colored
!     with hue(6), sat(6), and cint(6).  IF kolr(k) = 0, the colors of
!     individual cells of object k will be colored according to numbers
!     in last 9 bits of the y coordinates of object k as follows.
!     IF the last 9 bits of y(i,j) in object k represent the number n,
!     THEN cell (i,j) will be colored hue(n), sat(n), cint(n).  this can
!     be set up by the following fortran statement,
!          y(i,j) = y(i,j).AND..NOT.777b.or.n
!     IF kolr is negative, object will be colored according to distance
!     from a point (see rad and center below).  this does not work
!     for mesh plots.  example: IF kolr(3) = -5, THEN the color on
!     object 3 within radius rad(5) from a point defined by the
!     x,y,z coordinates center(1,3), center(2,3), and center(3,3) will
!     be hue(5), sat(5),cint(5).  points on the object that are farther
!     than rad(5) but less than rad(6) are colored hue(6), sat(6), etc.
!     any number of monotonically increasing radii (rad) may be used for
!     for an object but the last one must be greater than the distance
!     from center to the farthest point on the object.  otherwise
!     when the program looks for the next larger rad, it may go into
!     an endless loop.  different sets of radii may be defined for
!     different objects, or the same set may be used.
!DATA klrplt/4/
!----------------------------------------------------------------
! sat is an array corresponding to hue which gives saturation
!     of the colors.  sat = 1. gives pure colors.  sat = 0. gives
!     white.  intermediate values give pastel colors.
!DATA sat/1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0, 0.5d0,0.0d0/
!----------------------------------------------------------------
! cint is an array corresponding to hue that gives
!     relative intensity of reflected light.  cint must be in
!     the range 0.0 to 1.  normally it should be 1.
!DATA cint/1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,0.5d0, 0.d0/
!----------------------------------------------------------------
! rad is a table of radii that is used only when kolr is negative.
!DATA radii/0.0/
!----------------------------------------------------------------
! center is a table that defines the coordinates of the center
!     for distances rad.  center(1,n), center(2,n), and center(3,n)
!     are the x,y, and z coordinates of the point for object n.
!     center is used only when kolr(n) is negative.
!DATA cent/0.d0,0.d0,0.d0/
!----------------------------------------------------------------
! igloss is a flag that tells whether objects will have specualr
!     reflection (gloss).  IF igloss(n) = 0, object n will have
!     diffuse reflection.  IF igloss(n) = 1, object n will also
!     have specular reflection.
!DATA kgloss/1/
!----------------------------------------------------------------
! shad determines shadowing.  set shad(5) = 0.0 for no shadowing.
!     set shad(5) = 1.0 for darkest shadows.  intermediate values
!     give intermediate darkness.  shadows will be determined by
!     temporarily setting up viewpoints at the sun positions and
!     calculating hidden portions of the picture.  a window must
!     be set for the shadowing.  the x, y, and z coordinates of the
!     center of this window are given by shad(1), shad(2), shad(3),
!     respectively.  the width of the window is given by shad(4).
!     points outside this window (as viewed from the sun position) will
!     not be shadowed.  make the window large enough to DO the job, but
!     IF it is too large, the resolution becomes too course.
!     the roughness around the edges of the shadows is due to
!     the courseness of the mesh.
!     IF shad(5) = 0.0, the other values of the shad array need not
!     be defined.  multiple shadows are possible for multiple light
!     sources.  more time is required for more suns.
!DATA shad/0.d0, 0.d0, 0.d0, 0.d0, 0.d0/
!----------------------------------------------------------------
! nl specifies the number of scan lines.
!     nl=300 is best for a shaded surface on the tektronix.
!     nl=750 is good for the fr-80.
!     when iplot = 2, no scan lines are drawn, but nl defines
!     the resolution for the drawing of the mesh plot.
!DATA nl/300/
!----------------------------------------------------------------
! kolb is a pointer that defines the background color on the palette
!     hue, saturation, and intensity.  kolb(1) defines the color at the
!     top of the screen, kolb(2) gives the color at bottom, and
!     intermediate points have intermediate colors.  example, kolb = 7
!     gives background color hue(7), sat(7), and cint(7).  IF kolr(1) =
!     0, no background will be given.  background is not given for mesh
!     plots.
!DATA kolb/10 ,11/
!----------------------------------------------------------------
! iplot specifies the type of plot.  iplot = 1 gives shaded surfaces.
!     iplot = 2 gives mesh plot.  iplot=3 gives both shaded surfaces
!     and mesh plot.
!DATA iplt/2/
!DATA iproj/2/
!DATA denplot/.false./
!DATA kase/1/
!DATA chij/0.d0/
!DATA ifirst/.true./
!DATA range/1.0d0/
!DATA logpot/.false./
!DATA icnt/9999/
IF(icnt==9999)THEN
   OPEN(Unit=Pictre_Bin_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'Pictre.bin',form='unformatted', status='unknown')
   icnt=0
ENDIF
WRITE(Pictre_Bin_Unit)nelem, numnp, nelemx, ntheta, nchi, nidm, naph, mass(1), mass(2), mass(3), PES_Name, scheme, betaw, chncof, seq
WRITE(Pictre_Bin_Unit)chivals
WRITE(Pictre_Bin_Unit)thetaval
WRITE(Pictre_Bin_Unit)ndiv
WRITE(Pictre_Bin_Unit)autoev, amutoau, angtoau, radtodeg, cmm1toau
WRITE(Pictre_Bin_Unit)rhosurf
WRITE(Pictre_Bin_Unit)Etot, parity, irho, jtot, mtot, megamin, megamax
WRITE(Pictre_Bin_Unit)chif1
!WRITE(Pictre_Bin_Unit)view, xsn, hue, sat, cint, shad, kolb, cent, win, logpot,        &
!        nthmod, nl, nsurf, nil, irdsurf, denplot, iproj, klrplt,          &
!        iplt, kgloss, cutoffl, radii, chij, kase, range, cutoffh, npltskp
WRITE(Pictre_Bin_Unit)fourthpi, halfpi, pi, twopi
WRITE(Pictre_Bin_Unit)dthe, dchi
WRITE(Pictre_Bin_Unit)nlt, nlc, nfreq, nq, nmode, rtol, nodmax, nblock, maxj, nlancz
WRITE(Pictre_Bin_Unit)id
WRITE(Pictre_Bin_Unit)th
WRITE(Pictre_Bin_Unit)ch
WRITE(Pictre_Bin_Unit)phifun
WRITE(Pictre_Bin_Unit)phirms
CLOSE(Unit=Pictre_Bin_Unit)
RETURN
ENDSUBROUTINE picture
