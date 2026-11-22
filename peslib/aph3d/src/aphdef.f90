SUBROUTINE aphdef(nx, iy, ic, nxdim, nydim, x, snthd, thaph, chi1, chif1, betqf, thetd)
!
!this routine is called by: AphDel_Old or AphDel_New
!
!-----------------------------------------------------------------------
USE FileUnits_Module
USE Numeric_Kinds_Module
USE Numbers_Module
USE Narran_Module
IMPLICIT NONE
INTEGER nx,iy,ic,nxdim,nydim
REAL(Kind=WP_Kind) x(*), snthd(nydim,3), thaph(nxdim,nydim,3)
REAL(Kind=WP_Kind) chi1(nxdim,nydim,3), chif1(*)
REAL(Kind=WP_Kind) betqf(nxdim,nydim,3)
REAL(Kind=WP_Kind) thetd
REAL(Kind=WP_Kind) ththd,tanthd,snapth,qtemp,xx,chif,xy
INTEGER ix

!----------------------------------------------------------------
snthd(iy,ic)=sin(2.d0*thetd)
ththd=tan(2.d0*thetd)
tanthd=tan(thetd)
!     calculate aph theta and chi and also betaqf at each quadrature point.
DO  ix=1, nx
   snapth=sqrt(1.d0-x(ix)**2)
   qtemp=snthd(iy,ic)
   xx=qtemp*snapth                  !     xx is a working variable
   thaph(ix,iy,ic)=HalfPi-1.d-10
   IF(xx == 0.d0)GOTO 141
   thaph(ix,iy,ic)=atan(sqrt(1.d0-xx**2)/xx)
 141     CONTINUE
   xx=ththd*x(ix)
   chif=0.5d0*acos(1.0d0/sqrt(1.d0+xx**2))
   chif=sign(chif,x(ix))
   chi1(ix,iy,ic)=chif+chif1(ic)
   xy=chi1(ix,iy,ic)
   IF(ABS(xy) > pi)WRITE(Out_Unit,663)ix, iy, ic, xy
 663     FORMAT(1x, ' at ix,iy,ic=', 3i5, ' chi1=', e16.9)
   xx=tanthd*tan(chif)
   betqf(ix,iy,ic)=atan(xx*snapth/(1.d0+xx*x(ix)))
   xx=thaph(ix,iy,ic)

   IF(xx < 0.d0 .or. xx > HalfPi)THEN
      WRITE(Out_Unit,662) ix, iy, ic, xx, chi1(ix,iy,ic), betqf(ix,iy,ic)
   ENDIF

 662     FORMAT(1x, ' at ix,iy,ic=', 3i5, ' thaph,chi1,betqf=', 3e13.6)
ENDDO
RETURN
ENDSUBROUTINE aphdef
