SUBROUTINE potvb (n,r,u,up,upp,nstrt,nend)
USE FileUnits_MODULE
USE Narran_Module
USE parms_MODULE
USE nstate_MODULE
USE Arrch_Module
USE new_MODULE
USE GaussQuady_Module
USE pow_MODULE
USE Qall_Module
USE VFunc_Module
USE Masses_Module
USE melm_MODULE
USE tlogd_MODULE
USE totj_MODULE
USE approx_MODULE
IMPLICIT NONE
!#include <Inc.unitno>
!-----------------------------------------------------------------------
! calculates the matrix elements of the potential.
! written by g. a. parker
!-----------------------------------------------------------------------

LOGICAL wrtpot, redpot
LOGICAL iwrote
LOGICAL little, medium, full
INTEGER ithcll, ithsub, ixi, n, inpx, npow1, npt, npt2, ntwo, i
INTEGER kpx1, kpx, jj, lp1, l, index, nqp, nstrt, nend, jp, lbp, megp
INTEGER nujmp, nq, nujm, imeg, lmin1, j, lb, meg, lmax1, kp, in, inp
INTEGER iii, k, ik, in1, in2, in3, in4, ms, npwpt2, kp1, lengs
REAL(Kind=WP_Kind) alphai, r1, r2, r2mr1, a, b, ang, x, pn0, pn1, pn2, xp
REAL(Kind=WP_Kind) f, ccw, coef, cleb, percvl, cgama, cgamap, rb, fac2, r
REAL(Kind=WP_Kind) summ1, summ2, summ3, sum1, sum2, sum3, rr
REAL(Kind=WP_Kind) xptg(2000), whtg(2000)   ! PORT: enlarged from 96 (large nglegn)
REAL(Kind=WP_Kind) pv(1000), pvxp(4000), p(1000), pnwv(8000), pnw(8000)   ! PORT: enlarged
REAL(Kind=WP_Kind) u(1), up(1), upp(1)
REAL(Kind=WP_Kind) xptv(2000), whtv(2000)   ! PORT: enlarged from 96
REAL(Kind=WP_Kind) v(4000), vp(4000), vpp(4000)   ! PORT: enlarged from 500 (large nfreq)
REAL(Kind=WP_Kind) pot(4000), potp(4000), potpp(4000)   ! PORT: enlarged from 500
REAL(Kind=WP_Kind) s(5000)
DATA lengs /5000/, wrtpot /.false./
DATA redpot /.false./
DATA ithcll/0/,ithsub/0/
DATA little/.false./,medium/.false./,full/.false./
CALL popt('potvb   ',little,medium,full,ithcll,ithsub)
IF(full)THEN
WRITE(Out_Unit,*)' kchan mvib jrot lorb, nuj'
DO ixi=1,n
   inpx=1+(nuj3(ixi)-1)*nhermt(karran)
   WRITE(Out_Unit,3738)kchan(ixi),mvib3(ixi),jrot3(ixi),lorb3(ixi), nuj3(ixi),chinuj(inpx)
ENDDO
3738 FORMAT(1x,5i5,e14.7)
ENDIF
IF(.NOT.inext(karran)) GOTO 10
IF(redpot) GOTO 260
GOTO 200
10 iwrote=.false.
alphai=1.0d0/alpha(karran)
r1=re(karran)*(alphai*xpth(1,karran)+rx(karran))
r2=re(karran)*(alphai*xpth(nhermt(karran),karran)+rx(karran))
IF(nhermt(karran)==1) r1=0.d0
IF(nhermt(karran)==1) r2=2.d0*re(karran)
r2mr1=r2-r1
b=1.0d0
a=-1.0d0
npow1=npow(karran)+1
npt=nglegn(karran)+1
npt2=(npt+1)/2
IF(.NOT.symmetry) npt2=npt
ntwo=2
IF(.NOT.symmetry) ntwo=1
CALL glegen (npt,xptg,whtg,a,b)
CALL glegen (npow1,xptv,whtv,a,b)
IF(.NOT.symmetry) GOTO 20
IF((npt/2)*2/=npt) whtg(npt2)=0.5d0*whtg(npt2)
20 CONTINUE
IF(medium)THEN
WRITE(Out_Unit,340) nglegn(karran),npt,npt2
WRITE(Out_Unit,360)
DO i=1,npt2
   ang=acos(xptg(i))
   WRITE(Out_Unit,350) i,ang,xptg(i),whtg(i)
ENDDO
ENDIF
REWIND sb_Unit
!-----------------------------------------------------------------------
! determine constant integrals of vibrational &
! angular legendre polynomials in spherical
! harmonic-vibrational eigenfunction basis
!-----------------------------------------------------------------------
! determine vibrational legendre polynomials at all
! Gauss_Hermite integration points
!----------------------------------------------------------------------
kpx1=1
DO i=1,nhermt(karran)
kpx=kpx1
kpx1=kpx1+1
x=re(karran)*(alphai*xpth(i,karran)+rx(karran))
x=(2.d0*x-r1-r2)/r2mr1
pv(1)=1.d0
pv(2)=x
IF(npow1<=2) GOTO 50
pn0=1.d0
pn1=x
DO jj=3,npow1
   xp=x*pn1
   pn2=xp-pn0+xp-(xp-pn0)/REAL(jj-1,WP_Kind)
   pn0=pn1
   pn1=pn2
   pv(jj)=pn1
ENDDO
50 DO lp1=1,npow1
   l=lp1-1
   pvxp(kpx)=pv(lp1)*sqrt((2*l+1)/2.d0)
   kpx=kpx+nhermt(karran)
ENDDO
ENDDO
!------------------------------------------------------------------
!  loop over channels
!------------------------------------------------------------------
index=0
DO nqp=nstrt,nend
jp=jrot3(nqp)
lbp=lorb3(nqp)
megp=mega3(nqp)
nujmp=(nuj3(nqp)-1)*nhermt(karran)
DO nq=nstrt,nqp
j=jrot3(nq)
lb=lorb3(nq)
meg=mega3(nq)
nujm=(nuj3(nq)-1)*nhermt(karran)
IF(melmnt=='probf')THEN
IF(meg/=megp) GOTO 110
imeg=iabs(meg)
IF(imeg>j) GOTO 110
IF(imeg>jp) GOTO 110
lmin1=iabs(j-jp)+1
lmax1=min0(j+jp,npt-1)+1
ENDIF
IF(melmnt=='space')THEN
lmin1=max0(iabs(lb-lbp),iabs(j-jp))+1
lmax1=min0(j+jp,npt-1,lb+lbp)+1
ENDIF
IF(melmnt=='molbf')THEN
IF(meg/=megp) GOTO 110
imeg=iabs(meg)
IF(imeg>lb) GOTO 110
IF(imeg>lbp) GOTO 110
lmin1=iabs(lb-lbp)+1
lmax1=min0(lb+lbp,npt-1)+1
ENDIF
!-----------------------------------------------------------------------
! loop over vibrational legendre polynomials
!-----------------------------------------------------------------------
kpx=1
DO kp=1,npow1
f=0.d0
!-----------------------------------------------------------------------
! loop over the  Gauss_Hermite integration points.
!-----------------------------------------------------------------------
DO i=1,nhermt(karran)
   inp=i+nujmp
   in=i+nujm
   ccw=chinuj(inp)*chinuj(in)*wpth(i,karran)
   f=f+pvxp(kpx)*ccw
   kpx=kpx+1
ENDDO
IF(ABS(f)<1.d-10) f=0.d0
!-----------------------------------------------------------------------
! loop over the legendre polynomials.
!-----------------------------------------------------------------------
DO lp1=lmin1,lmax1,ntwo
l=lp1-1
index=index+1
IF(index<=lengs) GOTO 90
IF(.NOT.iwrote) REWIND sb_Unit
iwrote=.true.
WRITE(sb_Unit) s
index=1
90 IF(melmnt=='probf') coef=sqrt(REAL(2*j+1,WP_Kind)/REAL(2*jp+1,WP_Kind)) &
*cleb(j,l,jp,meg,0,meg)*cleb(j,l,jp,0,0,0)*sqrt((2*l+1)/2.d0)
IF(melmnt=='space') coef=percvl(l,j,lb,jp,lbp,jtot)*sqrt((2*l+1)/2.d0)
IF(melmnt=='molbf') coef=sqrt(REAL(2*lb+1,WP_Kind)/REAL(2*lbp+1,WP_Kind)) &
*cleb(lb,l,lbp,meg,0,meg)*cleb(lb,l,lbp,0,0,0)*sqrt((2*l+1)/2.d0)
IF(ABS(coef)<1.d-10) coef=0.d0
s(index)=f*coef
ENDDO
ENDDO
110 CONTINUE
ENDDO
ENDDO
IF(full)THEN
WRITE(Out_Unit,*)' The elements of S are:  ',(s(iii),iii=1,index)
ENDIF
IF(iwrote) WRITE(sb_Unit) s
IF(redpot) GOTO 260
!-----------------------------------------------------------------------
! calculate vibrational & angular legendre polynomials
! for expansion of potential
!-----------------------------------------------------------------------
in4=1
DO k=1,npow1
in2=in4
in4=in4+1
cgama=xptv(k)
pv(1)=1.d0
pv(2)=cgama
IF(npow1<=2) GOTO 130
pn0=1.d0
pn1=cgama
DO j=3,npow1
   cgamap=cgama*pn1
   pn2=cgamap-pn0+cgamap-(cgamap-pn0)/REAL(j-1,WP_Kind)
   pn0=pn1
   pn1=pn2
   pv(j)=pn1
ENDDO
130 CONTINUE
DO lp1=1,npow1
   l=lp1-1
   pv(lp1)=pv(lp1)*sqrt((2*l+1)/2.d0)
ENDDO
DO lp1=1,npow1
   pnwv(in2)=pv(lp1)*whtv(k)
   in2=in2+npow1
ENDDO
ENDDO
IF(full) WRITE(Out_Unit,440)
IF(full) CALL MxOutD(pnwv,npow1,npow1,0)
fac2=1.d0
IF(symmetry) fac2=2.d0
in4=1
DO k=1,npt2
in2=in4
in4=in4+1
cgama=xptg(k)
p(1)=1.d0
p(2)=cgama
IF(npt<=2) GOTO 170
pn0=1.d0
pn1=cgama
DO j=3,npt
   cgamap=cgama*pn1
   pn2=cgamap-pn0+cgamap-(cgamap-pn0)/REAL(j-1,WP_Kind)
   pn0=pn1
   pn1=pn2
   p(j)=pn1
ENDDO
170 CONTINUE
DO lp1=1,npt,ntwo
   l=lp1-1
   p(lp1)=p(lp1)*sqrt((2*l+1)/2.d0)*fac2
ENDDO
DO lp1=1,npt,ntwo
   pnw(in2)=p(lp1)*whtg(k)
   in2=in2+npt2
ENDDO
ENDDO
IF(full) WRITE(Out_Unit,450)
IF(full) CALL MxOutD(pnw,npt2,npt2,0)
!----------------------------------------------------------------------
! project out vibrational & angular legendre polynomials from potential
! expand the potential in a vibrational legendre polynomial basis set
!-----------------------------------------------------------------------
200 CONTINUE
ik=0
DO i=1,npow1
   rb=(r2mr1*xptv(i)+r2+r1)*0.5d0
   DO k=1,npt2
      ik=ik+1
      cgama=xptg(k)
      CALL potabc (r,rb,cgama,v(ik),vp(ik),vpp(ik))
   ENDDO
ENDDO
IF(ik>500) GOTO 320
in1=1-npow1
index=0
DO kp=1,npow1
in1=in1+npow1
!-----------------------------------------------------------------------
! expand the potential in an angular legendre polynomial basis set.
!-----------------------------------------------------------------------
in2=1-npt2
DO lp1=1,npt,ntwo
in2=in2+npt2
l=lp1-1
index=index+1
!-----------------------------------------------------------------------
! loop over the vibrational Gauss_Legendre integration points.
!-----------------------------------------------------------------------
summ1=0.0d0
summ2=0.0d0
summ3=0.0d0
ik=0
DO i=1,npow1
in3=in1+i-1
!-----------------------------------------------------------------------
! loop over the angular Gauss_Legendre integration points.
!-----------------------------------------------------------------------
sum1=0.0d0
sum2=0.0d0
sum3=0.0d0
DO k=1,npt2
   ik=ik+1
   in4=in2+k-1
   sum1=sum1+v(ik)*pnw(in4)
   IF(numder) GOTO 220
   sum2=sum2+vp(ik)*pnw(in4)
   sum3=sum3+vpp(ik)*pnw(in4)
220 CONTINUE
ENDDO
summ1=summ1+pnwv(in3)*sum1
IF(numder) GOTO 230
summ2=summ2+pnwv(in3)*sum2
summ3=summ3+pnwv(in3)*sum3
230 CONTINUE
ENDDO
IF(index>4000)THEN
    WRITE(Out_Unit,*)'error: index>4000',index
    STOP 'potvb'
ENDIF
pot(index)=summ1
IF(numder) GOTO 240
potp(index)=summ2
potpp(index)=summ3
240 CONTINUE
ENDDO
ENDDO
!-----------------------------------------------------------------------
! IF desired print out the expansion coefficients
! of the potential in the vib-rot basis
!-----------------------------------------------------------------------
IF(full)THEN
WRITE(Out_Unit,370) n,r
ms=0
WRITE(Out_Unit,380)
CALL MxOutD(pot,npt2,npow1,ms)
IF(numder) GOTO 250
WRITE(Out_Unit,390)
CALL MxOutD(potp,npt2,npow1,ms)
WRITE(Out_Unit,400)
CALL MxOutD(potpp,npt2,npow1,ms)
ENDIF
250 npwpt2=npow1*npt2
inext(karran)=.true.
IF(wrtpot) WRITE(Diatomic_Pot_Unit) npwpt2,(pot(i),potp(i),potpp(i),i=1,npwpt2),r,npow1,npt,ntwo,npt2
!-----------------------------------------------------------------------
! determine matrix elements of the potential and its
! derivatives.
!-----------------------------------------------------------------------
260 IF(.NOT.redpot) GOTO 270
READ(Diatomic_Pot_Unit) npwpt2,(pot(i),potp(i),potpp(i),i=1,npwpt2),rr,npow1,npt,ntwo,npt2
inext(karran)=.true.
IF(rr/=r) GOTO 320
270 IF(iwrote) REWIND sb_Unit
in2=0
in=0
DO nqp=nstrt,nend
jp=jrot3(nqp)
megp=mega3(nqp)
lbp=lorb3(nqp)
DO nq=nstrt,nqp
j=jrot3(nq)
meg=mega3(nq)
lb=lorb3(nq)
in2=nq+nqp*(nqp-1)/2
u(in2)=0.d0
up(in2)=0.d0
upp(in2)=0.d0
IF(melmnt=='probf')THEN
IF(meg/=megp) GOTO 300
imeg=iabs(meg)
IF(imeg>j) GOTO 300
IF(imeg>jp) GOTO 300
lmin1=iabs(j-jp)+1
lmax1=min0(j+jp,npt-1)+1
ENDIF
IF(melmnt=='space')THEN
lmin1=max0(iabs(lb-lbp),iabs(j-jp))+1
lmax1=min0(j+jp,npt-1,lb+lbp)+1
ENDIF
IF(melmnt=='molbf')THEN
IF(meg/=megp) GOTO 300
imeg=iabs(meg)
IF(imeg>lb) GOTO 300
IF(imeg>lbp) GOTO 300
lmin1=iabs(lb-lbp)+1
lmax1=min0(lb+lbp,npt-1)+1
ENDIF
sum1=0.d0
sum2=0.d0
sum3=0.d0
DO kp1=1,npow1
   DO lp1=lmin1,lmax1,ntwo
      in1=(lp1+1)/ntwo+(kp1-1)*npt2-1/ntwo
      in=in+1
      IF(.NOT.iwrote) GOTO 280
      IF(in>lengs) in=1
      IF((in==1)) READ(sb_Unit) s
  280 sum1=sum1+pot(in1)*s(in)
      IF(numder) GOTO 290
      sum2=sum2+potp(in1)*s(in)
      sum3=sum3+potpp(in1)*s(in)
  290 CONTINUE
   ENDDO
ENDDO
u(in2)=usys2*sum1
IF(numder) GOTO 300
up(in2)=usys2*sum2
upp(in2)=usys2*sum3
300 CONTINUE
ENDDO
ENDDO
!-----------------------------------------------------------------------
! IF desired print the matrix elements of the potential
! and its derivatives.
!-----------------------------------------------------------------------
IF(full)THEN
ms=1
WRITE(Out_Unit,410)
CALL MxOutL(u,n,n,ms,melmnt,melmnt)
IF(numder) GOTO 310
WRITE(Out_Unit,420)
CALL MxOutL(up,n,n,ms,melmnt,melmnt)
WRITE(Out_Unit,430)
CALL MxOutL(upp,n,n,ms,melmnt,melmnt)
ENDIF
310 RETURN
320 WRITE(Out_Unit,330) rr,r,ik
STOP 'potvb'
!-----------------------------------------------------------------------
!----------------***end-potvb***---------------------------------------
!-----------------------------------------------------------------------
!
330 FORMAT(//,"*****error*****","rr=",e14.7,5x,"r=",e14.7,5x,"ik=",i5)
340 FORMAT(///,"nglegn=",i5///,"npt=",i5///,"npt2=",i5)
350 FORMAT(//,i5,5x,e15.7,5x,e15.7,5x,e15.7)
360 FORMAT(//,4x,"i",8x,"theta",12x,"cos(theta)",12x,"weight")
370 FORMAT(///,"no. of channels=",i5///,"r=",f10.5)
380 FORMAT(//,"   vib-rot expansion of the potential")
390 FORMAT(//,"            vib-rot expansion of the first derivative of the potential")
400 FORMAT(//,"            vib-rot expansion of the second derivative of the potential")
410 FORMAT(//,"potential energy matrix elements")
420 FORMAT(//,"first derivative  of the potential energy matrix elements")
430 FORMAT(//,"second derivative of the potential energy matrix elements")
440 FORMAT(//,"vibrational legendre poly. x g.h. weight")
450 FORMAT(//,"legendre poly. x g.l. weight")
END
