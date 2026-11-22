SUBROUTINE clebp (ia, ib, ic, id, ie, IFinal, kselc, rac)
USE FileUnits_Module
USE racah_Module
IMPLICIT REAL(Kind=WP_Kind)(a-h, o-z)
INTEGER ia, ib, ic, id, ie, IFinal, Kselc
INTEGER i1, i2, i3, i4, i5, i6, i7, i8, n, ig2, ig4
INTEGER iabc, icab, iabcp, ibca, iapd, iamd, ibpe, ibme, icpf, icmf
INTEGER ig, nzmic2, nzmic3, nzmi, nzmx, nzt1, nzt2, nzt3, nzt4, nzt5, nz, nzm1
REAL(Kind=WP_Kind) Check, fn, fb,f1, f2, s1, fc2, rac, sqfclg, termlg, ssterm
!-----------------------------------------------------------------------
!  calculates clebsch-gordan coeffs.  PARAMETERs ia thru IFinal are 2j1, 2j2
!  2j3, 2m1, 2m2, 2m3, resp.  rac is the value returned.
!  kselc=1 checks the selection rules.  kselc=0 does not.
!  t tamura, comp phys comm 1 (1970) 337.
!  to increase size of angular momenta that can be handled, just increas
!  the dimension on faclog and the initial DO loop which calculates it.
!-----------------------------------------------------------------------
DATA check /0.0d0/
!-----------------------------------------------------------------------
!  construct logs of factorials the first time thru.
!-----------------------------------------------------------------------
IF(check /= 0.0d0) GOTO 20
faclog(1)=0.0d0
faclog(2)=0.0d0
fn=1.0d0
DO 10 n=3, 500
   fn=fn+1.0d0
   faclog(n)=faclog(n-1)+log(fn)
 10   CONTINUE
check=1.0d0
 20   CONTINUE
rac=0.0d0
ig4=ia+ib+ic
IF(kselc == 0) GOTO 30
IF(id+ie /= IFinal) GOTO 90
IF(mod(ig4,2) /= 0) GOTO 90
IF(ia+ib-ic < 0 .or. ic-iabs(ia-ib) < 0) GOTO 90
IF(min0(ia-iabs(id), ib-iabs(ie), ic-iabs(IFinal)) < 0) GOTO 90
IF(mod(ib+ie, 2) /= 0 .or. mod(ic+IFinal, 2) /= 0) GOTO 90
 30   IF(ia == 0 .or. ib == 0) GOTO 40
IF(IC<0)GOTO 90
IF(IC==0)GOTO 50
IF(IC>0)GOTO 60
 40   rac=1.0d0
GOTO 90
 50   fb=ib+1
rac=((-1.0d0)**((ia-id)/2))/sqrt(fb)
GOTO 90
 60   IF(id /= 0 .or. ie /= 0) GOTO 70
ig2=ig4/2
IF(mod(ig2, 2) /= 0) GOTO 90
ig=ig2/2
i1=(ia+ib-ic)/2+1
i2=(ic+ia-ib)/2+1
i3=(ic+ib-ia)/2+1
i4=ig2+2
i5=ig+1
i6=(ig2-ia)/2+1
i7=(ig2-ib)/2+1
i8=(ig2-ic)/2+1
f1=exp(.5d0*(faclog(i1)+faclog(i2)+faclog(i3)-faclog(i4))&
   +(faclog(i5)-faclog(i6)-faclog(i7)-faclog(i8)))
f2=ic+1
f2=sqrt(f2)
s1=1-2*mod((ig2+ic)/2,2)
rac=s1*f1*f2
GOTO 90
 70   fc2=ic+1
iabcp=ig4/2+1
iabc=iabcp-ic
icab=iabcp-ib
ibca=iabcp-ia
iapd=(ia+id)/2+1
iamd=iapd-id
ibpe=(ib+ie)/2+1
ibme=ibpe-ie
icpf=(ic+IFinal)/2+1
icmf=icpf-IFinal
sqfclg=0.5d0*(log(fc2)-faclog(iabcp+1)+faclog(iabc)+faclog(icab)&
   +faclog(ibca)+faclog(iapd)+faclog(iamd)+faclog(ibpe)&
   +faclog(ibme)+faclog(icpf)+faclog(icmf))
nzmic2=(ib-ic-id)/2
nzmic3=(ia-ic+ie)/2
nzmi=max0(0, nzmic2, nzmic3)+1
nzmx=min0(iabc, iamd, ibpe)
s1=(-1.0d0)**(nzmi-1)
DO 80 nz=nzmi, nzmx
   nzm1=nz-1
   nzt1=iabc-nzm1
   nzt2=iamd-nzm1
   nzt3=ibpe-nzm1
   nzt4=nz-nzmic2
   nzt5=nz-nzmic3
   termlg=sqfclg-faclog(nz)-faclog(nzt1)-faclog(nzt2)-faclog(nzt3)&
     -faclog(nzt4)-faclog(nzt5)
   ssterm=s1*exp(termlg)
   rac=rac+ssterm
   s1=-s1
 80   CONTINUE
IF(ABS(rac) < 0.000001d0) rac=0.0d0
 90   RETURN
END
