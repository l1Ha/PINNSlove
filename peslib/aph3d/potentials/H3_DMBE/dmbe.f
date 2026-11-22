      subroutine surface_dmbe(v,rij)
c-----------------------------------------------------------------------
c This routine supplies the potential energy surface for the
c reactive scattering codes.
c
c On entering: The arrary rij contains the internuclear distances
c              in units of bohr (NOT mass scaled).
c
c              The common /pes/ must contain a character string
c              of length 25 that properly describes the potential.
c
c On Exit:     The variable v is the value of the potential energy
c              surface in Hartree atomic units.
c
c
c-----------------------------------------------------------------------
      USE PES_MODULE
      implicit real*8 (a-h,o-z)
      dimension rij(3)
      character(LEN=25)potname
c      common /pes/ pesname
c-----------------------------------------------------------------------
c  The character string potname must contain a string that uniquely
c  specifies the potential energy surface to be used.
c-----------------------------------------------------------------------
      data potname/'DMBE HHH         PES     '/
c-----------------------------------------------------------------------
c Check for correct potential energy surface specification.
c-----------------------------------------------------------------------
      if(pes_name.ne.potname)then
      WRITE(0,*)'The potential energy surface name does not match',
     >          ' the name in the potential energy surface routine.'
      WRITE(0,304)pesname,potname
      WRITE(0,*)'***error***: Execution stopped in routine surface'    
          
      WRITE(6,*)'The potential energy surface name does not match',
     >          ' the name in the potential energy surface routine.'
      WRITE(6,304)pesname,potname
304   format(1h ,2a25)
      WRITE(6,*)'***error***: Execution stopped in routine surface'
      stop
      end if
c-----------------------------------------------------------------------
c This call to surfac supplies the potential energy surface
c for the particular system of interest.
c-----------------------------------------------------------------------
      call dmbe(v,rij)
      RETURN
      ENDSUBROUTINE Surface_dmbe
*dk dmbe
      subroutine dmbe(v,r1)
      implicit real*8 (a-h,o-z)
      dimension r1(3)
      common/potcm/r(3),pote1,dedr(3)
      common/pot2cm/pote2,dedr2(3),coup(3)
      r(1)=r1(1)
      r(2)=r1(2)
      r(3)=r1(3)
      call h3dmbe
      v=pote1+.1744741119d0
      return
      end
      subroutine h3dmbe
c***********************************************************************
c     calculates double many-body expansion of the h3 potential
c        pote1 is energy of surface 1.
c        dedr contains derivatives for surface 1.
c        pote2 is energy of surface 2.
c        dedr2 contains derivatives for surface 2.
c        coup contains nonadiabatic coupling.
c        note:  dedr,dedr2,coup not yet coded.
c***********************************************************************
c  refs.:  b. liu & p. siegbahn, j. chem. phys., 68, 2457, 1978;
c          d.g. truhlar & c.j. horowitz, j. chem. phys., 68, 2466, 1978;
c          a.j.c. varandas, mol. phys., 53, 1303, 1984;
c          a.j.c. varandas, j. mol. struct. (theochem), 21, 401, 1985;
c          a.j.c. varandas, d.g. truhlar, c.a. mead, & n.c. blais, to be
c             published.
c***********************************************************************
      implicit real*8 (a-h,o-z)
      common /acdcom/alph2,az2,ce0,cd1
      common /crdcom/ r12,r22,r32,qcoord,rho,s,cphi3,s2,s3
      common /h2com/rm,r0,rmt
      common /fixcom/alpha5,beta1,beta2,beta3
      common /picom/pi,sqrt3
      common /potcm/r(3),pote1,dedr(3)
      common /pot2cm/ pote2,dedr2(3),coup(3)
      common /tercom/ vhf,hflinc,bendt,hfadd,vhfc,ce2,ce3,corr
      common /lepcom/ vleps2,dlep(3),dlep2(3)
      common /xcocom/x(15)
      common /dampp/damps(10,3,2)
      common /corr/corr2s(3,2)
      external h3data
      data rhol /2.4848d0/
      r1=r(1)
      r2=r(2)
      r3=r(3)
      call coordx(r1,r2,r3)
      per=r1+r2+r3
      per2=per*per
      per3=per2*per
      exp1=exp(-beta1*per)
      exp2=exp(-beta2*per2)
      exp3=exp(-beta3*per)
      hflinc=0.0d0
      do 1 j=1,4
         hflinc=hflinc+x(j)*alin(r1,r2,r3)**(2*j)
    1 continue
      hflinc=hflinc*exp(-alpha5*per3)
c  construct damping function matrix damps(power,distance, state)
      call damp(r)
c  construct 2-body correlation fcn matrix corr2s(distance, state)
      call corr2(r)
      vhf=vleps(r1,r2,r3)
      b1=b1fn(r1,r2,r3)
      b12=b1**2
      b13=b12*b1
      b14=b13*b1
      bend1=b1*(x(5)+x(6)*per)*exp1
      bend2=(x(7)*b12+x(8)*b13+x(9)*b14)*exp2
      bend3=(1.d0/r1+1.d0/r2+1.d0/r3)*(x(10)*b1*exp1+x(11)*
     1      b12*exp2)
      bend4=((r2-r1)**2+(r3-r1)**2+(r3-r2)**2)*
     1      b1*(x(12)*exp1+x(13)*exp2)
      bend5=b1*(x(14)+x(15)*per2)*exp3
      bendt=bend1+bend2+bend3+bend4+bend5
   11 hfadd=s2*(1.d0+s3*cphi3)*(ce0+cd1*rho)*exp(-alph2*
     *      (rho-rhol)**2)
   45 vhfc=vhf+hflinc+bend1+bend2+bend3+bend4+bend5+hfadd
      call corr3(r1,r2,r3,rm,r0,ce3)
      ce2=corr2s(1,1)+corr2s(2,1)+corr2s(3,1)
      corr=ce2+ce3
      pote1=vhfc+corr
      pote2=pote1-vhf+vleps2
      return
      end
      subroutine prepot
      implicit real*8 (a-h,o-z)
c   dummy subprogram
      return
      end
      function alin(r1,r2,r3)
c***********************************************************************
c     computes base function for the second fit to the hartree-fock
c     energy at the linear symmetrical geometries.
c***********************************************************************
      implicit real*8 (a-h,o-z)
      alin=(r1-r2)*(r2-r3)*(r3-r1)
      return
      end
      function b1fn(r1,r2,r3)
c***********************************************************************
c     computes b1 function.
c***********************************************************************
      implicit real*8 (a-h,o-z)
      t1=(r1**2-r2**2-r3**2)/(2.d0*r2*r3)
      t2=(r2**2-r1**2-r3**2)/(2.d0*r1*r3)
      t3=(r3**2-r1**2-r2**2)/(2.d0*r1*r2)
      b1fn=1.d0+t1+t2+t3
      return
      end
      subroutine corr2(r)
c***********************************************************************
c     calculates 2-body correlation energy.
c***********************************************************************
      implicit real*8 (a-h,o-z)
      common /noncom/chh(10),chhe(10),chh2(10)
      common /dampp/damps(10,3,2)
      common /corr/corr2s(3,2)
      dimension ri(10),r(3)
      do 1 id=1,3
      ri(2)=1.0d0/(r(id)*r(id))
      ri(4)=ri(2)*ri(2)
      ri(6)=ri(4)*ri(2)
      ri(8)=ri(6)*ri(2)
      ri(10)=ri(8)*ri(2)
      do 1 is=1,2
      corr2s(id,is)=0.0d0
      do 1 i=6,10,2
         corr2s(id,is)=corr2s(id,is)-damps(i,id,is)*chh(i)*ri(i)
    1 continue
      return
      end
      function corr2d(r,rm,r0)
c***********************************************************************
c     calculates 1st derivative of the 2-body correlation energy.
c***********************************************************************
      implicit real*8 (a-h,o-z)
      common /noncom/chh(10),chhe(10),chh2(10)
      corr2d=0.0d0
      do 1 i=6,10,2
         corr2d=corr2d-(chh(i)/r**i)*(damdamp(r,i,rm,r0)-float(i)*
     *          damdamp(r,i,rm,r0)/r)
    1 continue
      return
      end
      subroutine damp(r)
c***********************************************************************
c  calculates matrix damps(power,dist,state) of dispersion damping
c  functions.
c***********************************************************************
      implicit real*8 (a-h,o-z)
      common /discom/alph0,alph1,bet0,bet1
      common /dampp/damps(10,3,2)
      common /h2com/rm,r0,rmt
      dimension r(3), a(10),b(10)
      do 1 i=6,10,2
      a(i)=alph0/float(i)**alph1
      b(i)=bet0*exp(-bet1*float(i))
1     continue
      do 2 is=1,2
      xm=rm
      if(is.eq.2)xm=rmt
      do 2 id=1,3
      x=2.0d0*r(id)/(xm+2.5d0*r0)
      xs=x*x
      do 2 n=6,10,2
      pol=a(n)*x+b(n)*xs
2     damps(n,id,is)=(1.0d0-exp(-pol))**n
      return
      end
      function damdamp(r,n,rm,r0)
c***********************************************************************
c     calculates 1st derivative of damping function.
c***********************************************************************
      implicit real*8 (a-h,o-z)
      common /discom/alph0,alph1,bet0,bet1
      a=alph0/float(n)**alph1
      b=bet0*exp(-bet1*float(n))
      denom=rm+2.5d0*r0
      x=2.0d0*r/denom
      pol=a*x+b*x**2
      dpol=a+2.d0*b*x
      damdamp=(float(2*n)/denom)*dpol*exp(-pol)*(1.d0-exp(-pol))**(n-1)
      return
      end
      function g(r,ck0,ck1)
c***********************************************************************
c     calculates g function.
c***********************************************************************
      implicit real*8 (a-h,o-z)
      common /h2com/rm,r0,rmt
      g=1.d0+ck0*exp(-ck1*(r-rm))
      return
      end
      function gd(r,ck0,ck1)
c***********************************************************************
c     calculates 1st derivative of g function.
c***********************************************************************
      implicit real*8 (a-h,o-z)
      common /h2com/rm,r0,rmt
      gd=-ck1*ck0*exp(-ck1*(r-rm))
      return
      end
      function h(r,et)
c***********************************************************************
c     calculates h function.
c***********************************************************************
      implicit real*8 (a-h,o-z)
      h=(tanh(et*r))**6
      return
      end
      function hd(r,et)
c***********************************************************************
c     calculates 1st derivative of h function.
c***********************************************************************
      implicit real*8 (a-h,o-z)
      hd=6.0d0*tanh(et*r)**5*et/cosh(et*r)**2
      return
      end
      function h2hf(r)
c***********************************************************************
c     calculates extended-hartree-fock curve for ground singlet state
c     of h2 [a.j.c. varandas & j.d. silva, j. chem. soc. faraday ii
c     (submitted)].
c***********************************************************************
      implicit real*8 (a-h,o-z)
      common /diacom/d,a1,a2,a3,gamma
      common /h2com/rm,r0,rmt
      dr=r-rm
      h2hf=-d*(1.0d0+a1*dr+a2*dr**2+a3*dr**3)*exp(-gamma*dr)
      return
      end
      function h2hfd(r)
c***********************************************************************
c     calculates 1st derivative of 2-body extended-hartree-fock curve
c     for the ground-singlet state of h2.
c***********************************************************************
      implicit real*8 (a-h,o-z)
      common /diacom/d,a1,a2,a3,gamma
      common /h2com/rm,r0,rmt
      dr=r-rm
      h2hfd=-gamma*h2hf(r)-d*(a1+2*a2*dr+3.0d0*a3*dr**2)*exp(-gamma*dr)
      return
      end
      function switch(r1,r2,r3)
c***********************************************************************
c     calculates switching term for the hartree-fock component of the diatomic
c     triplet state function.  the notation of thompson et al. (j. chem. phys.,
c     82,5597,1985; and references therein) is used throughout.
c***********************************************************************
      implicit real*8 (a-h,o-z)
      common/acdcom/alph2,az2,ce0,cd1
      common /picom/ pi,sqrt3
      common /crdcom/ r12,r22,r32,qcoord,rho,s,cphi3,s2,s3
      common /swfcom/sf
      if(s.ne.0.0d0)go to 1
      switch=exp(-az2*qcoord**2)
      go to 11
    1 switch=exp(-az2*(1.d0+s3*cphi3)*qcoord**2)
   11 sf=switch
      return
      end
      function trip(r,rsquar,al0,al1,al2,al3)
c***********************************************************************
c     calculates effective diatomic triplet state curve.
c***********************************************************************
      implicit real*8 (a-h,o-z)
      trip=(al0+al1*r+al2*rsquar)*exp(-al3*r)
      return
      end
      function vleps(r1,r2,r3)
c***********************************************************************
c     calculates 3-body extended-hartree-fock energy defined by a
c     leps-type function.
c***********************************************************************
      implicit real*8 (a-h,o-z)
      common /crdcom/ r12,r22,r32,qcoord,rho,s,cphi3,s2,s3
      common /h2com/rm,r0,rmt
      common /trncom/al0,al1,al2,al3
      common /lepcom/ vleps2,dlep(3),dlep2(3)
      common /corr/corr2s(3,2)
c   vleps2 is leps upper surface.
c   dlep(3) are leps lower surface derivatives (not yet coded)
c   dlep2(3) "   "   upper    "       "          "    "    "
      sng1=h2hf(r1)
      sng2=h2hf(r2)
      sng3=h2hf(r3)
      if(1)11,1,11
   11 f=switch(r1,r2,r3)
      call acct(r1,hft1,cet1,at1,1)
      cesng1=corr2s(1,1)
      t1=trip(r1,r12,al0,al1,al2,al3)*f+(1.d0-f)*(at1-cesng1)
      call acct(r2,hft2,cet2,at2,2)
      cesng2=corr2s(2,1)
      t2=trip(r2,r22,al0,al1,al2,al3)*f+(1.d0-f)*(at2-cesng2)
      call acct(r3,hft3,cet3,at3,3)
      cesng3=corr2s(3,1)
      t3=trip(r3,r32,al0,al1,al2,al3)*f+(1.d0-f)*(at3-cesng3)
      go to 2
    1 t1=trip(r1,r12,al0,al1,al2,al3)
      t2=trip(r2,r22,al0,al1,al2,al3)
      t3=trip(r3,r32,al0,al1,al2,al3)
    2 q1=0.5d0*(sng1+t1)
      q2=0.5d0*(sng2+t2)
      q3=0.5d0*(sng3+t3)
      ex1=0.5d0*(sng1-t1)
      ex2=0.5d0*(sng2-t2)
      ex3=0.5d0*(sng3-t3)
      f1=(ex1-ex2)**2
      f2=(ex2-ex3)**2
      f3=(ex3-ex1)**2
      qsum=q1+q2+q3
      exch=sqrt(0.5d0*(f1+f2+f3))
      vleps=qsum-exch
      vleps2=qsum+exch
      return
      end
      subroutine acct(r,hft,cet,t,i)
c***********************************************************************
c     computes the hface (hartree-fock-approximate correlation energy)
c     potential for the lowest triplet state of h2 [a.j.c. varandas &
c     j. brandao mol. phys.,45,1982,857].
c***********************************************************************
      implicit real*8 (a-h,o-z)
      common /h2com/rm,r0,rmt
      common /corr/corr2s(3,2)
      data a,b1,b2,b3,b4/0.448467d0,-0.056687d0,0.246831d0,-0.018419d0,
     *                   0.000598d0/
      r2=r*r
      r3=r2*r
      r4=r3*r
      hft=a*exp(-(b1*r+b2*r2+b3*r3+b4*r4))/r
      cet=corr2s(i,2)
      t=hft+cet
      return
      end
      subroutine coordx(r1,r2,r3)
c**********************************************************************
c   calculates d3h symmetry coordinates
c**********************************************************************
      implicit real*8 (a-h,o-z)
      common /crdcom/ r12,r22,r32,qcoord,rho,s,cphi3,s2,s3
      common /picom/pi,sqrt3
      r12=r1**2
      r22=r2**2
      r32=r3**2
      qcoord=r12+r22+r32
      rho=sqrt(qcoord/3.d0)
      gamma=2.d0*r12-r22-r32
      beta=sqrt3*(r22-r32)
      s2=(beta**2+gamma**2)/qcoord**2
      s=sqrt(s2)
      s3=s2*s
      if(s.ne.0.d0)go to 11
c   if(s.eq.0.0d0) arg,theta,phi3,cphi3 are not computed.
      go to 45
   11 arg=gamma/qcoord/s
      if(arg.le.-1.d0)arg=-1.d0
      if(arg.ge.1.d0)arg=1.d0
      theta=acos(arg)
      phi3=3.d0*(pi-theta)
      cphi3=cos(phi3)
   45 return
      end
      subroutine corr3(r1,r2,r3,rm,r0,ce3)
c***********************************************************************
c     calculates 3-body correlation energy.
c***********************************************************************
      implicit real*8 (a-h,o-z)
      dimension t1(10),t2(10),t3(10),co3(10)
      common /coecom/ck0(10),ck1(10)
      common /noncom/chh(10),chhe(10),chh2(10)
      common /dampp/damps(10,3,2)
c  define two statement functions
      g(x,c0,c1)=1.d0+c0*exp(-c1*(x-rm))
      h(x,c1)=(tanh(x*c1))**6
      do 1 i=6,10,2
         g1=g(r1,ck0(i),ck1(i))
         g2=g(r2,ck0(i),ck1(i))
         g3=g(r3,ck0(i),ck1(i))
         h1=h(r1,ck1(i))
         h2=h(r2,ck1(i))
         h3=h(r3,ck1(i))
         t1(i)=chh(i)*(1.0d0-0.5d0*(g2*h3+g3*h2))*damps(i,1,1)/r1**i
         t2(i)=chh(i)*(1.0d0-0.5d0*(g3*h1+g1*h3))*damps(i,2,1)/r2**i
         t3(i)=chh(i)*(1.0d0-0.5d0*(g1*h2+g2*h1))*damps(i,3,1)/r3**i
         co3(i)=t1(i)+t2(i)+t3(i)
    1 continue
      ce3=co3(6)+co3(8)+co3(10)
      return
      end
      subroutine corr3d(r1,r2,r3,rm,r0,dc3r1,dc3r2,dc3r3)
c***********************************************************************
c     calculates first-order partial derivatives of 3-body correlation
c     energy with respect to interatomic distances.
c***********************************************************************
      implicit real*8 (a-h,o-z)
      dimension d1(10),d2(10),d3(10)
      common /coecom/ck0(10),ck1(10)
      common /noncom/chh(10),chhe(10),chh2(10)
      dc3r1=0.d0
      dc3r2=0.d0
      dc3r3=0.d0
      do 1 i=6,10,2
         g1=g(r1,ck0(i),ck1(i))
         g2=g(r2,ck0(i),ck1(i))
         g3=g(r3,ck0(i),ck1(i))
         gd1=gd(r1,ck0(i),ck1(i))
         gd2=gd(r2,ck0(i),ck1(i))
         gd3=gd(r3,ck0(i),ck1(i))
         h1=h(r1,ck1(i))
         h2=h(r2,ck1(i))
         h3=h(r3,ck1(i))
         hd1=hd(r1,ck1(i))
         hd2=hd(r2,ck1(i))
         hd3=hd(r3,ck1(i))
         dd1=damdamp(r1,i,rm,r0)
         dd2=damdamp(r2,i,rm,r0)
         dd3=damdamp(r3,i,rm,r0)
         gh23=g2*h3+g3*h2
         gh31=g3*h1+g1*h3
         gh12=g1*h2+g2*h1
         d1(i)=chh(i)*dd1*(1.d0-0.5d0*gh23)/r1**i-
     1         float(i)*chh(i)*damdamp(r1,i,rm,r0)*(1.d0-0.5d0*gh23
     2         )/r1**(i+1)+
     3     chh(i)*damdamp(r2,i,rm,r0)*(-0.5d0*(g3*hd1+gd1*h3))/r2**i+
     4         chh(i)*damdamp(r3,i,rm,r0)*(-0.5d0*(gd1*h2+g2*hd1))/r3**i
         d2(i)=chh(i)*dd2*(1.d0-0.5d0*gh31)/r2**i-
     1         float(i)*chh(i)*damdamp(r2,i,rm,r0)*(1.d0-0.5d0*gh31
     2         )/r2**(i+1)+
     3     chh(i)*damdamp(r1,i,rm,r0)*(-0.5d0*(gd2*h3+g3*hd2))/r1**i+
     4         chh(i)*damdamp(r3,i,rm,r0)*(-0.5d0*(g1*hd2+h1*gd2))/r3**i
         d3(i)=chh(i)*dd3*(1.d0-0.5d0*gh12)/r3**i-
     1         float(i)*chh(i)*damdamp(r3,i,rm,r0)*(1.d0-0.5d0*gh12
     2         )/r3**(i+1)+
     3     chh(i)*damdamp(r1,i,rm,r0)*(-0.5d0*(g2*hd3+gd3*h2))/r1**i+
     4         chh(i)*damdamp(r2,i,rm,r0)*(-0.5d0*(gd3*h1+g1*hd3))/r2**i
         dc3r1=dc3r1+d1(i)
         dc3r2=dc3r2+d2(i)
         dc3r3=dc3r3+d3(i)
    1 continue
      return
      end
      block data h3data
c***********************************************************************
c     data for the dmbe h3 surface.
c***********************************************************************
      implicit real*8 (a-h,o-z)
      common /acdcom/alph2,az2,ce0,cd1
      common /barcom/rsp,vspkc
      common /diacom/d,a1,a2,a3,gamma
      common /discom/alph0,alph1,bet0,bet1
      common /fixcom/alpha5,beta1,beta2,beta3
      common /h2com/rm,r0,rmt
      common /noncom/chh(10),chhe(10),chh2(10)
      common /coecom/ck0(10),ck1(10)
      common /picom/pi,sqrt3
      common /trncom/al0,al1,al2,al3
      common /xcocom/x(15)
      data alph0,alph1,bet0,bet1/25.9528d0,1.1868d0,15.7381d0,0.09729d0/
      data alpha5,beta1,beta2,beta3/ 8.2433033d-03,
     2 0.53302897d+00,0.39156612d-01,0.69996945d+00/
      data alph2/4.735364d-1/
      data al0,al1,al2,al3/
     2 0.45024472d+01,-0.62467617d+01,0.40966542d+01,0.21813012d+01/
      data az2/4.453649d-4/
      data ce0,cd1/6.333404d-03,-1.726839d-03/
      data chh(6),chh(8),chh(10)/6.499027d0,124.3991d0,3285.8d0/
      data chhe(6),chhe(8),chhe(10)/2.82d0,41.8d0,865.6d0/
      data chh2(6),chh2(8),chh2(10)/8.813d0,163.87d0,4023.2d0/
      data ck0,ck1/       0.000000d+00,      0.000000d+00,
     *      0.000000d+00,      0.000000d+00,      0.000000d+00,
     *     -1.372843d-01,      0.000000d+00,     -1.638459d-01,
     *      0.000000d+00,     -1.973814d-01,      0.000000d+00,
     *      0.000000d+00,      0.000000d+00,      0.000000d+00,
     *      0.000000d+00,      1.011204d+00,      0.000000d+00,
     *      9.988099d-01,      0.000000d+00,      9.399411d-01/
      data d,a1,a2,a3,gamma/0.15796326d0,2.1977034d0,1.2932502d0,
     *                      0.64375666d0,2.1835071d0/
      data pi/3.14159265358979d0/
      data rm,r0,rmt/1.401d0,6.928203d0,7.82d0/
      data rsp,vspkc/1.757d0,9.65d0/
      data sqrt3/1.73205080756887d0/
      data x/-0.9286579d-02,  0.2811592d-03, -0.4665659d-05,
     2        0.2069800d-07,  0.2903613d+02, -0.2934824d+01,
     3        0.7181886d+00, -0.3753218d+00, -0.1114538d+01,
     4        0.2134221d+01, -0.4164343d+00,  0.2022584d+00,
     5       -0.4662687d-01, -0.4818623d+02,  0.2988468d+00/
      end
