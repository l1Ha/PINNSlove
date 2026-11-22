!     File created at Fri Jun  5 21:58:55 PDT 2020
!     Original source code: direct.f
module mod_direct 
      private
      public :: direct
      contains
        subroutine direct (rho,cvi,cro,ilev,jlev,klev,h,v,ia)
           use mod_common
           use mod_potenl, only: potenl 
      use mod_pes_multi, only: NPES, active_npes, potenl_matrix
          use mod_arrang , only: arrang 
          use mod_qrot , only: qrot
          use mod_qvib , only: qvib
          use mod_prot , only: prot
          use mod_pvibv , only: pvibv
        implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine fills in the direct matrix elements of h
!     within arrangement ia at hyperradius rho. Last modified 9/6/99.
!     -----------------------------------------------------------------
!
!     common blocks
!
      ! COMMON blocks replaced by `mod_common`
!
!     input arrays
!
      dimension cvi(nvi,n),cro(3,n)
      dimension ilev(n),jlev(n),klev(n)
      dimension h(nchan,nchan),v(mro,mvi)
!
!     local arrays
!
!     dimension wro(mro),xro(mro)
      allocatable wro(:),xro(:)

!     dimension wvi(mvi),xvi(mvi)
      allocatable wvi(:),xvi(:)

!     dimension pro(0:jmax,0:kmax,mro)
      allocatable pro(:,:,:)

!     dimension pvi(mvi,n)
      allocatable pvi(:,:)

   !     potential matrices for multi-surface
      allocatable Vcur(:,:), Vref(:,:), pqq(:,:)
      integer pp,qq,iaf,ibf
      integer ch,i,base,pes,npes_loc,nw,ilev_val
      double precision vab_pq, sabs

!     dimension hro(0:jmax,0:jmax,0:kmax)
      allocatable hro(:,:,:,:,:)

      allocate (wro(mro),xro(mro))
      allocate (wvi(mvi),xvi(mvi))
      allocate (pro(0:jmax,0:kmax,mro))
      allocate (pvi(mvi,n))
      npes_max = max(1,NPES)
      npes_act = active_npes(rho)
      allocate (hro(0:jmax,0:jmax,0:kmax,1:npes_max,1:npes_max))
      ! hro(ja,jb,k,p,q) stores rotational integrals weighted by potential element (p,q)

!
!     arrangement indices
!
      call arrang (ilev,n,jpar,ia,nla,nha,na)
      if (na .lt. 1) return
!
!     quadrature rules
!
      call qrot (mro,wro,xro)
      tmin = 0.d0
      tmax = asin(min(1.d0,smax/rho))
      call qvib (tmin,tmax,mvi,wvi,xvi)
!
!     rotational basis functions
!
      do kro = 1,mro
         cosa = xro(kro)
         call prot (pro(0,0,kro),jmax,kmax,cosa)
      enddo
!
!     vibrational basis functions
!
      call pvibv (tmin,xvi,tmax,cvi(1,nla),nvi,na,pvi(1,nla),mvi,0)

   !     allocate potential scratch matrices once and reuse
      allocate(Vref(1:npes_max,1:npes_max))
      allocate(pqq(1:npes_max,1:npes_max))

   !     loop over quadrature points

      do kvi = 1,mvi
         ta = xvi(kvi)
         ra = rho*cos(ta)
         sa = rho*sin(ta)
         ! compute reference potential matrix at large ra for this ta
         call potenl_matrix(100.d0,sa,0.d0,Vref,ia)
          do pp = 1, npes_act
             do qq = 1, npes_act
                Vref(pp,qq) = rmlmda * Vref(pp,qq)
             enddo
          enddo
!
!        rotational integrals
!
         do k = kmin,kmax
            do jb = k,jmax
               do ja = k,jb
                  do pp = 1, max(1,NPES)
                     do qq = 1, max(1,NPES)
                        hro(ja,jb,k,pp,qq) = 0.d0
                     enddo
                  enddo
               enddo
            enddo
         enddo
         do kro = 1,mro
            cosa = xro(kro)
            ! potential matrix at this quadrature point
             allocate(Vcur(1:npes_max,1:npes_max))
             call potenl_matrix(ra,sa,cosa,Vcur,ia)
             do pp = 1, npes_act
                do qq = 1, npes_act
                   Vcur(pp,qq) = rmlmda * Vcur(pp,qq)
                   pqq(pp,qq) = Vcur(pp,qq) - Vref(pp,qq)
                enddo
            enddo
            sab = wro(kro)*wvi(kvi)
            nj = jmax+1
            do k = kmin,kmax
               nk = nj-k
                do pp = 1, npes_act
                   do qq = 1, npes_act
                      vab_pq = sab * pqq(pp,qq)
                      call dsyr ('u',nk,vab_pq,pro(k,k,kro),1,hro(k,k,k,pp,qq),nj)
                   enddo
               enddo
            enddo
            deallocate(Vcur)
         enddo
         do k = kmin,kmax
            do jb = k,jmax
               do ja = k,jb
                   do pp = 1, npes_act
                      do qq = 1, npes_act
                         hro(jb,ja,k,pp,qq) = hro(ja,jb,k,pp,qq)
                      enddo
                   enddo
               enddo
            enddo
         enddo
!
!        vibrational integrals and Coriolis coupling
!
         cent = wvi(kvi)/ra**2
         do lb = nla,nha
            jb = jlev(lb)
            k = klev(lb)
            cmb = cent*cro(1,lb)*pvi(kvi,lb)
            cpb = cent*cro(3,lb)*pvi(kvi,lb)
            do la = nla,lb
               ja = jlev(la)
               ka = klev(la)
               if (ka .eq. k) then
                   do pp = 1, npes_act
                      do qq = 1, npes_act
                         iaf = (la-1)*npes_max + pp
                         ibf = (lb-1)*npes_max + qq
                         h(iaf,ibf) = h(iaf,ibf) + pvi(kvi,la) * hro(ja,jb,k,pp,qq) * pvi(kvi,lb)
                      enddo
                   enddo
                else if (ja .eq. jb) then
                   if (ka .eq. k-1) then
                      do pp = 1, npes_act
                         iaf = (la-1)*npes_max + pp
                         ibf = (lb-1)*npes_max + pp
                         h(iaf,ibf) = h(iaf,ibf) + pvi(kvi,la)*cmb
                      enddo
                   else if (ka .eq. k+1) then
                      do pp = 1, npes_act
                         iaf = (la-1)*npes_max + pp
                         ibf = (lb-1)*npes_max + pp
                         h(iaf,ibf) = h(iaf,ibf) + pvi(kvi,la)*cpb
                      enddo
                  endif
               endif
            enddo
         enddo
      enddo
      deallocate(Vref)
      if (allocated(pqq)) deallocate(pqq)
      ! symmetrise full h matrix across nchan
      do lb = 1,nchan
         do la = 1,lb-1
            h(lb,la) = h(la,lb)
         enddo
      enddo
      return
      end
   end module mod_direct 
